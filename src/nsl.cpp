#include <vector>
#include <string>

#include "nsl.hpp"
#include<thread>
#include "hapmap.hpp"
#include "utils.hpp"
using namespace std;

NSL::NSL () {
}

NSL::~NSL () {
	//output_filename = "";
	delete[] nsl1;
	delete[] nsl0;
}

void NSL::cli_prepare(CLI::App * app) {
	this->subapp = app->add_subcommand("nsl", "Calculate nSL.");

    CLI::Option * opt;
	opt = subapp->add_option("-p, --hap", input_filename_hap, "A hapfile with one row per haplotype, \n"
                                                                "and one column per variant.\n"
                                                                "Variants should be coded 0/1")->capture_default_str();
	//opt->required();
	opt->check(CLI::ExistingFile);


    opt = subapp->add_option("--vcf", input_filename_vcf, "A VCF file with one row per haplotype, \n"
                                                                "and one column per variant.\n"
                                                                "Variants should be coded 0/1")->capture_default_str();
    opt->check(CLI::ExistingFile);
    
    


	numThread = std::thread::hardware_concurrency(); //may return 0 when not able to detect
	opt = subapp->add_option("-t, --threads", numThread, "The number of threads to spawn during the calculation. Partitions loci across threads. (by default uses all available cores)")->capture_default_str();
	
    // CLI::Option * in_option2 = subapp->add_option("-m, --map", input_filename_map, "Input MAP file");
	// in_option2->required();
	// in_option2->check(CLI::ExistingFile);
        
    CLI::Option_group * group=subapp->add_option_group("Set Locus/All");
    group->add_flag("--all", calc_all, "Calculate nsl for all loci.");
    CLI::Option* optt = group->add_option<unsigned int>("-l, --loc", locus, "Locus")->capture_default_str();
    optt->check( [](const std::string &str) {
                if(stoi(str) < 0)
                    return std::string("Locus must be non-negative.");
                else
                    return std::string();
                });
    group->require_option(1); 

	opt = subapp->add_option<double>("--cutoff", cutoff, "The EHH decay cutoff")->capture_default_str();
	opt->check(CLI::Range(0.0,1.0));
    opt->option_text("<float> [0.05]");


    opt = subapp->add_option<double>("--maf", min_maf, "If a site has a MAF below this value, the program \n"
                                                        "will not use it as a core snp.")->capture_default_str();
	opt->check(CLI::Range(0.0,1.0));
    opt->option_text("<float> [0.05]");

	opt = subapp->add_option<string>("-o, --out", output_filename, "The basename for all output files reporting stats. \n"
                                        "Formatted: <locusID> <physicalPos> <1 freq> <sl1> <sl0> <unstandardized nSL>")->capture_default_str();;
    //worry about it later!
    
    opt = subapp->add_option<unsigned int>("--max-extend-nsl", max_extend_nsl, "The maximum distance an nSL haplotype \n" 
                                                                                "is allowed to extend from the core. \n"
                                                                                "Set = 0 for no restriction. ")->capture_default_str();;
    

    //FLAGS
    subapp->add_flag("--alt", alt, "Set this flag to calculate homozygosity based on\n" 
                                    "the sum of the squared haplotype frequencies in the \n"
                                    "observed data instead of using binomial coefficients.");

    

	app->get_formatter()->column_width(10);
}

void NSL::init() {  
    cout<<"v3.0.0"<<endl;
    cout<<"Calculating nSL"<<endl;
    cout<<"Input filename vcf: "<<input_filename_vcf<<endl;
    cout<<"Input filename: "<<input_filename_hap<<endl;
    cout<<"Map filename: "<<input_filename_map<<endl;
    cout<<"Output file: "<<output_filename<<endl;
    cout<<"Threads: "<<numThread<<endl;  
    cout<<"Scale parameter: "<<20000<<endl;  
    cout<<"Max gap parameter: "<<max_extend_nsl<<endl;  
    cout<<"EHH cutoff value: "<<cutoff<<endl;  
    cout<<"Phased: "<<"Yes"<<endl;  
    cout<<"Alt flag set: "<<"No"<<endl;  

    useVCF=true;
    if(useVCF){
        cout<<"Loading "<<input_filename_vcf<<endl;
    	hm.loadVCF(input_filename_vcf.c_str(), min_maf); //populate the matrix
    }else{
        cout<<"Loading "<<input_filename_hap<<endl;
    	hm.loadHapMap(input_filename_hap.c_str(), input_filename_map.c_str(), min_maf); //populate the matrix
    }
    

    cout<<"Loading "<<input_filename_hap<<endl;
	//hm.loadHap(input_filename_hap.c_str(), min_maf, this->all_positions); //populate the hapmap
	
    
    //hm.loadHapMap(input_filename_hap.c_str(), input_filename_map.c_str(), min_maf, this->all_positions, mentries); //populate the hapmap
	
    //hm.loadHapMapVCF(input_filename_hap.c_str(), input_filename_map.c_str(), min_maf, this->all_positions, mentries); //populate the hapmap
	
    this->numHaps = hm.numHaps();
	this->numSnps = hm.numSnps();

	nsl1 = new long[numSnps];
	nsl0 = new long[numSnps];

    //calc_nSL(locus);
}

void NSL::exec() {
	this->init(); //initialize
    if(calc_all){
        calc_nSL_all();
    }else{
        calc_nSL(locus);
    }
}

void NSL::calc_nSL_upstream(int locus, map<int, vector<int> > & m){
        if(locus >= numHaps){
            return;
        }
        //if group_count[location] == 1 already unique
        //stop when totgc = N
        //nsl0[locus] = 0;
        //nsl1[locus] = 0;
        
        int n_c0= 0;
        int n_c1=0;
        int n_c0_squared_minus = 0;
        int n_c1_squared_minus = 0;

        int group_count[numHaps];
        int group_id[numHaps];
        bool isDerived[numHaps]; 

        //will be vectorized
        for(int i = 0; i<numHaps; i++){
            group_count[i] = 0;
            group_id[i] = 0;
            isDerived[i] = false;
        }

        int totgc=0;
        vector<unsigned int> v = hm.all_positions[locus];

        if(v.size()==0){
            n_c0 = numHaps;

            group_count[0] = numHaps;
            totgc+=1;

        }else if (v.size()==numHaps){ // all set
            group_count[0] = numHaps;
            totgc+=1;
            n_c1 = numHaps;
            
            for (int set_bit_pos : v){
                isDerived[set_bit_pos] = true;
            }

        }else{
            group_count[1] = v.size();
            group_count[0] = numHaps - v.size();
            n_c1 = v.size();
            n_c0 = numHaps - v.size();

            for (int set_bit_pos : v){
                group_id[set_bit_pos] = 1;
                isDerived[set_bit_pos] = true;

            }
            
            totgc+=2;

        }

        
        n_c1_squared_minus =  n_c1* n_c1 -  n_c1;
        n_c0_squared_minus =  n_c0* n_c0 -  n_c0;
        

        //print_a<bool>(isDerived, "derived");
        //print_a(group_count, "GC");
        //print_a(group_id, "GID");

        for ( int i = locus+1; i<numSnps; i++ ){
            if(i-locus > max_extend_nsl){
                break;
            }
            for (int set_bit_pos : hm.all_positions[i])
            {
                int old_group_id = group_id[set_bit_pos];
                m[old_group_id].push_back(set_bit_pos);
            }

            for (const auto &ele : m)
            {
                int old_group_id = ele.first;
                int newgroup_size = ele.second.size();

                if (group_count[old_group_id] == newgroup_size || newgroup_size == 0)
                {
                    continue;
                }

                for (int v : ele.second)
                {
                    group_id[v] = totgc;
                }

                int del_update = 0;
                del_update=  (group_count[old_group_id]-newgroup_size) * newgroup_size * (i-locus);
                
                group_count[old_group_id] -= newgroup_size;
                
                bool isDerivedGroup =  isDerived[ele.second[0]];
                
                if(isDerivedGroup)// just check first element if it is derived or not, 
                {
                    nsl1[locus] += del_update;
                }else{
                    nsl0[locus] += del_update;
                }
                group_count[totgc] += newgroup_size;
                totgc += 1;
            }
            m.clear();
        

            if(n_c1_squared_minus==0){
                nsl1[locus] = 0;
            }
            if(n_c0_squared_minus==0){
                nsl0[locus] = 0;
            }
           
            
            if(i==numSnps-1){
                set<int> derivedGroups;
                for(int i: hm.all_positions[locus]){
                    derivedGroups.insert(group_id[i]);
                }
                //for all groups having a non-1 group count do nc2
                for(int k = 0; k<totgc; k++){
                    int count = group_count[k];
                    if(count > 1){ //TRY TO OPT
                        if(derivedGroups.count(k)){
                           nsl1[locus]+=0.5*count*(count - 1)*(i-locus+1);
                        }else{
                           nsl0[locus]+=0.5*count*(count - 1)*(i-locus+1);
                        }
                    }
                }
            }
            // if(true)
              //  std::cout<<"Iter "<<i<<": nsl[1,0]["<<locus<<"]="<<nsl1[locus]*2.0/n_c1_squared_minus<<","<<nsl0[locus]*2.0/n_c0_squared_minus<<endl;
            //std::cout<<"Iter "<<i<<": nsl[1,0]["<<locus<<"]="<<nsl1[locus]<<","<<nsl0[locus]<<endl;


            if(totgc==numHaps){
                return;
            }
        }
    }


void NSL::calc_nSL_downstream(int locus, map<int, vector<int> > & m){
    
        //if group_count[location] == 1 already unique
        //stop when totgc = N
        
        if(locus<=0) return;

        //nsl0[locus] = 0;
        //nsl1[locus] = 0;
        
        int n_c0= 0;
        int n_c1=0;
        int n_c0_squared_minus = 0;
        int n_c1_squared_minus = 0;

        int group_count[numHaps];
        int group_id[numHaps];
        bool isDerived[numHaps]; 

        //will be vectorized
        for(int i = 0; i<numHaps; i++){
            group_count[i] = 0;
            group_id[i] = 0;
            isDerived[i] = false;
        }

        int totgc=0;
        vector<unsigned int> v = hm.all_positions[locus];

        if(v.size()==0){
            n_c0 = numHaps;

            group_count[0] = numHaps;
            totgc+=1;

        }else if (v.size()==numHaps){ // all set
            group_count[0] = numHaps;
            totgc+=1;
            n_c1 = numHaps;
            
            for (int set_bit_pos : v){
                isDerived[set_bit_pos] = true;
            }

        }else{
            group_count[1] = v.size();
            group_count[0] = numHaps - v.size();
            n_c1 = v.size();
            n_c0 = numHaps - v.size();

            for (int set_bit_pos : v){
                group_id[set_bit_pos] = 1;
                isDerived[set_bit_pos] = true;

            }
            
            totgc+=2;

        }

        
        n_c1_squared_minus =  n_c1* n_c1 -  n_c1;
        n_c0_squared_minus =  n_c0* n_c0 -  n_c0;
        

        // print_a<bool>(isDerived, "derived");
        // print_a(group_count, "GC");
        // print_a(group_id, "GID");

        for ( int i = locus-1; i>=0; i-- ){
            if(locus-i > max_extend_nsl){
                break;
            }

            for (int set_bit_pos : hm.all_positions[i])
            {
                int old_group_id = group_id[set_bit_pos];
                m[old_group_id].push_back(set_bit_pos);
            }

            for (const auto &ele : m)
            {
                int old_group_id = ele.first;
                int newgroup_size = ele.second.size();

                if (group_count[old_group_id] == newgroup_size || newgroup_size == 0)
                {
                    continue;
                }

                for (int v : ele.second)
                {
                    group_id[v] = totgc;
                }

                int del_update = 0;
                del_update=  (group_count[old_group_id]-newgroup_size) * newgroup_size * (locus-i);
                
                group_count[old_group_id] -= newgroup_size;
                
                bool isDerivedGroup =  isDerived[ele.second[0]];
                
                if(isDerivedGroup)// just check first element if it is derived or not, 
                {
                    nsl1[locus] += del_update;
                }else{
                    nsl0[locus] += del_update;
                }
                group_count[totgc] += newgroup_size;
                totgc += 1;
            }
            m.clear();
        

            if(n_c1_squared_minus==0){
                nsl1[locus] = 0;
            }
            if(n_c0_squared_minus==0){
                nsl0[locus] = 0;
            }
           
            
            if(i==0){
                set<int> derivedGroups;
                for(int i: hm.all_positions[locus]){
                    derivedGroups.insert(group_id[i]);
                }
                //for all groups having a non-1 group count do nc2
                for(int k = 0; k<totgc; k++){
                    int count = group_count[k];
                    if(count > 1){ //TRY TO OPT
                        if(derivedGroups.count(k)){
                           nsl1[locus]+=0.5*count*(count - 1)*(locus-i+1);
                        }else{
                           nsl0[locus]+=0.5*count*(count - 1)*(locus-i+1);
                        }
                    }
                }
            }
            
            
            //std::cout<<"Iter "<<i<<": nsl[1,0]["<<locus<<"]="<<nsl1[locus]<<","<<nsl0[locus]<<endl;


            if(totgc==numHaps){
                return;
            }
        }
    }

void NSL::calc_nSL(int locus){
    if(!calc_all)
        out_fp = fopen(("outbin."+to_string(locus)+".out").c_str(),"w");
	map<int, vector<int> > m;

	nsl0[locus] = 0;
	nsl1[locus] = 0;

	calc_nSL_upstream(locus, m);
	calc_nSL_downstream(locus, m);

	float nsl_1 = 0.0;
	float nsl_0 = 0.0;
	
	int numDerivedAtLocus = hm.all_positions[locus].size();
	int numAncestralAtLocus = numHaps - numDerivedAtLocus;2.0 / (numDerivedAtLocus*numDerivedAtLocus - numDerivedAtLocus);
	int numDerivedPairs = (numDerivedAtLocus*numDerivedAtLocus - numDerivedAtLocus) / 2;
	int numAncestralPairs = (numAncestralAtLocus*numAncestralAtLocus - numAncestralAtLocus) / 2;


	if(numDerivedAtLocus>1){
		nsl_1 = nsl1[locus] * 2.0 / (numDerivedAtLocus*numDerivedAtLocus - numDerivedAtLocus)  ;
		if(locus!=0 and locus!=numSnps-1){
			nsl_1 -= 1;
		}
	}
	if(numAncestralAtLocus>1){
		nsl_0 = nsl0[locus] * 2.0 / (numAncestralAtLocus*numAncestralAtLocus - numAncestralAtLocus)  ;
		if(locus!=0 and locus!=numSnps-1){
			nsl_0 -= 1;
		}
	}

	//std::cout<<" nsl1["<<mentries[locus].locId<<"]="<<nsl_1<<",nsl0["<<mentries[locus].locId<<"]="<<nsl_0<<endl;

    //< locusID > < physicalPos > < ’1 ’ freq > <sl1 > <sl0 > < unstandardized nSL >
	//std::cout<<<< <<nsl_1<<",nsl0["<<mentries[locus].locId<<"]="<<nsl_0<<endl;
    if(nsl_1==0){nsl_1+=1;};
    if(nsl_0==0){nsl_0+=1;};

    if(!calc_all)
        fprintf(out_fp, "%d %d %f %f %f %f\n",hm.mentries[locus].locId, hm.mentries[locus].phyPos, hm.all_positions[locus].size() * 1.0/this->numHaps, nsl_1, nsl_0, log10(nsl_1/nsl_0)); 
}

void NSL::calc_nSL_all(){
    out_fp = fopen("outbin.nsl.out","w");
    for(int i = 0; i<numSnps; i++){
        calc_nSL(i);
    }
    fclose(out_fp);
	//return;
}