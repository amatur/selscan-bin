#include "ehh.hpp"
#include "benchmark.hpp"
#include "hapmap.hpp"
#include "utils.hpp"
#include <omp.h>

EHH::EHH(){
}
EHH::~EHH(){
    delete[] iHH0;
    delete[] iHH1;
    delete[] logg;
}


void EHH::loadMPHB(string input_mphb_file, int numHaps, int numSnps, map< int, map <int, priority_queue<pair<int, int>  > > > & mphbs){
    //vector< map<int, priority_queue<pair<int, int>  > > > mphbs(numHaps);
    //mphbs[0] initilize with numHaps 
    //for each hap index by left

    //struct hapblock

    //map< int, vector<priority_queue<pair<int, int>  > > >  mphbs;

    // for (int i = 0 ; i < numHaps; i++){
    //     //map<int, queue<pair<int, int> > > v;
    //     vector<priority_queue<pair<int, int> > > v;

    //     mphbs.push_back(v);
    //     // for (int j = 0 ; j < numSnps; j++){
    //     //     queue<pair<int, int> > q;
    //     //     mphbs[i].push_back(q);
    //     // }
    // }
     

    std::ifstream inputFile(input_mphb_file);
    if (!inputFile.is_open()) {
        std::cerr << "Error opening file." << std::endl;
        return;
    }

    std::string line;

    while (std::getline(inputFile, line)) {
        std::istringstream iss(line);
        std::string data;

        //struct hapblock hb; 
        int left, seqId, right, numSeq;;
        
        // Read four columns delimited by space
        (iss >> left) ;
        (iss >> right) ;
        (iss >> seqId) ;
        (iss >> numSeq) ;

        if(!mphbs.count(left)){
            priority_queue<pair<int, int> > q;
            map<int, priority_queue<pair<int, int> > > v;
            
            v[seqId] = q;
            mphbs[left] = v;
        }
        if(right-left>2){

            if(!mphbs[left].count(seqId)){
                priority_queue<pair<int, int> > q2;
                mphbs[left][seqId] = q2;
            }
            mphbs[left][seqId].push(make_pair(numSeq, right));
        }
            

        //cout<<seqId<<" "<<left<<" "<<right<<" "<<numSeq<<endl;
        //verify that it is larger than stack top
    }

    /*
    for (int i = 0 ; i < numHaps; i++){
        for (int j = 0 ; j < numSnps; j++){
            if(mphbs.count(j)){
                priority_queue<pair<int, int> > q = mphbs[j][i];
                if (q.size()!=0 ){
                    cout<<i << " "<< j << " "<< q.size()<<": ";
                    while (!q.empty()){
                        int numSeq = q.top().first;
                        int right = q.top().second;
                        q.pop();
                        cout<<numSeq<<" ("<<right<<")"<<" ";
                    }
                    
                    cout<<endl;
                }
            }
            
        }
    }
    */
}


void EHH::loadMPHB(string input_mphb_file, int numHaps, int numSnps, vector< map< int, priority_queue<pair<int, int>  > > > & mphbs){
    //vector< map<int, priority_queue<pair<int, int>  > > > mphbs(numHaps);
    //mphbs[0] initilize with numHaps 
    //for each hap index by left

    //struct hapblock

    for (int i = 0 ; i < numHaps; i++){
        //map<int, queue<pair<int, int> > > v;
        map<int, priority_queue<pair<int, int> > > v;

        mphbs.push_back(v);
        // for (int j = 0 ; j < numSnps; j++){
        //     queue<pair<int, int> > q;
        //     mphbs[i].push_back(q);
        // }
    }
     

    std::ifstream inputFile(input_mphb_file);
    if (!inputFile.is_open()) {
        std::cerr << "Error opening file." << std::endl;
        return;
    }

    std::string line;

    while (std::getline(inputFile, line)) {
        std::istringstream iss(line);
        std::string data;

        //struct hapblock hb; 
        int left, seqId, right, numSeq;;
        
        // Read four columns delimited by space
        (iss >> left) ;
        (iss >> right) ;
        (iss >> seqId) ;
        (iss >> numSeq) ;

        if(!mphbs[seqId].count(left)){
            priority_queue<pair<int, int> > q;
            mphbs[seqId][left] = q;
        }
        if(right-left>2)
            mphbs[seqId][left].push(make_pair(numSeq, right));

        //cout<<seqId<<" "<<left<<" "<<right<<" "<<numSeq<<endl;
        //verify that it is larger than stack top
    }

    /*
    for (int i = 0 ; i < numHaps; i++){
        for (int j = 0 ; j < numSnps; j++){
            if(mphbs[i].count(j)){
                priority_queue<pair<int, int> > q = mphbs[i][j];
                //if (q.size()!=0 )
                    cout<<i << " "<< j << " "<< q.size()<<": ";
                    while (!q.empty()){
                        int numSeq = q.top().first;
                        int right = q.top().second;
                        q.pop();
                        cout<<numSeq<<" ("<<right<<")"<<" ";
                    }
                    
                    cout<<endl;
            }
            
        }
    }
    */

}

void EHH::cli_prepare(CLI::App * app) {
	this->subapp = app->add_subcommand("ehh", "Calculate EHH.");

    CLI::Option_group * group_input = subapp->add_option_group("Input: ");
    CLI::Option* group_opt = group_input->add_option("--vcf", input_filename_vcf, "A VCF file with one row per haplotype, \n"
                                                                "and one column per variant.\n"
                                                                "Variants should be coded 0/1")->capture_default_str();
    group_opt->check(CLI::ExistingFile);
    
    
    group_opt = group_input->add_option("-p, --hap", input_filename_hap, "A hapfile with one row per haplotype, \n"
                                                                "and one column per variant.\n"
                                                                "Variants should be coded 0/1")->capture_default_str();
    group_opt->check(CLI::ExistingFile);
    //group_input->require_option(1); 

    useVCF = true;
    CLI::Option * opt;
	// opt = subapp->add_option("-p, --hap", input_filename_hap, "A hapfile with one row per haplotype, \n"
    //                                                             "and one column per variant.\n"
    //                                                             "Variants should be coded 0/1")->capture_default_str();
	// //opt->required();
	// opt->check(CLI::ExistingFile);


	numThread = std::thread::hardware_concurrency(); //may return 0 when not able to detect
	opt = subapp->add_option("-t, --threads", numThread, "The number of threads to spawn during the calculation. Partitions loci across threads. (by default uses all available cores)")->capture_default_str();
	
    // CLI::Option * in_option2 = subapp->add_option("-m, --map", input_filename_map, "Input MAP file");
	// in_option2->required();
	// in_option2->check(CLI::ExistingFile);
        
    CLI::Option_group * group=subapp->add_option_group("Set Locus/All");
    group->add_flag("--all", calc_all, "Calculate EHH for all loci.");
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
    //opt->option_text("<float> [0.05]");


    opt = subapp->add_option<double>("--maf", min_maf, "If a site has a MAF below this value, the program \n"
                                                        "will not use it as a core snp.")->capture_default_str();
	opt->check(CLI::Range(0.0,1.0));
    opt->option_text("<float> [0.05]");

	opt = subapp->add_option<string>("-o, --out", output_filename, "The basename for all output files reporting stats. \n"
                                        "Formatted: <locusID> <physicalPos> <1 freq> <sl1> <sl0> <unstandardized nSL>")->capture_default_str();;
    //worry about it later!
    
    opt = subapp->add_option<string>("--mphb", input_filename_mphb, "MPHB_FILE");;
    //worry about it later!

    opt = subapp->add_option<unsigned int>("--max-extend", max_extend, "The maximum distance an nSL haplotype \n" 
                                                                                "is allowed to extend from the core. \n"
                                                                                "Set = 0 for no restriction. ")->capture_default_str();;
    

    //FLAGS
    subapp->add_flag("--alt", alt, "Set this flag to calculate homozygosity based on\n" 
                                    "the sum of the squared haplotype frequencies in the \n"
                                    "observed data instead of using binomial coefficients.");

    subapp->add_flag("--om", openmp_enabled, "Set this flag to use openmp thread");



	app->get_formatter()->column_width(10);

    if(subapp->count("--vcf")>0){
        useVCF = true;
    }else{
        useVCF = false;
    }

    // if(subapp->count("--mphb")>0){
    //     haplo_enabled = true;
    // }else{
    //     haplo_enabled = false;
    // }
    
    

}

void EHH::init(){
    if(input_filename_mphb==""){
        haplo_enabled = false;
    }else{
        haplo_enabled = true;
    }
    cout<<"MPHB usage enabled is" <<haplo_enabled<<endl;

    useVCF=true;
    if(useVCF){
        cout<<"Loading "<<input_filename_vcf<<endl;
    	hm.loadVCF(input_filename_vcf.c_str(), min_maf); //populate the matrix
    }else{
        cout<<"Loading "<<input_filename_hap<<endl;
    	hm.loadHapMap(input_filename_hap.c_str(), input_filename_map.c_str(), min_maf); //populate the matrix
    }

    this->numHaps = hm.numHaps();
	this->numSnps = hm.numSnps();

    if(haplo_enabled)
        loadMPHB(input_filename_mphb, numHaps, numSnps, mphbs);

    iHH0 = new double[numSnps];
    iHH1 = new double[numSnps];

    logg = new string[numThread];
}

void EHH::calc_EHH(int locus){
    map<int, vector<int> > m;
    iHH0[locus] = 0;
    iHH1[locus] = 0;


    // if(hm.getMAF(locus) < min_maf || 1-hm.getMAF(locus) < min_maf){
    //     return;
    // }
    // ehh1[locus] = 0;
    // ehh0[locus] = 0;
    calc_EHH2(locus, m, false); //upstream
    //calc_EHH_downstream(locus, m);
    calc_EHH2(locus, m, true);
    //std::cout<<" ehh1["<<locus<<"]="<<ehh1[locus]<<",ehh0["<<locus<<"]="<<ehh0[locus]<<endl;
    
    //handle all corner cases
    if(hm.all_positions[locus].size()==0){
        iHH1[locus] = 1;
    }
    if(hm.all_positions[locus].size()==numHaps){ //iHH0[locus]==0
        iHH0[locus] = 1;
    }
    if(hm.all_positions[locus].size()==1){
        if(locus == 0 or locus == numSnps-1){
            iHH1[locus] = 0.5;
        }else{
            iHH1[locus] = 1;
        }
    }
    if(hm.all_positions[locus].size()==numHaps-1){
        if(locus == 0 or locus == numSnps-1){
            iHH0[locus] = 0.5;
        }else{
            iHH0[locus] = 1;
        }
    }
}

void EHH::exec() {
	this->init(); //initialize
    cout<<"Number of valid loci: "<<numSnps<<endl;
    cout<<"Number of haplotypes: "<<numHaps<<endl;
    if(calc_all){
        calc_iHS();
    }else{
        calc_EHH(locus);
    }
}

inline unsigned int twice_num_pair(int n){
    return n*n - n;
}

inline unsigned int num_pair(int n){
    return (n*n - n)/2;
}

void EHH::calc_EHH2(int locus, map<int, vector<int> > & m, bool downstream){

    int total_iteration_of_m = 0;
    uint64_t ehh0_before_norm = 0;
    uint64_t ehh1_before_norm = 0;

    bool gap_skip = false;

    uint64_t curr_ehh0_before_norm = 0;
    uint64_t curr_ehh1_before_norm = 0;

    uint64_t n_c0=0;
    uint64_t n_c1=0;

    uint64_t core_n_c0=0;
    uint64_t core_n_c1=0;

    int group_count[numHaps];
    int group_id[numHaps];
    bool isDerived[numHaps]; 

    //NEW
    int RE[numHaps];


    //will be vectorized
    for(int i = 0; i<numHaps; i++){
        group_count[i] = 0;
        group_id[i] = 0;
        isDerived[i] = false;
        RE[i] = -1;
    }

    int totgc=0;
    vector<unsigned int> v = hm.all_positions[locus];

    if(v.size()==0){
        n_c0 = numHaps;
        group_count[0] = numHaps;
        totgc+=1;
        ehh0_before_norm = twice_num_pair(n_c0);
    }else if (v.size()==numHaps){ // all set
        group_count[0] = numHaps;
        totgc+=1;
        n_c1 = numHaps;
        
        for (int set_bit_pos : v){
            isDerived[set_bit_pos] = true;
        }
        ehh1_before_norm = twice_num_pair(n_c1);
    }else{
        group_count[1] = v.size();
        group_count[0] = numHaps - v.size();
        n_c1 = v.size();
        n_c0 = numHaps - v.size();

        for (int set_bit_pos : v){
            isDerived[set_bit_pos] = true;
            group_id[set_bit_pos] = 1;
        }
        
        totgc+=2;
        ehh0_before_norm = twice_num_pair(n_c0);
        ehh1_before_norm = twice_num_pair(n_c1);
    }
    if(downstream){
        if(!calc_all)
            /*
            cout<<"Iter "<<0<<": EHH1["<<locus<<","<<locus<<"]="<<1<<" "<<1<<endl;
            */

        if(twice_num_pair(n_c1)!=0){
            iHH1[locus] += (curr_ehh1_before_norm + ehh1_before_norm) * 0.5 / twice_num_pair(n_c1);
            //cout<<"Summing "<<1.0*curr_ehh1_before_norm/n_c1_squared_minus<<" "<<1.0*ehh1_before_norm/n_c1_squared_minus<<endl;
        }
        if(twice_num_pair(n_c0)!=0){
            iHH0[locus] += (curr_ehh0_before_norm + ehh0_before_norm) * 0.5 / twice_num_pair(n_c0);
        }
    }
    
    curr_ehh1_before_norm = ehh1_before_norm;
    curr_ehh0_before_norm = ehh0_before_norm;

    
    int i = locus;  
    while(true){ // Upstream: for ( int i = locus+1; i<all_positions.size(); i++ )
        if(downstream){
            if (--i < 0) break;
            //if (hm.mentries[locus].phyPos - hm.mentries[i].phyPos > max_extend) break;
        }else{
            if (++i >= numSnps) break;
            //if (hm.mentries[i].phyPos -hm.mentries[locus].phyPos > max_extend) break;
        }
        // if(curr_ehh1_before_norm*1.0/n_c1_squared_minus < cutoff and curr_ehh0_before_norm*1.0/n_c0_squared_minus < cutoff){
        //     break;
        // }
        
        
        //if(curr_ehh1_before_norm*1.0/n_c1_squared_minus < cutoff and curr_ehh0_before_norm*1.0/n_c0_squared_minus < cutoff){
        if(curr_ehh1_before_norm*1.0/twice_num_pair(n_c1) < cutoff or curr_ehh0_before_norm*1.0/twice_num_pair(n_c0)  < cutoff){
        
            //std::cout<<"breaking"<<endl;
            break;
        }
        int distance;
        
        if(downstream){
            distance = hm.mentries[i+1].phyPos - hm.mentries[i].phyPos;
        }else{
            distance = hm.mentries[i].phyPos - hm.mentries[i-1].phyPos;
        }
        assert(distance>=0);
        // if(distance> max_gap){
        //     gap_skip = true;
        //     break;
        // }
        // if(distance> gap_scale){
        //     distance /= gap_scale;
        // }
        
        //distance = 1;
        
        //OPT IDEA: INSTEAD OF MAP JUST USE ARR OF VECTOR
        v = hm.all_positions[i];

        assert(!(hm.all_positions[i].size()==0 or hm.all_positions[i].size()==numHaps));
        if(hm.all_positions[i].size()==0 or hm.all_positions[i].size()==numHaps){
            //cout<<"SKIPPING MONOMORPHIC SITE"<<endl;
            // if(!calc_all)
            //     cout<<"Mono: Iter "<<i-locus<<": EHH1["<<locus<<","<<i<<"]="<<curr_ehh1_before_norm*1.0/n_c1_squared_minus<<","<<curr_ehh0_before_norm*1.0/n_c0_squared_minus<<endl;
        
            
            if(twice_num_pair(n_c1)!=0){
                iHH1[locus] += (curr_ehh1_before_norm + ehh1_before_norm) * distance * 0.5 / twice_num_pair(n_c1) ;
                //cout<<"Summing "<<1.0*curr_ehh1_before_norm/n_c1_squared_minus<<" "<<1.0*ehh1_before_norm/n_c1_squared_minus<<endl;
            }
            if(twice_num_pair(n_c0)!=0){
                iHH0[locus] += (curr_ehh0_before_norm + ehh0_before_norm) * distance * 0.5 / twice_num_pair(n_c0)  ;
            }
            continue;
        }

        for (int set_bit_pos : v){
            int old_group_id = group_id[set_bit_pos];
            
            if(haplo_enabled && not downstream ){
                if(i!=0){
                    //right extreme is reached
                    if(i-1==RE[old_group_id] && RE[old_group_id]!=-1){
                        RE[old_group_id] = -1;
                    }
                }
                
                if(RE[old_group_id] == -1){
                    if(mphbs[set_bit_pos][locus].empty()){ // if that locus is empty for all the set bit pos then this could have been invoked
                        //RE[old_group_id] == -2; //come back later
                    }else{
                        int numSeqInGroup = mphbs[set_bit_pos][locus].top().first; 
                        
                        if(group_count[old_group_id] == numSeqInGroup){
                            int right = mphbs[set_bit_pos][locus].top().second; 
                            
                            if(!calc_all){
                                cout<<old_group_id<<" "<<numSeqInGroup<<" "<<right<<" "<<set_bit_pos<<endl;
                            }

                            if( right+1  != i){
                                //numSeqInGroup= mphbs[set_bit_pos][locus].top().second; 
                                RE[old_group_id] = right;


                                m[old_group_id].clear();
                                //m.erase (m.find (old_group_id));
                                //m.erase (m.find (old_group_id));
                                //m[old_group_id].clear(); // because there is only one hap per group giving this information, it is possible we did not find out until this position that this group should not be checked
                            }
                            mphbs[set_bit_pos][locus].pop();    
                            
                            
                        }
                        while(group_count[old_group_id] < numSeqInGroup ){
                            mphbs[set_bit_pos][locus].pop();
                            numSeqInGroup = mphbs[set_bit_pos][locus].top().first; 
                            
                            if(mphbs[set_bit_pos][locus].empty())
                                break;
                        }
                        
                    }
                }else{
                                //++total_iteration_of_m;
                    
                    //skip
                    // if(i < RE[old_group_id] ){
                    // //Skip calculation this group
                    // }else{
                    //     RE[old_group_id] = -1;
                    //     //mphbs[set_bit_pos][locus].pop();
                    //     //Update RightExtreme
                    // }                 
                }
                
                // still -1 even after doing the abve step
                if(RE[old_group_id]==-1){ //opt:  group_count[old_group_id] > 1
                    m[old_group_id].push_back(set_bit_pos);
                    //++total_iteration_of_m;

                }
                
            }

            if(!haplo_enabled){
                //++total_iteration_of_m;

                m[old_group_id].push_back(set_bit_pos);
            }
        }

        for (const auto &ele : m) {
            
            int old_group_id = ele.first;
            int newgroup_size = ele.second.size() ;
                            
            total_iteration_of_m += newgroup_size;
                            
            if(group_count[old_group_id] == newgroup_size || newgroup_size == 0){
                continue;
            }

            //cout<<old_group_id<<" (";
            for(int v: ele.second){
                //++total_iteration_of_m;
                //cout<<v<<" ";

                group_id[v] = totgc;
            }
            //cout<<endl;
            
            int del_update = -twice_num_pair(group_count[old_group_id]) + twice_num_pair(newgroup_size) + twice_num_pair(group_count[old_group_id] - newgroup_size);
            group_count[old_group_id] -= newgroup_size;
            group_count[totgc] += newgroup_size;
            

            
            // if(haplo_enabled && not downstream && RE[old_group_id] == -1){
            //     int numSeqInGroup = mphbs[ele.second[0]][locus].top().first; 
            //     if(group_count[totgc] == numSeqInGroup){
            //         int right = mphbs[ele.second[0]][locus].top().second; 
            //         RE[old_group_id] = right;
            //     }
            // }

            totgc+=1;

            bool isDerivedGroup =  isDerived[ele.second[0]]; // just check first element to know if it is derived. 
            if(isDerivedGroup)
            {
                ehh1_before_norm += del_update;
            }else{
                ehh0_before_norm += del_update;
            }
            
        }

        if(twice_num_pair(n_c1)!=0){
            iHH1[locus] += (curr_ehh1_before_norm + ehh1_before_norm) * distance * 0.5 / twice_num_pair(n_c1);
            //cout<<"Summing "<<1.0*curr_ehh1_before_norm/n_c1_squared_minus<<" "<<1.0*ehh1_before_norm/n_c1_squared_minus<<endl;
        }

        if(twice_num_pair(n_c0)!=0){
            iHH0[locus] += (curr_ehh0_before_norm + ehh0_before_norm) * distance * 0.5 / twice_num_pair(n_c0);
        }

        curr_ehh1_before_norm = ehh1_before_norm;
        curr_ehh0_before_norm = ehh0_before_norm;

        // ehh1[locus] = 1.0*ehh1_before_norm/n_c1_squared_minus;
        // ehh0[locus] = 1.0*ehh0_before_norm/n_c0_squared_minus;

        //this shouldn't execute
        if(twice_num_pair(n_c1)==0){
            curr_ehh1_before_norm = 0;
        }
        if(twice_num_pair(n_c0)==0){
            curr_ehh0_before_norm = 0;
        }


        if(!calc_all){
            std::cout<<"Iter "<<i-locus<<": EHH1["<<locus<<","<<i<<"]="<<curr_ehh1_before_norm*1.0/twice_num_pair(n_c1)<<","<<curr_ehh0_before_norm*1.0/twice_num_pair(n_c0)<<" ";
            
            for (int x = 0 ; x < totgc;  x++){
                cout<<group_count[x]<<"("<<x<<")";
            }
            
           std::cout<<endl;
        }
        
        //logg[tid]+="map size before"+to_string(m.size())+"\n";
        //cout<< logg[tid];
        m.clear();
        //logg[tid]+="map size after"+to_string(m.size())+"\n";
        //cout<< logg[tid];
        // if(gap_skip==true)
        //     break;
    }
    if(!downstream) cout<<"Total_iter_m "<<total_iteration_of_m<<endl;
}


void EHH::calc_EHH2_v2(int locus, map<int, vector<int> > & m, bool downstream){

    map<int, priority_queue<pair<int, int> > >& mphbs_l = mphbs[locus];
    int total_iteration_of_m = 0;
    uint64_t ehh0_before_norm = 0;
    uint64_t ehh1_before_norm = 0;

    bool gap_skip = false;

    uint64_t curr_ehh0_before_norm = 0;
    uint64_t curr_ehh1_before_norm = 0;

    uint64_t n_c0=0;
    uint64_t n_c1=0;

    uint64_t core_n_c0=0;
    uint64_t core_n_c1=0;

    int group_count[numHaps];
    int group_id[numHaps];
    bool isDerived[numHaps]; 

    //NEW
    int RE[numHaps];


    //will be vectorized
    for(int i = 0; i<numHaps; i++){
        group_count[i] = 0;
        group_id[i] = 0;
        isDerived[i] = false;
        RE[i] = -1;
    }

    int totgc=0;
    vector<unsigned int> v = hm.all_positions[locus];

    if(v.size()==0){
        n_c0 = numHaps;
        group_count[0] = numHaps;
        totgc+=1;
        ehh0_before_norm = twice_num_pair(n_c0);
    }else if (v.size()==numHaps){ // all set
        group_count[0] = numHaps;
        totgc+=1;
        n_c1 = numHaps;
        
        for (int set_bit_pos : v){
            isDerived[set_bit_pos] = true;
        }
        ehh1_before_norm = twice_num_pair(n_c1);
    }else{
        group_count[1] = v.size();
        group_count[0] = numHaps - v.size();
        n_c1 = v.size();
        n_c0 = numHaps - v.size();

        for (int set_bit_pos : v){
            isDerived[set_bit_pos] = true;
            group_id[set_bit_pos] = 1;
        }
        
        totgc+=2;
        ehh0_before_norm = twice_num_pair(n_c0);
        ehh1_before_norm = twice_num_pair(n_c1);
    }
    if(downstream){
        if(!calc_all)
            /*
            cout<<"Iter "<<0<<": EHH1["<<locus<<","<<locus<<"]="<<1<<" "<<1<<endl;
            */

        if(twice_num_pair(n_c1)!=0){
            iHH1[locus] += (curr_ehh1_before_norm + ehh1_before_norm) * 0.5 / twice_num_pair(n_c1);
            //cout<<"Summing "<<1.0*curr_ehh1_before_norm/n_c1_squared_minus<<" "<<1.0*ehh1_before_norm/n_c1_squared_minus<<endl;
        }
        if(twice_num_pair(n_c0)!=0){
            iHH0[locus] += (curr_ehh0_before_norm + ehh0_before_norm) * 0.5 / twice_num_pair(n_c0);
        }
    }
    
    curr_ehh1_before_norm = ehh1_before_norm;
    curr_ehh0_before_norm = ehh0_before_norm;

    
    int i = locus;  
    while(true){ // Upstream: for ( int i = locus+1; i<all_positions.size(); i++ )
        if(downstream){
            if (--i < 0) break;
            //if (hm.mentries[locus].phyPos - hm.mentries[i].phyPos > max_extend) break;
        }else{
            if (++i >= numSnps) break;
            //if (hm.mentries[i].phyPos -hm.mentries[locus].phyPos > max_extend) break;
        }
        // if(curr_ehh1_before_norm*1.0/n_c1_squared_minus < cutoff and curr_ehh0_before_norm*1.0/n_c0_squared_minus < cutoff){
        //     break;
        // }
        
        
        //if(curr_ehh1_before_norm*1.0/n_c1_squared_minus < cutoff and curr_ehh0_before_norm*1.0/n_c0_squared_minus < cutoff){
        if(curr_ehh1_before_norm*1.0/twice_num_pair(n_c1) < cutoff or curr_ehh0_before_norm*1.0/twice_num_pair(n_c0)  < cutoff){
        
            //std::cout<<"breaking"<<endl;
            break;
        }
        int distance;
        
        if(downstream){
            distance = hm.mentries[i+1].phyPos - hm.mentries[i].phyPos;
        }else{
            distance = hm.mentries[i].phyPos - hm.mentries[i-1].phyPos;
        }
        assert(distance>=0);
        // if(distance> max_gap){
        //     gap_skip = true;
        //     break;
        // }
        // if(distance> gap_scale){
        //     distance /= gap_scale;
        // }
        
        //distance = 1;
        
        //OPT IDEA: INSTEAD OF MAP JUST USE ARR OF VECTOR
        v = hm.all_positions[i];

        assert(!(hm.all_positions[i].size()==0 or hm.all_positions[i].size()==numHaps));
        if(hm.all_positions[i].size()==0 or hm.all_positions[i].size()==numHaps){
            //cout<<"SKIPPING MONOMORPHIC SITE"<<endl;
            // if(!calc_all)
            //     cout<<"Mono: Iter "<<i-locus<<": EHH1["<<locus<<","<<i<<"]="<<curr_ehh1_before_norm*1.0/n_c1_squared_minus<<","<<curr_ehh0_before_norm*1.0/n_c0_squared_minus<<endl;
        
            
            if(twice_num_pair(n_c1)!=0){
                iHH1[locus] += (curr_ehh1_before_norm + ehh1_before_norm) * distance * 0.5 / twice_num_pair(n_c1) ;
                //cout<<"Summing "<<1.0*curr_ehh1_before_norm/n_c1_squared_minus<<" "<<1.0*ehh1_before_norm/n_c1_squared_minus<<endl;
            }
            if(twice_num_pair(n_c0)!=0){
                iHH0[locus] += (curr_ehh0_before_norm + ehh0_before_norm) * distance * 0.5 / twice_num_pair(n_c0)  ;
            }
            continue;
        }

        for (int set_bit_pos : v){
            int old_group_id = group_id[set_bit_pos];
            
            if(haplo_enabled && not downstream ){
                if(i!=0){
                    //right extreme is reached
                    if(i-1==RE[old_group_id] && RE[old_group_id]!=-1){
                        RE[old_group_id] = -1;
                    }
                }
                
                if(RE[old_group_id] == -1){
                    if(!mphbs_l.count(set_bit_pos)){ // if that locus is empty for all the set bit pos then this could have been invoked
                        //RE[old_group_id] == -2; //come back later
                    }else{
                        int numSeqInGroup = mphbs_l[set_bit_pos].top().first; 
                        
                        if(group_count[old_group_id] == numSeqInGroup){
                            int right = mphbs_l[set_bit_pos].top().second; 
                            
                            if(!calc_all){
                                cout<<old_group_id<<" "<<numSeqInGroup<<" "<<right<<" "<<set_bit_pos<<endl;
                            }

                            if( right+1  != i){
                                //numSeqInGroup= mphbs[set_bit_pos][locus].top().second; 
                                RE[old_group_id] = right;


                                m[old_group_id].clear();
                                m.erase (m.find (old_group_id));
                                //m.erase (m.find (old_group_id));
                                //m[old_group_id].clear(); // because there is only one hap per group giving this information, it is possible we did not find out until this position that this group should not be checked
                            }
                            mphbs_l[set_bit_pos].pop();    
                            
                            
                        }
                        while(group_count[old_group_id] < numSeqInGroup ){
                            mphbs_l[set_bit_pos].pop();
                            numSeqInGroup = mphbs_l[set_bit_pos].top().first; 
                            
                            if(mphbs[set_bit_pos].empty())
                                break;
                        }
                        
                    }
                }else{
                                //++total_iteration_of_m;
                    
                    //skip
                    // if(i < RE[old_group_id] ){
                    // //Skip calculation this group
                    // }else{
                    //     RE[old_group_id] = -1;
                    //     //mphbs[set_bit_pos][locus].pop();
                    //     //Update RightExtreme
                    // }                 
                }
                
                // still -1 even after doing the abve step
                if(RE[old_group_id]==-1){ //opt:  group_count[old_group_id] > 1
                    m[old_group_id].push_back(set_bit_pos);
                    //++total_iteration_of_m;

                }
                
            }

            if(!haplo_enabled){
                //++total_iteration_of_m;

                m[old_group_id].push_back(set_bit_pos);
            }
        }

        for (const auto &ele : m) {
            
            int old_group_id = ele.first;
            int newgroup_size = ele.second.size() ;
                            
            total_iteration_of_m += newgroup_size;
                            
            if(group_count[old_group_id] == newgroup_size || newgroup_size == 0){
                continue;
            }

            //cout<<old_group_id<<" (";
            for(int v: ele.second){
                //++total_iteration_of_m;
                //cout<<v<<" ";

                group_id[v] = totgc;
            }
            //cout<<endl;
            
            int del_update = -twice_num_pair(group_count[old_group_id]) + twice_num_pair(newgroup_size) + twice_num_pair(group_count[old_group_id] - newgroup_size);
            group_count[old_group_id] -= newgroup_size;
            group_count[totgc] += newgroup_size;
            

            
            // if(haplo_enabled && not downstream && RE[old_group_id] == -1){
            //     int numSeqInGroup = mphbs[ele.second[0]][locus].top().first; 
            //     if(group_count[totgc] == numSeqInGroup){
            //         int right = mphbs[ele.second[0]][locus].top().second; 
            //         RE[old_group_id] = right;
            //     }
            // }

            totgc+=1;

            bool isDerivedGroup =  isDerived[ele.second[0]]; // just check first element to know if it is derived. 
            if(isDerivedGroup)
            {
                ehh1_before_norm += del_update;
            }else{
                ehh0_before_norm += del_update;
            }
            
        }

        if(twice_num_pair(n_c1)!=0){
            iHH1[locus] += (curr_ehh1_before_norm + ehh1_before_norm) * distance * 0.5 / twice_num_pair(n_c1);
            //cout<<"Summing "<<1.0*curr_ehh1_before_norm/n_c1_squared_minus<<" "<<1.0*ehh1_before_norm/n_c1_squared_minus<<endl;
        }

        if(twice_num_pair(n_c0)!=0){
            iHH0[locus] += (curr_ehh0_before_norm + ehh0_before_norm) * distance * 0.5 / twice_num_pair(n_c0);
        }

        curr_ehh1_before_norm = ehh1_before_norm;
        curr_ehh0_before_norm = ehh0_before_norm;

        // ehh1[locus] = 1.0*ehh1_before_norm/n_c1_squared_minus;
        // ehh0[locus] = 1.0*ehh0_before_norm/n_c0_squared_minus;

        //this shouldn't execute
        if(twice_num_pair(n_c1)==0){
            curr_ehh1_before_norm = 0;
        }
        if(twice_num_pair(n_c0)==0){
            curr_ehh0_before_norm = 0;
        }


        if(!calc_all){
            std::cout<<"Iter "<<i-locus<<": EHH1["<<locus<<","<<i<<"]="<<curr_ehh1_before_norm*1.0/twice_num_pair(n_c1)<<","<<curr_ehh0_before_norm*1.0/twice_num_pair(n_c0)<<" ";
            
            for (int x = 0 ; x < totgc;  x++){
                cout<<group_count[x]<<"("<<x<<")";
            }
            
           std::cout<<endl;
        }
        
        //logg[tid]+="map size before"+to_string(m.size())+"\n";
        //cout<< logg[tid];
        m.clear();
        //logg[tid]+="map size after"+to_string(m.size())+"\n";
        //cout<< logg[tid];
        // if(gap_skip==true)
        //     break;
    }
    if(!downstream) cout<<"Total_iter_m "<<total_iteration_of_m<<endl;
}


void EHH::thread_ihs(int tid, map<int, vector<int> >& m, map<int, vector<int> >& md, EHH* ehh_obj){
    int elem_per_block = floor(ehh_obj->numSnps/ehh_obj->numThread);
    int start = tid*elem_per_block ;
    int end = start + elem_per_block  ;
    if(tid == ehh_obj->numThread-1 ){
        end = ehh_obj->numSnps;
    }

    //#pragma omp parallel 
    for(int locus = start; locus< end; locus++){
        
            ehh_obj->calc_EHH(locus);
        
        // ehh_obj->iHH0[locus] = 0;
        // ehh_obj->iHH1[locus] = 0;
        // ehh_obj->calc_EHH2(locus, m, false); //upstream
        // ehh_obj->calc_EHH2(locus, md, true);
        // //std::cout<<" ehh1["<<locus<<"]="<<ehh1[locus]<<",ehh0["<<locus<<"]="<<ehh0[locus]<<endl;
        // if(ehh_obj->all_positions[locus].size()==0){
        //     ehh_obj->iHH1[locus] = 1;
        // }
        // if(ehh_obj->all_positions[locus].size()==ehh_obj->ADVANCED_N){ //iHH0[locus]==0
        //     ehh_obj->iHH0[locus] = 1;
        // }
        // if(ehh_obj->all_positions[locus].size()==1){
        //     if(locus -= 0 or locus == ehh_obj->ADVANCED_D-1){
        //         ehh_obj->iHH1[locus] = 0.5;
        //     }else{
        //         ehh_obj->iHH1[locus] = 1;
        //     }
        // }
        // if(ehh_obj->all_positions[locus].size()==ehh_obj->ADVANCED_N-1){
        //     if(locus == 0 or locus == ehh_obj->ADVANCED_D-1){
        //         ehh_obj->iHH0[locus] = 0.5;
        //     }else{
        //         ehh_obj->iHH0[locus] = 1;
        //     }
        // }
    }
    
    ehh_obj->logg[tid]+="finishing thread #"+to_string(tid)+"\n"; 
    cout<<"finishing thread # "+to_string(tid)+" at "+to_string(readTimer())+"\n";
}

void EHH::calc_iHS(){
    std::map<int, std::vector<int> > map_per_thread[numThread];
    std::map<int, std::vector<int> > mapd_per_thread[numThread];

    if (!openmp_enabled)
    {
        int total_calc_to_be_done = numSnps;
        std::thread *myThreads = new std::thread[numThread];
        for (int i = 0; i < numThread; i++)
        {
            myThreads[i] = std::thread(thread_ihs, i, std::ref(map_per_thread[i]),  std::ref(mapd_per_thread[i]), this);
        }
        for (int i = 0; i < numThread; i++)
        // Join will block our main thread, and so the program won't exit until
        // everyone comes home.
        {
            myThreads[i].join();
        }
        delete[] myThreads;
        cout << "all threads finished. now calculating ihh..." << endl;
    }else{
        // #pragma clang loop unroll_count(8) // 
        // #pragma clang loop vectorize(assume_safety)
        cout<<"open mp enabled. "<<omp_get_max_threads()<<"threads"<<endl;

        // // // #pragma omp parallel shared(hm)
        {
            // // // #pragma omp for schedule(dynamic,10)
            for(int i = 0 ; i< numSnps; i++){
                if(i>50)
                    haplo_enabled = false;
                calc_EHH(i);
                
                //cout<<"open mp enabled. "<<omp_get_num_threads()<<"threads"<<endl;
            }
        }
         cout<<"finishing all threads # at "+to_string(readTimer())+"\n";
    }
   
    out_fp = fopen(output_filename.c_str(),"w");
    for (int i = 0; i < numSnps; i++){
         if(hm.getMAF(i) >= min_maf && 1-hm.getMAF(i) >= min_maf){
            fprintf(out_fp, "%d %d %f %f %f %f\n", hm.mentries[i].phyPos, hm.mentries[i].locId, hm.all_positions[i].size()*1.0/numHaps, iHH1[i], iHH0[i], log10(iHH1[i]/ iHH0[i]));
         }
    }
    fclose(out_fp);
   
    return;
}


