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
    if(hm.unphased)
        delete[] iHH2;

    delete[] logg;
}

void EHH::cli_prepare(CLI::App * app) {
	this->subapp = app->add_subcommand("ehh", "Calculate EHH. version ZERO-ONE.");

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

    CLI::Option * map_opt = subapp->add_option("-m, --map", input_filename_map, "Input MAP file");
	// in_option2->required();
	map_opt->check(CLI::ExistingFile);
    group_opt->needs(map_opt);

    useVCF = true;
    CLI::Option * opt;

	numThread = std::thread::hardware_concurrency(); //may return 0 when not able to detect
	opt = subapp->add_option("-t, --threads", numThread, "The number of threads to spawn during the calculation. Partitions loci across threads. (by default uses all available cores)")->capture_default_str();
        
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

}

void EHH::init(){
    useVCF=true;
    if(useVCF){
        cout<<"Loading "<<input_filename_vcf<<endl;
    	hm.loadVCF(input_filename_vcf.c_str(), min_maf); //populate the matrix
    }else{
        cout<<"Loading "<<input_filename_hap<<endl;
    	hm.loadHapMap(input_filename_hap.c_str(), min_maf); //populate the matrix
    }

    this->numHaps = hm.numHaps();
	this->numSnps = hm.numSnps();

    iHH0 = new double[numSnps];
    iHH1 = new double[numSnps];

    if(hm.unphased)
        iHH2 = new double[numSnps];


    logg = new string[numThread];
}

void EHH::calc_EHH(int locus){
    //map<int, vector<int> > m;
    //map<int, vector<int> > m2;

    iHH0[locus] = 0;
    iHH1[locus] = 0;


    // if(hm.getMAF(locus) < min_maf || 1-hm.getMAF(locus) < min_maf){
    //     return;
    // }
    // ehh1[locus] = 0;
    // ehh0[locus] = 0;


    calc_EHH2(locus, false); //upstream
    //calc_EHH_downstream(locus, m);
    calc_EHH2(locus, true);
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
    std::cout<<"Number of valid loci: "<<numSnps<<endl;
    std::cout<<"Number of haplotypes: "<<numHaps<<endl;
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

//void EHH::calc_EHH2(int locus, map<int, vector<int> > & m, bool downstream){ 
void EHH::calc_EHH2(int locus, bool downstream){
    map<int, vector<int> > m;
    map<int, vector<int> > m2;

    uint64_t ehh0_before_norm = 0;
    uint64_t ehh1_before_norm = 0;
    uint64_t ehh2_before_norm = 0;

    bool gap_skip = false;

    uint64_t curr_ehh0_before_norm = 0;
    uint64_t curr_ehh1_before_norm = 0;
    uint64_t curr_ehh2_before_norm = 0;

    uint64_t n_c0=0;
    uint64_t n_c1=0;
    uint64_t n_c2=0;

    uint64_t core_n_c0=0;
    uint64_t core_n_c1=0;
    uint64_t core_n_c2=0;

    int group_count[numHaps];
    int group_id[numHaps];
    bool isDerived[numHaps];  
    bool isTwoDerived[numHaps];  
    bool isAncestral[numHaps]; 

    //will be vectorized
    for(int i = 0; i<numHaps; i++){
        group_count[i] = 0;
        group_id[i] = 0;
        isDerived[i] = false;
        isTwoDerived[i] = false;
        isAncestral[i] = false;

    }

    int totgc=0;
    vector<unsigned int> v = hm.all_positions[locus];
    vector<unsigned int> v2 = hm.all_positions_two[locus];

    if(v.size()==0){
        cerr<<"MAC cant be zero"<<endl;
        exit(2);
        n_c0 = numHaps;
        group_count[0] = numHaps;
        totgc+=1;
        ehh0_before_norm = twice_num_pair(n_c0);
    }else if (v.size()==numHaps){ // all set
        cerr<<"MAC cant be zero"<<endl;
        exit(2);
        group_count[0] = numHaps;
        totgc+=1;
        n_c1 = numHaps;
        
        for (int set_bit_pos : v){
            isDerived[set_bit_pos] = true;
        }
        ehh1_before_norm = twice_num_pair(n_c1);

    }else{
        if(hm.data[locus].flipped){
            group_count[1] = v.size();
            group_count[0] = numHaps - v.size();
            n_c0 = v.size();
            n_c1 = numHaps - v.size();

            for (int set_bit_pos : v){
                isAncestral[set_bit_pos] = true;
                group_id[set_bit_pos] = 1;
            }
            // for(int i=0; i<numHaps; i++){
            //     if(!isAncestral[i]){
            //         isDerived[i] = true;
            //     }
            // }
        }else{
            if(hm.unphased){
                group_count[2] = v2.size();
                group_count[1] = v.size();
                group_count[0] = numHaps - v.size();
                n_c2 = v2.size();
                n_c1 = v.size();
                n_c0 = numHaps - n_c1 - n_c2;

                for (int set_bit_pos : v){
                    isDerived[set_bit_pos] = true;
                    group_id[set_bit_pos] = 1;
                }

                for (int set_bit_pos : v2){
                    isTwoDerived[set_bit_pos] = true;
                    group_id[set_bit_pos] = 2;
                }
                totgc+=1;

            }else{
                group_count[1] = v.size();
                group_count[0] = numHaps - v.size();
                n_c1 = v.size();
                n_c0 = numHaps - v.size();

                for (int set_bit_pos : v){
                    isDerived[set_bit_pos] = true;
                    group_id[set_bit_pos] = 1;
                }
            }
            
        }
        
        totgc+=2;
        ehh0_before_norm = twice_num_pair(n_c0);
        ehh1_before_norm = twice_num_pair(n_c1);
        ehh2_before_norm = twice_num_pair(n_c1);

    }


    if(downstream){
        if(!calc_all)
            std::cout<<"Iter "<<0<<": EHH1["<<locus<<","<<locus<<"]="<<1<<" "<<1<<endl;
        

        if(twice_num_pair(n_c2)!=0){
            iHH2[locus] += (curr_ehh2_before_norm + ehh2_before_norm) * 0.5 / twice_num_pair(n_c2);
            //cout<<"Summing "<<1.0*curr_ehh1_before_norm/n_c1_squared_minus<<" "<<1.0*ehh1_before_norm/n_c1_squared_minus<<endl;
        }

        if(twice_num_pair(n_c1)!=0){
            iHH1[locus] += (curr_ehh1_before_norm + ehh1_before_norm) * 0.5 / twice_num_pair(n_c1);
            //cout<<"Summing "<<1.0*curr_ehh1_before_norm/n_c1_squared_minus<<" "<<1.0*ehh1_before_norm/n_c1_squared_minus<<endl;
        }
        if(twice_num_pair(n_c0)!=0){
            iHH0[locus] += (curr_ehh0_before_norm + ehh0_before_norm) * 0.5 / twice_num_pair(n_c0);
        }
    }
    
    curr_ehh2_before_norm = ehh2_before_norm;
    curr_ehh1_before_norm = ehh1_before_norm;
    curr_ehh0_before_norm = ehh0_before_norm;

    
    int i = locus;  
    while(true){ // Upstream: for ( int i = locus+1; i<all_positions.size(); i++ )
        if(downstream){
            if (--i < 0) break;
            //if (hm.data[locus].phyPos - hm.data[i].phyPos > max_extend) break;
        }else{
            if (++i >= numSnps) break;
            //if (hm.data[i].phyPos -hm.data[locus].phyPos > max_extend) break;
        }
        // if(curr_ehh1_before_norm*1.0/n_c1_squared_minus < cutoff and curr_ehh0_before_norm*1.0/n_c0_squared_minus < cutoff){
        //     break;
        // }
        
        
        //if(curr_ehh1_before_norm*1.0/n_c1_squared_minus < cutoff and curr_ehh0_before_norm*1.0/n_c0_squared_minus < cutoff){
        
        if(hm.unphased){
            if(curr_ehh2_before_norm*1.0/twice_num_pair(n_c2) < cutoff or curr_ehh1_before_norm*1.0/twice_num_pair(n_c1) < cutoff or curr_ehh0_before_norm*1.0/twice_num_pair(n_c0)  < cutoff){
                break;
            }
        }else{
            if(curr_ehh1_before_norm*1.0/twice_num_pair(n_c1) < cutoff or curr_ehh0_before_norm*1.0/twice_num_pair(n_c0)  < cutoff){
                break;
            }
        }
        
        int distance;
        
        if(downstream){
            distance = hm.data[i+1].phyPos - hm.data[i].phyPos;
        }else{
            distance = hm.data[i].phyPos - hm.data[i-1].phyPos;
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

        if(hm.unphased){
            if(n_c0 == numHaps or n_c1 == numHaps or n_c2 == numHaps){
                cerr<<"WARNING: Monomorphic site."<<endl;
                //exit(2);
            
                // if(!calc_all)
                //     cout<<"Mono: Iter "<<i-locus<<": EHH1["<<locus<<","<<i<<"]="<<curr_ehh1_before_norm*1.0/n_c1_squared_minus<<","<<curr_ehh0_before_norm*1.0/n_c0_squared_minus<<endl;
                if(twice_num_pair(n_c2)!=0){
                    iHH2[locus] += (curr_ehh2_before_norm + ehh2_before_norm) * distance * 0.5 / twice_num_pair(n_c2) ;
                    //cout<<"Summing "<<1.0*curr_ehh1_before_norm/n_c1_squared_minus<<" "<<1.0*ehh1_before_norm/n_c1_squared_minus<<endl;
                }
                if(twice_num_pair(n_c1)!=0){
                    iHH1[locus] += (curr_ehh1_before_norm + ehh1_before_norm) * distance * 0.5 / twice_num_pair(n_c1) ;
                    //cout<<"Summing "<<1.0*curr_ehh1_before_norm/n_c1_squared_minus<<" "<<1.0*ehh1_before_norm/n_c1_squared_minus<<endl;
                }
                if(twice_num_pair(n_c0)!=0){
                    iHH0[locus] += (curr_ehh0_before_norm + ehh0_before_norm) * distance * 0.5 / twice_num_pair(n_c0)  ;
                }
                continue;
            }
        }else{
            if(hm.all_positions[i].size()==0 or hm.all_positions[i].size()==numHaps){
                cerr<<"WARNING: Monomorphic site."<<endl;
                if(twice_num_pair(n_c1)!=0){
                    iHH1[locus] += (curr_ehh1_before_norm + ehh1_before_norm) * distance * 0.5 / twice_num_pair(n_c1) ;
                    //cout<<"Summing "<<1.0*curr_ehh1_before_norm/n_c1_squared_minus<<" "<<1.0*ehh1_before_norm/n_c1_squared_minus<<endl;
                }
                if(twice_num_pair(n_c0)!=0){
                    iHH0[locus] += (curr_ehh0_before_norm + ehh0_before_norm) * distance * 0.5 / twice_num_pair(n_c0)  ;
                }
                continue;
            }
        }
        

        for (int set_bit_pos : v){
            int old_group_id = group_id[set_bit_pos];
            m[old_group_id].push_back(set_bit_pos);

            //if(MPHAPLOBLOCKS)

        }
        if(hm.unphased){
            set<int> seen_groups;
            for (int set_bit_pos : v2){
                int old_group_id = group_id[set_bit_pos];
                m2[old_group_id].push_back(set_bit_pos);
             }
             for (const auto &ele : m) {
                int old_group_id = ele.first;
                seen_groups.insert(old_group_id);
                int newgroup_size = ele.second.size() ;
                
                if(group_count[old_group_id] == newgroup_size  ){ // all 1
                    continue;
                }

                //handle groups 120, 12, and 10 ::: remaining 20, 

                // this block is pointless
                if( newgroup_size == 0 ){
                    cerr<<"Should never execute"<<endl;
                    exit(2);
                    if( m2[old_group_id].size() == 0 || group_count[old_group_id] == m2[old_group_id].size() ){  //all 0, all 2
                        continue;
                    }
                }

                if(m2[old_group_id].size() != 0 ){  // both 1 and 2 has something
                    // if there is no 0, then dont do ++totgc (case 12)
                    if( m2[old_group_id].size() + newgroup_size == group_count[old_group_id] ){ //case 12
                        //DO nothing
                    }else{  // case 120
                        for(int v: m2[old_group_id]){
                            group_id[v] = totgc;
                        }
                        ++totgc;
                    }
                }   
                    

                for(int v: ele.second){
                    group_id[v] = totgc;
                }
                
                int del_update = -twice_num_pair(group_count[old_group_id]) + twice_num_pair(newgroup_size) + twice_num_pair(group_count[old_group_id] - newgroup_size);
                //int del_update = -(newgroup_size)*(group_count[old_group_id]-newgroup_size);
                
                group_count[old_group_id] -= newgroup_size;
                group_count[totgc] += newgroup_size;
                totgc+=1;

                

                bool isDerivedGroup =  isDerived[ele.second[0]];
                bool isTwoDerivedGroup = isTwoDerived[ele.second[0]];
                if(isDerivedGroup)
                {
                    ehh2_before_norm += del_update;
                }else if (isTwoDerivedGroup){
                    ehh1_before_norm += del_update;

                }else{
                    ehh0_before_norm += del_update;

                }
            }

             for (const auto &ele : m2) {
                    int old_group_id = ele.first;
                    if(seen_groups.count(old_group_id)!=0){
                        continue;
                    }
                    int newgroup_size = ele.second.size() ;
                    
                    if(group_count[old_group_id] == newgroup_size  ){ // all 2
                        continue;
                    }

                    for(int v: ele.second){
                        group_id[v] = totgc;
                    }
                    
                    int del_update = -twice_num_pair(group_count[old_group_id]) + twice_num_pair(newgroup_size) + twice_num_pair(group_count[old_group_id] - newgroup_size); 
                    group_count[old_group_id] -= newgroup_size;
                    group_count[totgc] += newgroup_size;
                    totgc+=1;

                    bool isDerivedGroup =  isDerived[ele.second[0]];
                    bool isTwoDerivedGroup = isTwoDerived[ele.second[0]];
                    if(isDerivedGroup)
                    {
                        ehh2_before_norm += del_update;
                    }else if (isTwoDerivedGroup){
                        ehh1_before_norm += del_update;
                    }else{
                        ehh0_before_norm += del_update;
                    }
                }
        }else{
            for (const auto &ele : m) {
                int old_group_id = ele.first;
                int newgroup_size = ele.second.size() ;
                
                if(group_count[old_group_id] == newgroup_size || newgroup_size == 0){
                    continue;
                }
                
                for(int v: ele.second){
                    group_id[v] = totgc;
                }
                
                int del_update = -twice_num_pair(group_count[old_group_id]) + twice_num_pair(newgroup_size) + twice_num_pair(group_count[old_group_id] - newgroup_size);
                //int del_update = -(newgroup_size)*(group_count[old_group_id]-newgroup_size);
                
                group_count[old_group_id] -= newgroup_size;
                group_count[totgc] += newgroup_size;
                totgc+=1;

                bool isDerivedGroup =  (!hm.data[locus].flipped && isDerived[ele.second[0]]) || (hm.data[locus].flipped && !isAncestral[ele.second[0]]); // just check first element to know if it is derived. 
                if(isDerivedGroup)
                {
                    ehh1_before_norm += del_update;
                }else{
                    ehh0_before_norm += del_update;
                }
                
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


        if(!calc_all)
            std::cout<<"Iter "<<i-locus<<": EHH1["<<locus<<","<<i<<"]="<<curr_ehh1_before_norm*1.0/twice_num_pair(n_c1)<<","<<curr_ehh0_before_norm*1.0/twice_num_pair(n_c0)<<endl;
        
        //logg[tid]+="map size before"+to_string(m.size())+"\n";
        //cout<< logg[tid];
        m.clear();
        m2.clear();

        //logg[tid]+="map size after"+to_string(m.size())+"\n";
        //cout<< logg[tid];
        // if(gap_skip==true)
        //     break;
    }
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

        #pragma omp parallel shared(hm)
        {
            #pragma omp for schedule(dynamic,10)
            for(int i = 0 ; i< numSnps; i++){
                calc_EHH(i);
                //cout<<"open mp enabled. "<<omp_get_num_threads()<<"threads"<<endl;
            }
        }
         cout<<"finishing all threads # at "+to_string(readTimer())+"\n";
    }
   
    out_fp = fopen(output_filename.c_str(),"w");
    for (int i = 0; i < numSnps; i++){
         if(hm.get1F(i) >= min_maf && 1-hm.get1F(i) >= min_maf){
            fprintf(out_fp, "%d %d %f %f %f %f\n", hm.data[i].phyPos, hm.data[i].locId, hm.get1F(i), iHH1[i], iHH0[i], log10(iHH1[i]/ iHH0[i]));
         }
    }
    fclose(out_fp);
   
    return;
}


