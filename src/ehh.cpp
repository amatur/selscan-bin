#include "ehh.hpp"
#include "benchmark.hpp"
#include "hapmap.hpp"
#include "utils.hpp"

EHH::EHH(){
}
EHH::~EHH(){
    delete[] iHH0;
    delete[] iHH1;
    delete[] logg;
}

void EHH::cli_prepare(CLI::App * app) {
	this->subapp = app->add_subcommand("ehh", "Calculate EHH.");

    CLI::Option * opt;
	opt = subapp->add_option("-p, --hap", input_filename_hap, "A hapfile with one row per haplotype, \n"
                                                                "and one column per variant.\n"
                                                                "Variants should be coded 0/1")->capture_default_str();
	//opt->required();
	opt->check(CLI::ExistingFile);


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
    opt->option_text("<float> [0.05]");


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

    

	app->get_formatter()->column_width(10);
}

void EHH::init(){
    HapMap hm;
    cout<<"Loading "<<input_filename_hap<<endl;
	hm.loadHap(input_filename_hap.c_str(), min_maf,this->all_positions, this->loc_map); //populate the hapmap
	this->ADVANCED_N = hm.numHaps();
	this->ADVANCED_D = hm.numSnps();

    // ehh0 = new float[ADVANCED_D];
    // ehh1 = new float[ADVANCED_D];
    // ehh0_downstream  = new float[ADVANCED_D];
    // ehh1_downstream = new float[ADVANCED_D];

    iHH0 = new float[ADVANCED_D];
    iHH1 = new float[ADVANCED_D];
    // iHH0_downstream = new float[ADVANCED_D];
    // iHH1_downstream = new float[ADVANCED_D];

    logg = new string[numThread];

}

void EHH::calc_EHH(int locus){
    map<int, vector<int> > m;
    iHH0[locus] = 0;
    iHH1[locus] = 0;
    // ehh1[locus] = 0;
    // ehh0[locus] = 0;
    calc_EHH2(locus, m, false); //upstream
    //calc_EHH_downstream(locus, m);
    calc_EHH2(locus, m, true);
    //std::cout<<" ehh1["<<locus<<"]="<<ehh1[locus]<<",ehh0["<<locus<<"]="<<ehh0[locus]<<endl;
    if(all_positions[locus].size()==0){
        iHH1[locus] = 1;
    }
    if(all_positions[locus].size()==ADVANCED_N){ //iHH0[locus]==0
        iHH0[locus] = 1;
    }
    if(all_positions[locus].size()==1){
        if(locus -= 0 or locus == ADVANCED_D-1){
            iHH1[locus] = 0.5;
        }else{
            iHH1[locus] = 1;
        }
    }
    if(all_positions[locus].size()==ADVANCED_N-1){
        if(locus == 0 or locus == ADVANCED_D-1){
            iHH0[locus] = 0.5;
        }else{
            iHH0[locus] = 1;
        }
    }
}

void EHH::exec() {
	this->init(); //initialize
    cout<<"Number of valid loci: "<<ADVANCED_D<<endl;
    cout<<"Number of haplotypes: "<<ADVANCED_N<<endl;
    if(calc_all){
        calc_iHS();
    }else{
        calc_EHH(locus);
    }
}

//IDEA: memoize
unsigned int twice_num_pair(int n){
    return n*n - n;
}

void EHH::calc_EHH2(int locus, map<int, vector<int> > & m, bool downstream){
    int ehh0_before_norm = 0;
    int ehh1_before_norm = 0;

    int curr_ehh0_before_norm = 0;
    int curr_ehh1_before_norm = 0;

    int n_c0=0;
    int n_c1=0;
    int n_c0_squared_minus = 0;
    int n_c1_squared_minus = 0;

    int group_count[ADVANCED_N];
    int group_id[ADVANCED_N];
    bool isDerived[ADVANCED_N]; 

    //will be vectorized
    for(int i = 0; i<ADVANCED_N; i++){
        group_count[i] = 0;
        group_id[i] = 0;
        isDerived[i] = false;
    }

    int totgc=0;
    vector<unsigned int> v = all_positions[locus];

    if(v.size()==0){
        n_c0 = ADVANCED_N;

        group_count[0] = ADVANCED_N;
        totgc+=1;

        ehh0_before_norm = n_c0*n_c0 - n_c0;
    }else if (v.size()==ADVANCED_N){ // all set
        group_count[0] = ADVANCED_N;
        totgc+=1;
        n_c1 = ADVANCED_N;
        
        for (int set_bit_pos : v){
            isDerived[set_bit_pos] = true;
        }
        ehh1_before_norm = n_c1*n_c1 - n_c1;

    }else{
        group_count[1] = v.size();
        group_count[0] = ADVANCED_N - v.size();
        n_c1 = v.size();
        n_c0 = ADVANCED_N - v.size();

        for (int set_bit_pos : v){
            isDerived[set_bit_pos] = true;
            group_id[set_bit_pos] = 1;
        }
        
        totgc+=2;
        ehh0_before_norm = n_c0*n_c0 - n_c0;
        ehh1_before_norm = n_c1*n_c1 - n_c1;
    }

    if(downstream){
        if(!calc_all)
            cout<<"Iter "<<0<<": EHH1["<<locus<<","<<locus<<"]="<<1<<" "<<1<<endl;
        
        if(n_c1_squared_minus!=0){
            iHH1[locus] += (curr_ehh1_before_norm + ehh1_before_norm) * 0.5/n_c1_squared_minus;
            //cout<<"Summing "<<1.0*curr_ehh1_before_norm/n_c1_squared_minus<<" "<<1.0*ehh1_before_norm/n_c1_squared_minus<<endl;
        }
        if(n_c0_squared_minus!=0){
            iHH0[locus] += (curr_ehh0_before_norm + ehh0_before_norm) * 0.5/n_c0_squared_minus;
        }
    }
    
    curr_ehh1_before_norm = ehh1_before_norm;
    curr_ehh0_before_norm = ehh0_before_norm;
    
    n_c1_squared_minus =  n_c1* n_c1 -  n_c1;
    n_c0_squared_minus =  n_c0* n_c0 -  n_c0;
    
    int i = locus;  
    while(true){ // Upstream: for ( int i = locus+1; i<all_positions.size(); i++ )
        if(downstream){
            if (--i < 0) break;
            if (loc_map[locus]-i > max_extend) break;
        }else{
            if (++i >= ADVANCED_D) break;
            if (i-locus > max_extend) break;
        }
        if(curr_ehh1_before_norm*1.0/n_c1_squared_minus < cutoff and curr_ehh0_before_norm*1.0/n_c0_squared_minus < cutoff){
            break;
        }
        int distance;
        if(downstream){
            distance = loc_map[i+1] - loc_map[i];
        }else{
            distance = loc_map[i] - loc_map[i-1];
        }
        //distance = 1;
        
        //OPT IDEA: INSTEAD OF MAP JUST USE ARR OF VECTOR
        v = all_positions[i];
        if(all_positions[i].size()==0 or all_positions[i].size()==ADVANCED_N){
            //cout<<"SKIPPING MONOMORPHIC SITE"<<endl;
            // if(!calc_all)
            //     cout<<"Mono: Iter "<<i-locus<<": EHH1["<<locus<<","<<i<<"]="<<curr_ehh1_before_norm*1.0/n_c1_squared_minus<<","<<curr_ehh0_before_norm*1.0/n_c0_squared_minus<<endl;
        
            
            if(n_c1_squared_minus!=0){
                iHH1[locus] += (curr_ehh1_before_norm + ehh1_before_norm) * distance * 0.5/n_c1_squared_minus;
                //cout<<"Summing "<<1.0*curr_ehh1_before_norm/n_c1_squared_minus<<" "<<1.0*ehh1_before_norm/n_c1_squared_minus<<endl;
            }
            if(n_c0_squared_minus!=0){
                iHH0[locus] += (curr_ehh0_before_norm + ehh0_before_norm) * distance * 0.5/n_c0_squared_minus;
            }
            continue;
        }

        for (int set_bit_pos : v){
            int old_group_id = group_id[set_bit_pos];
            m[old_group_id].push_back(set_bit_pos);
        }

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
            group_count[old_group_id] -= newgroup_size;
            group_count[totgc] += newgroup_size;
            totgc+=1;

            bool isDerivedGroup =  isDerived[ele.second[0]]; // just check first element to know if it is derived. 
            if(isDerivedGroup)
            {
                ehh1_before_norm += del_update;
            }else{
                ehh0_before_norm += del_update;
            }
            
        }

        if(n_c1_squared_minus!=0){
            iHH1[locus] += (curr_ehh1_before_norm + ehh1_before_norm) * distance * 0.5/n_c1_squared_minus;
            //cout<<"Summing "<<1.0*curr_ehh1_before_norm/n_c1_squared_minus<<" "<<1.0*ehh1_before_norm/n_c1_squared_minus<<endl;
        }

        if(n_c0_squared_minus!=0){
            iHH0[locus] += (curr_ehh0_before_norm + ehh0_before_norm) * distance * 0.5/n_c0_squared_minus;
        }

        curr_ehh1_before_norm = ehh1_before_norm;
        curr_ehh0_before_norm = ehh0_before_norm;

        // ehh1[locus] = 1.0*ehh1_before_norm/n_c1_squared_minus;
        // ehh0[locus] = 1.0*ehh0_before_norm/n_c0_squared_minus;

        //this shouldn't execute
        if(n_c1_squared_minus==0){
            curr_ehh1_before_norm = 0;
        }
        if(n_c0_squared_minus==0){
            curr_ehh0_before_norm = 0;
        }


        if(!calc_all)
            cout<<"Iter "<<i-locus<<": EHH1["<<locus<<","<<i<<"]="<<curr_ehh1_before_norm*1.0/n_c1_squared_minus<<","<<curr_ehh0_before_norm*1.0/n_c0_squared_minus<<endl;
        
        //logg[tid]+="map size before"+to_string(m.size())+"\n";
        //cout<< logg[tid];
        m.clear();
        //logg[tid]+="map size after"+to_string(m.size())+"\n";
        //cout<< logg[tid];
    }
}



void EHH::thread_ihs(int tid, map<int, vector<int> >& m, map<int, vector<int> >& md, EHH* ehh_obj){
    int elem_per_block = floor(ehh_obj->ADVANCED_D/ehh_obj->numThread);
    int start = tid*elem_per_block ;
    int end = start + elem_per_block  ;
    if(tid == ehh_obj->numThread-1 ){
        end = ehh_obj->ADVANCED_D;
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

    bool thread_enabled = true;
    if (thread_enabled)
    {
        int total_calc_to_be_done = ADVANCED_D;
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
        #pragma clang loop unroll_count(8) // 
        #pragma clang loop vectorize(assume_safety)
        for(int i = 0 ; i< ADVANCED_D; i++){
            calc_EHH(i);
        }
    }
   
    out_fp = fopen(output_filename.c_str(),"w");
    for (int i = 0; i < ADVANCED_D; i++){
        fprintf(out_fp, "%d %d %f %f %f %f\n", loc_map[i], loc_map[i], all_positions[i].size()*1.0/ADVANCED_N, iHH1[i], iHH0[i], log10(iHH1[i]/ iHH0[i]));
    }
    fclose(out_fp);
   
    return;
}


