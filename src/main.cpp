// main.cpp
// version 1:
// sep 26
// upstream only
// dont store it at all
// just store EHH values 
#include <thread>                            
#include <iostream>
#include <algorithm>

#include <fstream>
//#include "selscan.h"

#include<set>
#include "map"
#include <cmath>

#include <memory>
#include <iostream>
#include <sstream>
#include <vector>

#include "cxxopts.hpp"


using namespace std;

// #define ADVANCED_N 4
// #define ADVANCED_D 6404


// #define FILENAME "data/out5000.imp#define FILENAME "data/out5000.impute.hap"
// #define ADVANCED_N 6404
// #define ADVANCED_D 5000ute.hap"
// #define ADVANCED_N 6404
// #define ADVANCED_D 5000


// #define FILENAME "/storage/home/aur1111/s/dataset/out20000.impute.hap"
// #define ADVANCED_N 6404
// #define ADVANCED_D 20000

// #define FILENAME "/storage/home/aur1111/s/msdir/ms5k6k.hap2"
// #define ADVANCED_N 6000
// #define ADVANCED_D 5000


// #define FILENAME "data/out500.impute.hap"
// #define ADVANCED_N 6404
// #define ADVANCED_D 500

// #define ADVANCED_N 5
// #define ADVANCED_D 5
// #define FILENAME "test.hap"

#define ADVANCED_N 4
#define ADVANCED_D 6
#define FILENAME "test4x6.hap"


// #define ADVANCED_N 3
// #define ADVANCED_D 6404

// #define ADVANCED_N 6404
// #define ADVANCED_D 20000

// #define ADVANCED_N 6404
// #define ADVANCED_D 5000
std::vector<std::vector<unsigned int> > all_positions;
float ehh0[ADVANCED_D];
float ehh1[ADVANCED_D];
float ehh0_downstream[ADVANCED_D];
float ehh1_downstream[ADVANCED_D];

float iHH0[ADVANCED_D];
float iHH1[ADVANCED_D];
float iHH0_downstream[ADVANCED_D];
float iHH1_downstream[ADVANCED_D];

long nsl0[ADVANCED_D];
long nsl1[ADVANCED_D];


double start_time=0;

#define NUM_THREAD 8
string logg[NUM_THREAD];

map<int, vector<int> > map_per_thread[NUM_THREAD];
map<int, vector<int> > mapd_per_thread[NUM_THREAD];


int N = 0;
int D = 0;

double readTimer() {
    return clock() / (double) CLOCKS_PER_SEC;
}
unsigned int countIntersectFromSortedLists(std::vector<unsigned>& a, std::vector<unsigned> & b){
    int counter1 = 0;
    int counter2 = 0;
    int n11 = 0;

    while (true){
        if(counter1 >= a.size() || counter2 >= b.size()){
            return n11;
        }
        if (a[counter1] == b[counter2]){
            ++n11;
            counter1 += 1;
            counter2 += 1;
        }else if (a[counter1] < b[counter2]){
            counter1 += 1;
        }else{
            counter2 += 1;
        }   
    }
}

template <typename T> void print_a(T* arr, string name="v"){
    cout<<"vector: "<<name<<": ";
    for (int i=0; i<ADVANCED_N ; i++){
        cout<<arr[i]<<" ";
    }
    cout<<endl;
}

void print_v(vector<unsigned int> v, string name="v"){
    cout<<"vector: "<<name<<": ";
    for (int i: v){
        cout<<i<<" ";
    }
    cout<<endl;
}

//std::vector<std::vector<unsigned int> > &
void readMatrixMethod2(const std::string& filename){
    std::ifstream file(filename);
    double MAF=0;
 
    if (!file.is_open()) {
        std::cerr << "Failed to open file: " << filename << std::endl;
        return;
    }

    std::string line;
    while (std::getline(file, line)) {
        cout<<line<<endl;
        std::vector<unsigned int> positions;
 
        std::istringstream rowStream(line);
        char value;
        int pos=0;
        while (rowStream >> value) {        
            if(value=='1'){
                positions.push_back(pos++);
            }else if(value=='0'){
                pos++;
            }   
        }

        if(pos!=ADVANCED_N){
            cout<<"ERROR";
            exit(2);
        }

        if(N==0){
            N = pos;
        }

    
        // for(int i=0; i<positions.size();i++){
        //     std::cout<<i<<" "<<arr[i]<<endl;
        // }
        if(positions.size()==0){
           //cout<<"Trigger error"<<endl;
           //exit(1);
        }
        all_positions.push_back(positions); //check if all 0
    }
    file.close();
}

void calc_EHH(map<int, vector<int> >& m, int locus=0, bool print=false){
    //calc IHH
    iHH0[locus] = 0;
    iHH1[locus] = 0;


    int ehh0_before_norm = 0;
    int ehh1_before_norm = 0;
    //ehh0[locus] = 0;
    //ehh1[locus] = 0;
    int n_c0= 0;
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
    vector<unsigned int> &v = all_positions[locus];

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
            group_id[set_bit_pos] = 1;
        }
        
        totgc+=2;

        ehh0_before_norm = n_c0*n_c0 - n_c0;

        ehh1_before_norm = n_c1*n_c1 - n_c1;

    }

    
    n_c1_squared_minus =  n_c1* n_c1 -  n_c1;
    n_c0_squared_minus =  n_c0* n_c0 -  n_c0;
    
    for ( int i = locus+1; i<all_positions.size(); i++ ){
         //OPT IDEA: INSTEAD OF MAP JUST USE ARR OF VECTOR
        v = all_positions[i];
        for (int set_bit_pos : v){
            //cout<<set_bit_pos<<" ";
            int old_group_id = group_id[set_bit_pos];
            m[old_group_id].push_back(set_bit_pos);
        }
        //cout<<endl;
        for (const auto &ele : m) {
            int old_group_id = ele.first;
            //cout<<"old id"<<old_group_id<<endl;
            int newgroup_size = ele.second.size() ;
            
            if(group_count[old_group_id] == newgroup_size || newgroup_size == 0){
                continue;
            }
            
            
            for(int v: ele.second){
                group_id[v] = totgc;
            }
            
            int del_update = -pow(group_count[old_group_id],2) + pow(newgroup_size,2) + pow(group_count[old_group_id] - newgroup_size,2);
            group_count[old_group_id] -= newgroup_size;
            group_count[totgc] += newgroup_size;
            totgc+=1;

            bool isDerivedGroup =  isDerived[ele.second[0]];
            if(isDerivedGroup)// just check first element if it is derived or not, 
            {
                ehh1_before_norm += del_update;
            }else{
                ehh0_before_norm += del_update;
            }
            
        }

        //cout<<"GCC:"<<i<<" "<<totgc<<endl;
        if(n_c1_squared_minus!=0){
        iHH1[locus] += (ehh1[locus] +  (1.0*ehh1_before_norm/n_c1_squared_minus))*0.5;
        }

        if(n_c0_squared_minus!=0){
        iHH0[locus] += (ehh0[locus] +  (1.0*ehh0_before_norm/n_c0_squared_minus))*0.5;
        }

        ehh1[locus] = 1.0*ehh1_before_norm/n_c1_squared_minus;
        ehh0[locus] = 1.0*ehh0_before_norm/n_c0_squared_minus;

                
        if(n_c1_squared_minus==0){
            ehh1[locus] = 0;
        }
        if(n_c0_squared_minus==0){
            ehh0[locus] = 0;
        }
        if(print)
            cout<<"Iter "<<i<<": EHH1["<<locus<<"]="<<ehh1[locus]<<","<<ehh0[locus]<<endl;
        m.clear();
    }
    
    // for(int i = 0; i< ADVANCED_N; i++){
    //     if(isDerived[i]){
    //         ehh1[i]/=n_c1_squared_minus;
    //     }else{
    //         ehh0[i]/=n_c0_squared_minus;
    //     }
    //     cout<<"EHH["<<i<<"]="<<group_count[i]<<endl;
        
    // }
    // for(int i = 0; i< ADVANCED_N; i++){
        //cout<<"GC["<<i<<"]="<<group_count[i]<<endl;
    //     cout<<"GID["<<i<<"]="<<group_id[i]<<endl;
    // }
}





void calc_EHH2(int locus, map<int, vector<int> > & m){
        iHH0[locus] = 0;
        iHH1[locus] = 0;
        ehh1[locus] = 0;
        ehh0[locus] = 0;


        int ehh0_before_norm = 0;
        int ehh1_before_norm = 0;
        //ehh0[locus] = 0;
        //ehh1[locus] = 0;
        int n_c0= 0;
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
                group_id[set_bit_pos] = 1;
            }
            
            totgc+=2;

            ehh0_before_norm = n_c0*n_c0 - n_c0;

            ehh1_before_norm = n_c1*n_c1 - n_c1;

        }

        
        n_c1_squared_minus =  n_c1* n_c1 -  n_c1;
        n_c0_squared_minus =  n_c0* n_c0 -  n_c0;
        
        for ( int i = locus+1; i<all_positions.size(); i++ ){
            //OPT IDEA: INSTEAD OF MAP JUST USE ARR OF VECTOR
            v = all_positions[i];
            for (int set_bit_pos : v){
                //cout<<set_bit_pos<<" ";
                int old_group_id = group_id[set_bit_pos];
                m[old_group_id].push_back(set_bit_pos);
            }
            //cout<<endl;
            for (const auto &ele : m) {
                int old_group_id = ele.first;
                //cout<<"old id"<<old_group_id<<endl;
                int newgroup_size = ele.second.size() ;
                
                if(group_count[old_group_id] == newgroup_size || newgroup_size == 0){
                    continue;
                }
                
                
                for(int v: ele.second){
                    group_id[v] = totgc;
                }
                
                int del_update = -pow(group_count[old_group_id],2) + pow(newgroup_size,2) + pow(group_count[old_group_id] - newgroup_size,2);
                group_count[old_group_id] -= newgroup_size;
                group_count[totgc] += newgroup_size;
                totgc+=1;

                bool isDerivedGroup =  isDerived[ele.second[0]];
                if(isDerivedGroup)// just check first element if it is derived or not, 
                {
                    //iHH1[locus]+=ehh1_before_norm;
                    ehh1_before_norm += del_update;
                    //iHH1[locus]+=ehh1_before_norm;
                    
                }else{
                    //iHH0[locus]+=ehh0_before_norm;
                    ehh0_before_norm += del_update;
                   // iHH0[locus]+=ehh0_before_norm;

                }
                
            }

            //cout<<"GCC:"<<i<<" "<<totgc<<endl;
            // iHH1[locus] += (ehh1[locus] +  (1.0*ehh1_before_norm/n_c1_squared_minus))*0.5;
            // iHH0[locus] += (ehh0[locus] +  (1.0*ehh0_before_norm/n_c0_squared_minus))*0.5;

            if(n_c1_squared_minus!=0){
            iHH1[locus] += (ehh1[locus] + ehh1_before_norm) * 0.5/n_c1_squared_minus;
            }
            if(n_c0_squared_minus!=0){
            iHH0[locus] += (ehh0[locus] + ehh0_before_norm) * 0.5/n_c0_squared_minus;
            }



            ehh1[locus] = ehh1_before_norm;
            ehh0[locus] = ehh0_before_norm;

            // ehh1[locus] = 1.0*ehh1_before_norm/n_c1_squared_minus;
            // ehh0[locus] = 1.0*ehh0_before_norm/n_c0_squared_minus;

            if(n_c1_squared_minus==0){
                ehh1[locus] = 0;
            }
            if(n_c0_squared_minus==0){
                ehh0[locus] = 0;
            }
            if(false)
                cout<<"Iter "<<i<<": EHH1["<<locus<<"]="<<ehh1[locus]*1.0/n_c1_squared_minus<<","<<ehh0[locus]*1.0/n_c0_squared_minus<<endl;
            
            //logg[tid]+="map size before"+to_string(m.size())+"\n";
            //cout<< logg[tid];
            m.clear();
            //logg[tid]+="map size after"+to_string(m.size())+"\n";
            //cout<< logg[tid];

        }
    }

    void calc_EHH_downstream(int locus, map<int, vector<int> > & m){
        iHH1_downstream[locus] = 0;
        iHH0_downstream[locus] = 0;
        ehh1_downstream[locus] = 0;
        ehh0_downstream[locus] = 0;

        if(locus == 0){
            ehh1_downstream[0] = 0;
            ehh0_downstream[0] = 0;
            return;
        }
        int ehh0_before_norm = 0;
        int ehh1_before_norm = 0;
        //ehh0[locus] = 0;
        //ehh1[locus] = 0;
        int n_c0= 0;
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
                group_id[set_bit_pos] = 1;
            }
            
            totgc+=2;

            ehh0_before_norm = n_c0*n_c0 - n_c0;

            ehh1_before_norm = n_c1*n_c1 - n_c1;

        }

        
        n_c1_squared_minus =  n_c1* n_c1 -  n_c1;
        n_c0_squared_minus =  n_c0* n_c0 -  n_c0;
        
        for ( int i = locus-1; i>=0; i-- ){
            //OPT IDEA: INSTEAD OF MAP JUST USE ARR OF VECTOR
            v = all_positions[i];
            for (int set_bit_pos : v){
                //cout<<set_bit_pos<<" ";
                int old_group_id = group_id[set_bit_pos];
                m[old_group_id].push_back(set_bit_pos);
            }
            //cout<<endl;
            for (const auto &ele : m) {
                int old_group_id = ele.first;
                //cout<<"old id"<<old_group_id<<endl;
                int newgroup_size = ele.second.size() ;
                
                if(group_count[old_group_id] == newgroup_size || newgroup_size == 0){
                    continue;
                }
                
                
                for(int v: ele.second){
                    group_id[v] = totgc;
                }
                
                int del_update = -pow(group_count[old_group_id],2) + pow(newgroup_size,2) + pow(group_count[old_group_id] - newgroup_size,2);
                group_count[old_group_id] -= newgroup_size;
                group_count[totgc] += newgroup_size;
                totgc+=1;

                bool isDerivedGroup =  isDerived[ele.second[0]];
                if(isDerivedGroup)// just check first element if it is derived or not, 
                {
                    //iHH1_downstream[locus] += ehh1_before_norm;
                    ehh1_before_norm += del_update;
                    //iHH1_downstream[locus] += ehh1_before_norm;
                    
                }else{
                    //iHH0_downstream[locus] += ehh0_before_norm;
                    ehh0_before_norm += del_update;
                    //iHH0_downstream[locus] += ehh0_before_norm;

                }
                
            }

            //cout<<"GCC:"<<i<<" "<<totgc<<endl;
            // iHH1_downstream[locus] += (ehh1_downstream[locus] +  (1.0*ehh1_before_norm/n_c1_squared_minus))*0.5;
            // iHH0_downstream[locus] += (ehh0_downstream[locus] +  (1.0*ehh0_before_norm/n_c0_squared_minus))*0.5;
            // iHH1_downstream[locus] *= 0.5/n_c1_squared_minus;
            // iHH0_downstream[locus] *= 0.5/n_c0_squared_minus;

            // ehh1_downstream[locus] = 1.0*ehh1_before_norm/n_c1_squared_minus;
            // ehh0_downstream[locus] = 1.0*ehh0_before_norm/n_c0_squared_minus;


            if(n_c1_squared_minus!=0){
            iHH1_downstream[locus] += (ehh1_downstream[locus] + ehh1_before_norm) * 0.5/n_c1_squared_minus;
            }
            if(n_c0_squared_minus!=0){
            iHH0_downstream[locus] += (ehh0_downstream[locus] + ehh0_before_norm) * 0.5/n_c0_squared_minus;
            }
            ehh1_downstream[locus] = ehh1_before_norm;
            ehh0_downstream[locus] = ehh0_before_norm;


            if(n_c1_squared_minus==0){
                ehh1_downstream[locus] = 0;
            }
            if(n_c0_squared_minus==0){
                ehh0_downstream[locus] = 0;
            }
            // if(print)
            //     cout<<"Iter "<<i<<": EHH1["<<locus<<"]="<<ehh1[locus]<<","<<ehh0[locus]<<endl;
            //
            //logg[tid]+="map size before"+to_string(m.size())+"\n";
            //cout<< logg[tid];
            m.clear();
            //logg[tid]+="map size after"+to_string(m.size())+"\n";
            //cout<< logg[tid];

        }
    }
void thread_ihs(int tid, map<int, vector<int> >& m, map<int, vector<int> >& md){
    int elem_per_block = floor(ADVANCED_D/NUM_THREAD);
    // int start = tid*elem_per_block + 1;
    // int end = start + elem_per_block -1 ;
    int start = tid*elem_per_block ;
    int end = start + elem_per_block  ;
    if(tid == NUM_THREAD-1 ){
        end = ADVANCED_D;
    }

    // for(int j = start; j< end; j++){
    //     //cout<<j<<" ";
    //     logg[tid]+=to_string(j)+" ";
    // }
    // logg[tid]+="\n";

    ///*
    //#pragma omp parallel 
    for(int j = start; j< end; j++){
        //#pragma omp task 
        calc_EHH2(j, m);
        //#pragma omp task 
        calc_EHH_downstream(j, md);
    }
    // for(int j = start; j< end; j++){
    //     calc_EHH_downstream(j, m);
    // }
    //*/thread_ihs
    
    logg[tid]+="finishing thread #"+to_string(tid)+"\n"; 
    
    cout<<"finishing thread # "+to_string(tid)+" at "+to_string(readTimer())+"\n";
}

float calc_iHS(){

    bool thread_enabled = true;
    if (thread_enabled)
    {
        int total_calc_to_be_done = ADVANCED_D;

        std::thread *myThreads = new std::thread[NUM_THREAD];
        for (int i = 0; i < NUM_THREAD; i++)
        {
            myThreads[i] = std::thread(thread_ihs, i, std::ref(map_per_thread[i]),  std::ref(mapd_per_thread[i]));
        }
        for (int i = 0; i < NUM_THREAD; i++)
        // Join will block our main thread, and so the program won't exit until
        // everyone comes home.
        {
            myThreads[i].join();
        }
        delete[] myThreads;
        cout << "all threads finished. now calculating ihh..." << endl;
    }else{
        map<int, vector<int> > m;
        #pragma clang loop unroll_count(8) // 
        #pragma clang loop vectorize(assume_safety)
        for(int i = 0 ; i< ADVANCED_D; i++){
            calc_EHH(m, i);
        }
    }

    //THREADED

   
  //sum all the ones with 0
  //sum all the ones with 1
    float ihh1=0;
    float ihh0=0;
    for (int i = 1; i < ADVANCED_D; i++){
        //cout<<i<<" "<<(ehh1[i-1] + ehh1[i]) << " " <<(ehh0[i-1] + ehh0[i])<<endl;
        // ihh1 += (ehh1[i-1]+ehh1_downstream[i-1] + ehh1[i]+ehh1_downstream[i])*0.5; 
        // ihh0 += (ehh0[i-1]+ehh0_downstream[i-1] + ehh0[i]+ehh0_downstream[i])*0.5; 
        ihh1 = iHH1[i] + iHH1_downstream[i] ;
        ihh0 = iHH0[i] + iHH0_downstream[i] ;
        cout<<i << " ihh1 "<<ihh1<<" ihh0 "<<ihh0<<endl;
    }
    //cout<<"ihh1, ihh0 = "<<ihh1<<" "<<ihh0<<endl;
    // for(int i =0; i< NUM_THREAD; i++){
    //     cout<<logg[i]<<endl;
    // }
    return log(ihh1/ihh0);
}


void calc_nSL(int locus, map<int, vector<int> > & m){

        //if group_count[location] == 1 already unique
        //stop when totgc = N
        nsl0[locus] = 0;
        nsl1[locus] = 0;
        
        int n_c0= 0;
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

        }else if (v.size()==ADVANCED_N){ // all set
            group_count[0] = ADVANCED_N;
            totgc+=1;
            n_c1 = ADVANCED_N;
            
            for (int set_bit_pos : v){
                isDerived[set_bit_pos] = true;
            }

        }else{
            group_count[1] = v.size();
            group_count[0] = ADVANCED_N - v.size();
            n_c1 = v.size();
            n_c0 = ADVANCED_N - v.size();

            for (int set_bit_pos : v){
                group_id[set_bit_pos] = 1;
                isDerived[set_bit_pos] = true;

            }
            
            totgc+=2;

        }

        
        n_c1_squared_minus =  n_c1* n_c1 -  n_c1;
        n_c0_squared_minus =  n_c0* n_c0 -  n_c0;
        
        bool left_extreme = false;
        bool right_extreme = false;
        if(locus == 0){
            left_extreme = true;
        }
        if(locus == ADVANCED_D-1 ){
            right_extreme = true;
        }

        print_a<bool>(isDerived, "derived");

        print_a(group_count, "GC");
        print_a(group_id, "GID");

        for ( int i = locus+1; i<ADVANCED_D; i++ ){
            
           
            int neg_i = 2 * locus - i;
           

            int n11 = 0;
            int n10 = 0;
            int n00 = 0;
            int n01 = 0;
            
            vector<unsigned int> ancestral_right1;
            vector<unsigned int> derived_right1;
            vector<unsigned int> ancestral_left1;
            vector<unsigned int> derived_left1;

            if(not right_extreme){
                for(unsigned int ii : all_positions[i]){
                    if(isDerived[ii]){
                        derived_right1.push_back(ii);
                    }else{
                        ancestral_right1.push_back(ii);
                    }
                }
            }


            
            if (not left_extreme){
                for(unsigned int ii :  all_positions[neg_i]){
                    if(isDerived[ii]){
                        derived_left1.push_back(ii);
                    }else{
                        ancestral_left1.push_back(ii);
                    }
                }
            }



            print_v(derived_left1, "DL1");
            print_v(derived_right1, "DR1");

            print_v(ancestral_left1, "AL1");
            print_v(ancestral_right1, "AR1");

            //ancestral
            
            if(left_extreme){
                n11 = ancestral_right1.size();
                n10 = n_c0 - n11;
                n01 = 0;
                n00 = 0;
                nsl0[locus] +=  n10*n11*(i);

            }else if(right_extreme){
                n11 = ancestral_left1.size();
                n01 = n_c0 - n11;
                n10 = 0;
                n00 = 0;
                nsl0[locus] +=  n01*n11*(i);
            }else{
                n11 = countIntersectFromSortedLists(ancestral_right1, ancestral_left1);
                n10 = ancestral_left1.size() - n11;
                n01 = ancestral_right1.size() - n11;
                n00 = n_c0 - ancestral_right1.size() - n10;
                nsl0[locus] += n00*(i+1)*(n10 + n01) + n00*(i)*n11 + n01*n10*(i) + n01*n11*(i+1) + n10*n11*(i+1);

            }
            
            std::cout<<"at "<<i<<" n11,n10,n01,n00="<<n11<<","<<n10<<","<<n01<<","<<n00<<endl;
            
           
            //derived
            if(left_extreme){
                n11 = derived_right1.size();
                n10 = n_c1 - n11;
                n01 = 0;
                n00 = 0;
                nsl1[locus] +=  n10*n11*(i);

            }else if(right_extreme){
                n11 = derived_left1.size();
                n01 = n_c1 - n11;
                n10 = 0;
                n00 = 0;
                nsl1[locus] +=  n01*n11*(i);

            }else{
                n11 = countIntersectFromSortedLists(derived_right1, derived_left1);
                n10 = derived_left1.size() - n11;
                n01 = derived_right1.size() - n11;
                n00 = n_c1 - derived_right1.size() - n10;
                nsl1[locus] += n00*(i+1)*(n10 + n01) + n00*(i)*n11 + n01*n10*(i) + n01*n11*(i+1) + n10*n11*(i+1);

            }
            std::cout<<"at "<<i<<" n11,n10,n01,n00="<<n11<<","<<n10<<","<<n01<<","<<n00<<endl;
            


            if(not right_extreme){
                for (int set_bit_pos : all_positions[i])
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

                    group_count[old_group_id] -= newgroup_size;
                    group_count[totgc] += newgroup_size;
                    totgc += 1;
                }
                m.clear();
            }
           

            if (not left_extreme){

                for (int set_bit_pos : all_positions[neg_i]){
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
                    
                    group_count[old_group_id] -= newgroup_size;
                    group_count[totgc] += newgroup_size;
                    totgc+=1;
                }
                m.clear();
            }
            

            if(n_c1_squared_minus==0){
                nsl1[locus] = 0;
            }
            if(n_c0_squared_minus==0){
                nsl0[locus] = 0;
            }
           
            
            if(left_extreme and right_extreme){
                set<int> derivedGroups;
                for(int i: all_positions[locus]){
                    derivedGroups.insert(group_id[i]);
                }
                //for all groups having a non-1 group count do nc2
                for(int k = 0; k<totgc; k++){
                    int count = group_count[k];
                    if(count > 1){ //TRY TO OPT
                        if(derivedGroups.count(k)){
                           nsl1[locus]+=0.5*( count*count - count);
                        }else{
                           nsl0[locus]+=0.5*( count*count - count);
                        }
                    }
                }
            }
             if(true)
                std::cout<<"Iter "<<i<<": nsl[1,0]["<<locus<<"]="<<nsl1[locus]*2.0/n_c1_squared_minus<<","<<nsl0[locus]*2.0/n_c0_squared_minus<<endl;
             if(true)
                std::cout<<"Iter "<<i<<": nsl[1,0]["<<locus<<"]="<<nsl1[locus]<<","<<nsl0[locus]<<endl;


            if (i>=all_positions.size()-1){
                right_extreme = true;
            }
            if (neg_i <= 0)
            {
                left_extreme = true;
            }

            if(totgc==ADVANCED_N){
                return;
            }
        }
    }



  int main(int argc, char **argv)
  {   
    start_time = readTimer();
    vector<string> args(argv + 1, argv + argc);


    // cxxopts::Options options("Selscan 3.0", "[TODO] Program to calculate haplotype homozygosity");
    // options.add_options()
    // ("d,debug", "Enable debugging") // a bool parameter
    // ("i,integer", "Int param", cxxopts::value<int>())
    // ("f,file", "File name", cxxopts::value<std::string>())
    // ("v,verbose", "Verbose output", cxxopts::value<bool>()->default_value("false"))
    // ;

    // auto result = options.parse(argc, argv);

    // int x = result["opt"].as<int>();
    // string xx = result["file"].as<string>();



    // string filename;
    // for (auto i = args.begin(); i != args.end(); ++i) {
    //     if (*i == "-h" || *i == "--help") {
    //         cout << "Syntax: tool -g debug -i <DE-DUP-bitmatrix> -d <dup-bitmatrix> -c <num-colors> -m <M> -k <num-kmers> -s <spss-bound> -x <max-run>" << endl;
    //         return 0;
    //     } else if (*i == "-i") {
    //         filename = *++i;
    //     } 
    // }
    
    string filename = FILENAME;
   // string filename = "data/out500.impute.hap";

    //string filename = "data/out20000.impute.hap";
    //string filename ="data/out3.txt";
    //string filename ="test.hap";

    //readMatrixFromFile("data/out20000.impute.hap");
    //readMatrixFromFile("data/out3.txt");
    // readMatrixFromFile("test.hap");
    
  readMatrixMethod2(filename); 
  //calc_EHH(4, true);
  //calc_EHH(5, true);
  //calc_EHH(6, true);

  //calc_EHH(7, true);


//   calc_EHH(7504);
    // float ihs = calc_iHS();
    // cout<<"iHS="<<ihs<<endl;

     map<int, vector<int> >  m;
    calc_nSL(0, m);

    cout<<"Finish time:"<<to_string(readTimer())<<endl;
    return 0;
}

              