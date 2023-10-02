// main.cpp
// version 1:
// sep 26
// upstream only
// dont store it at all
// just store EHH values 
#include <thread>                            
#include <iostream>
#include <algorithm>
#include "bm.h"
#include <fstream>
#include "bmserial.h"
#include "bmundef.h" 
#include "bm.h"
#include "bmaggregator.h"
#include "map"
#include <cmath>

#include <memory>
#include <iostream>
#include <sstream>
#include <vector>


using namespace std;

// #define ADVANCED_N 4
// #define ADVANCED_D 6404

// #define ADVANCED_N 5
// #define ADVANCED_D 5

// #define ADVANCED_N 3
// #define ADVANCED_D 6404

#define ADVANCED_N 6404
#define ADVANCED_D 20000
std::vector<std::vector<unsigned int> > all_positions;
float ehh0[ADVANCED_D];
float ehh1[ADVANCED_D];



int NUM_THREAD = 8;
int N = 0;
int D = 0;

int group_count[ADVANCED_D];

void bvenumerate(bm::bvector<> bv, string msg=""){
     //bm::bvector<> bv = bvs[0];
     std::cout<<msg<<" ";
    std::cout<<"Enumerating...";
    bm::bvector<>::enumerator en = bv.first();
    bm::bvector<>::enumerator en_end = bv.end();

    while (en < en_end)
    {
        std::cout << *en << " ";
        ++en;  // pre-increment - fastest way to increment enumerator
    }
    cout<<endl;
    en = bv.first();
}
static
void print_bvector(const bm::bvector<>& bv)
{
    bm::bvector<>::size_type cnt = 0;
    bm::bvector<>::enumerator en = bv.first();
    for (; en.valid() && cnt < 10; ++en, ++cnt)
        cout << *en << ", ";
    if (cnt == 10)
        cout << " ...";
    cout << "(size = "<< bv.size() << ")" << endl;
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

void calc_EHH(int locus=0){
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
        map<int, vector<int> > m;
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
        ehh1[locus] = 1.0*ehh1_before_norm/n_c1_squared_minus;
        ehh0[locus] = 1.0*ehh0_before_norm/n_c0_squared_minus;

        cout<<"Iter "<<i<<": EHH1["<<locus<<"]="<<ehh1[locus]<<","<<ehh0[locus]<<endl;

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

void thread_ihs(int tid){
    int elem_per_block = floor(ADVANCED_D/NUM_THREAD);
    int start = tid*elem_per_block;
    int end = start + elem_per_block;
    if(tid == NUM_THREAD-1 ){
        end = NUM_THREAD-1;
    }
    for(int j = start; j< end; j++){
        calc_EHH(j);
    }
    cout<<"finishing thread #"<<tid<<endl;
}

float calc_iHS(){
    float ihh1=0;
    float ihh0=0;

    int total_calc_to_be_done = ADVANCED_D;

    std::thread * myThreads = new std::thread[NUM_THREAD];
     for (int i = 0; i < NUM_THREAD; i++)
    {
      myThreads[i] = std::thread(thread_ihs, i);
    }
  for (int i = 0; i < NUM_THREAD; i++)
      // Join will block our main thread, and so the program won't exit until
      // everyone comes home.
    {
      myThreads[i].join();
    }
  delete [] myThreads;
  cout<<"all threads finished. now calculating ihh..."<<endl;
  //sum all the ones with 0
  //sum all the ones with 1
    for (int i = 1; i < ADVANCED_D; i++){
        ihh1 += (ehh1[i-1] + ehh1[i])*1; 
        ihh0 += (ehh0[i-1] + ehh0[i])*1; 
    }

    return log(ihh1/ihh0);
}
void readMatrixHalfTime(const std::string& filename) {
    bm::serializer<bm::bvector<> >::buffer sbuf[ADVANCED_D]; // declare serialization buffers
    bm::operation_deserializer<bm::bvector<> > od;
    std::ifstream file(filename);
    BM_DECLARE_TEMP_BLOCK(tb)
    bm::serializer<bm::bvector<> > bvs(tb);
    bvs.set_compression_level(4);
 
    if (!file.is_open()) {
        std::cerr << "Failed to open file: " << filename << std::endl;
        return;
    }

    std::string line;
    while (std::getline(file, line)) {
        std::vector<unsigned int> positions;
        bm::bvector<>   bv;    
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

        unsigned int* arr = new unsigned int[positions.size()];
        copy(positions.begin(),positions.end(),arr);
        // for(int i=0; i<positions.size();i++){
        //     std::cout<<i<<" "<<arr[i]<<endl;
        // }
        if(positions.size()!=0){
            bv.import_sorted(arr, positions.size(), true);
        }
        
        delete[] arr;
        bm::bvector<>::statistics st;
        bv.optimize(tb, bm::bvector<>::opt_compress, &st);
            
        bvs.serialize(bv, sbuf[D++], &st); // perform serialization
    }

    file.close();
    //------------------------------------------------------------
    //done reading

    bm::bvector<> bv_results[ADVANCED_N];
    // std::vector<std::unique_ptr<bm::bvector<> > > vect; //vect[i] entry gets free when u move to vect[i+1]
    // for (unsigned i = 0; i < ADVANCED_N; ++i)
    // {
    //     //bv->set(ADVANCED_N-1);
    //     //bvenumerate(std::move(bv));
         
    //     //plz dont delete
    //     //bm::bvector<> bvnew(*vect[i].get());

    //     bm::bvector<> bvnew;

    //     bm::deserialize(bvnew, sbuf[i].data());

    //     bv_results[i] = bvnew;
    //     bvenumerate(bv_results[i], "bv_results");
    //     std::unique_ptr<bm::bvector<> > bv(&bvnew);    
    // // for (unsigned i = 0; i < vect.size(); ++i)
    // // {
    // //     vect[i] = (vect[i].get());
    // // }

    //     vect.push_back(std::move(bv));

    //     bvenumerate(*vect[i].get());
        
    //     int cnt = od.deserialize(*vect[i].get(), sbuf[i].buf(), bm::set_COUNT_AND);
    //     //bvenumerate(bv_results[x]);
    //     std::cout<<"cnt: "<<cnt<<endl;

    // }

    group_count[0]=2;//assume 0 and 1

    int cnt = od.deserialize(bv_results[0], sbuf[0].buf(), bm::set_COUNT_OR);
    if(cnt==0 || cnt==N){
        bm::deserialize(bv_results[0], sbuf[0].buf()); //actually deserailize
        group_count[0]=1;
        if(cnt == 0){
            bv_results[0].resize(ADVANCED_N);
            bv_results[0].invert();
        }
    }else{
         group_count[0]=2;
         //needs checking
         bv_results[1].bit_or(bv_results[0]);
         bv_results[0].resize(ADVANCED_N);
         bv_results[0].invert();
    }

    bm::bvector<> bvnew;
    for(int i = 1 ; i< ADVANCED_D; i++){
        group_count[i] = group_count[i-1];
        //bm:: 
       // bvt2.bit_or(bv1, bv2, bvect::opt_compress);
       // bv_T.bit_sub(bv_results[i-1], bv_B, bm::bvector<>::opt_compress); // 3
        for(int j = 0; j< group_count[i-1]; j++){
            
            //read from sbuf_result[j]
            //perform and or with sbuf
            //update sbuf_result
            // update group count
            //print count
            
            int countAND = od.deserialize(bv_results[j], sbuf[i].buf(), bm::set_COUNT_AND); //sure its not 0
            int countSUB = od.deserialize(bv_results[j], sbuf[i].buf(), bm::set_COUNT_SUB_AB); //sure its not 0

            if(countAND==0){
                od.deserialize(bv_results[j], sbuf[i].buf(), tb, bm::set_SUB); //sure its not 0
                // cout<<"bv"<<i<<","<<j<<endl;
                // bvenumerate(bv_results[ j]);
            }else if(countSUB==0){
                od.deserialize(bv_results[j], sbuf[i].buf(), tb, bm::set_AND); //sure its not 0
                // cout<<"bv"<<i<<","<<j<<endl;
                // bvenumerate(bv_results[ j]);
            }else{
                bv_results[ group_count[i] ].bit_or(bv_results[j]);
                od.deserialize(bv_results[ group_count[i] ], sbuf[i].buf(), tb,bm::set_SUB); //sure its not 0
                od.deserialize(bv_results[j], sbuf[i].buf(), tb, bm::set_AND); //sure its not 0
                
                // cout<<"bv"<<i<<","<<j<<endl;
                // bvenumerate(bv_results[ j]);
                // bvenumerate(bv_results[ group_count[i]]);

                group_count[i]++;
                //IDEA-OPT
            }
           
        }
        //cout<<"GC["<<i<<"]"<<group_count[i]<<endl;
    }
    for(int i = 0; i< ADVANCED_N; i++){
        cout<<"GC["<<i<<"]="<<group_count[i]<<endl;
    }
}


void readMatrixFromFile(const std::string& filename) {
    bm::serializer<bm::bvector<> >::buffer sbuf[ADVANCED_D]; // declare serialization buffers
    
    bm::bvector<> allbvs[ADVANCED_D]; // declare serialization buffers
    
    bm::operation_deserializer<bm::bvector<> > od;
    std::ifstream file(filename);
    BM_DECLARE_TEMP_BLOCK(tb)
    bm::serializer<bm::bvector<> > bvs(tb);
    bvs.set_compression_level(4);
 
    if (!file.is_open()) {
        std::cerr << "Failed to open file: " << filename << std::endl;
        return;
    }

    std::string line;
    while (std::getline(file, line)) {
        std::vector<unsigned int> positions;
        //bm::bvector<>   bv;    
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

        unsigned int* arr = new unsigned int[positions.size()];
        copy(positions.begin(),positions.end(),arr);
        // for(int i=0; i<positions.size();i++){
        //     std::cout<<i<<" "<<arr[i]<<endl;
        // }
        if(positions.size()!=0){
            allbvs[D++].import_sorted(arr, positions.size(), true);
        }
        //bvenumerate(allbvs[D-1]);
        
        delete[] arr;
        //bm::bvector<>::statistics st;
        //bv.optimize(tb, bm::bvector<>::opt_compress, &st);
            
        //bvs.serialize(bv, sbuf[D++], &st); // perform serialization
    }


    file.close();

    //------------------------------------------------------------
    //done reading

    bm::bvector<> bv_results[ADVANCED_N];
    // std::vector<std::unique_ptr<bm::bvector<> > > vect; //vect[i] entry gets free when u move to vect[i+1]
    // for (unsigned i = 0; i < ADVANCED_N; ++i)
    // {
    //     //bv->set(ADVANCED_N-1);
    //     //bvenumerate(std::move(bv));
         
    //     //plz dont delete
    //     //bm::bvector<> bvnew(*vect[i].get());

    //     bm::bvector<> bvnew;

    //     bm::deserialize(bvnew, sbuf[i].data());

    //     bv_results[i] = bvnew;
    //     bvenumerate(bv_results[i], "bv_results");
    //     std::unique_ptr<bm::bvector<> > bv(&bvnew);    
    // // for (unsigned i = 0; i < vect.size(); ++i)
    // // {
    // //     vect[i] = (vect[i].get());
    // // }

    //     vect.push_back(std::move(bv));

    //     bvenumerate(*vect[i].get());
        
    //     int cnt = od.deserialize(*vect[i].get(), sbuf[i].buf(), bm::set_COUNT_AND);
    //     //bvenumerate(bv_results[x]);
    //     std::cout<<"cnt: "<<cnt<<endl;

    // }

    group_count[0]=2;//assume 0 and 1

    //int cnt = od.deserialize(bv_results[0], sbuf[0].buf(), bm::set_COUNT_OR);
    bm::bvector<> bv_T;
    bv_T.bit_or(bv_results[0], allbvs[0]);
    int cnt =  bv_T.count();

    if(cnt==0 || cnt==N){
        //bm::deserialize(bv_results[0], sbuf[0].buf()); //actually deserailize
        bv_results[0] =  allbvs[0];
        group_count[0]=1;
        if(cnt == 0){
            bv_results[0].resize(ADVANCED_N);
            bv_results[0].invert();
        }
    }else{
         group_count[0]=2;
         //needs checking
         bv_results[1].bit_or(bv_results[0]);
         bv_results[0].resize(ADVANCED_N);
         bv_results[0].invert();
    }

    bm::bvector<> bvnew;
    for(int i = 1 ; i< ADVANCED_D; i++){
        group_count[i] = group_count[i-1];
        //bm:: 
       // bvt2.bit_or(bv1, bv2, bvect::opt_compress);
       // bv_T.bit_sub(bv_results[i-1], bv_B, bm::bvector<>::opt_compress); // 3
        for(int j = 0; j< group_count[i-1]; j++){
            
            //read from sbuf_result[j]
            //perform and or with sbuf
            //update sbuf_result
            // update group count
            //print count
            
            ////int countAND = od.deserialize(bv_results[j], sbuf[i].buf(), bm::set_COUNT_AND); //sure its not 0
            ////int countSUB = od.deserialize(bv_results[j], sbuf[i].buf(), bm::set_COUNT_SUB_AB); //sure its not 0

            bm::bvector<> bvAND;
            bm::bvector<> bvSUB;

            bvAND.bit_and(bv_results[j], allbvs[i]);
            bvSUB.bit_sub(bv_results[j], allbvs[i]);

            int countAND =  bvAND.count();
            int countSUB =  bvSUB.count();



            if(countAND==0){
                bv_results[j] = bvSUB;
                ////od.deserialize(bv_results[j], sbuf[i].buf(), tb, bm::set_SUB); //sure its not 0
                // cout<<"bv"<<i<<","<<j<<endl;
                // bvenumerate(bv_results[ j]);
            }else if(countSUB==0){
                bv_results[j] = bvAND;
                
                //od.deserialize(bv_results[j], sbuf[i].buf(), tb, bm::set_AND); //sure its not 0
                // cout<<"bv"<<i<<","<<j<<endl;
                // bvenumerate(bv_results[ j]);
            }else{
                bv_results[ group_count[i] ] = bvSUB;
                bv_results[ j ] = bvAND;


                ////bv_results[ group_count[i] ].bit_or(bv_results[j]);
                ////od.deserialize(bv_results[ group_count[i] ], sbuf[i].buf(), tb,bm::set_SUB); //sure its not 0
                ////od.deserialize(bv_results[j], sbuf[i].buf(), tb, bm::set_AND); //sure its not 0
                
                // cout<<"bv"<<i<<","<<j<<endl;
                // bvenumerate(bv_results[ j]);
                // bvenumerate(bv_results[ group_count[i]]);

                group_count[i]++;
                //IDEA-OPT
            }
           
        }
        //cout<<"GC["<<i<<"]"<<group_count[i]<<endl;
    }
    for(int i = 0; i< ADVANCED_N; i++){
        cout<<"GC["<<i<<"]="<<group_count[i]<<endl;
    }
}


  int main(int argc, char **argv)
  {   
    vector<string> args(argv + 1, argv + argc);
    // string filename;
    // for (auto i = args.begin(); i != args.end(); ++i) {
    //     if (*i == "-h" || *i == "--help") {
    //         cout << "Syntax: tool -g debug -i <DE-DUP-bitmatrix> -d <dup-bitmatrix> -c <num-colors> -m <M> -k <num-kmers> -s <spss-bound> -x <max-run>" << endl;
    //         return 0;
    //     } else if (*i == "-i") {
    //         filename = *++i;
    //     } 
    // }
    
    string filename = "data/out20000.impute.hap";
    //string filename ="data/out3.txt";
    //string filename ="test.hap";

    //readMatrixFromFile("data/out20000.impute.hap");
    //readMatrixFromFile("data/out3.txt");
    // readMatrixFromFile("test.hap");
    
  readMatrixMethod2(filename); 
  calc_EHH(0);
  //cout<<"iHS="<<calc_iHS();


    return 0;
}

              