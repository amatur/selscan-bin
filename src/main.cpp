// main.cpp
// version 1:
// sep 26
// upstream only
// dont store it at all
// just store EHH values 
                            
#include <iostream>
#include <algorithm>
#include "bm.h"
#include <fstream>
#include "bmserial.h"
#include "bmundef.h" 
#include "bm.h"
#include "bmaggregator.h"
using namespace std;

#include <memory>
#include <iostream>
#include <sstream>
#include <vector>

// #define ADVANCED_N 4
// #define ADVANCED_D 6404

#define ADVANCED_N 5
#define ADVANCED_D 5

int N = 0;
int D = 0;

int group_count[ADVANCED_N];
int ehh_0_i[ADVANCED_N];


// static
// void make_BLOB(vector<unsigned char>& target_buf, bm::bvector<>& bv)
// {
//     BM_DECLARE_TEMP_BLOCK(tb)
//     bm::serializer<bm::bvector<> > bvs(tb);
//     bvs.set_compression_level(4);
    
//     bv.optimize(tb, bm::bvector<>::opt_compress); // memory compression
 
//     bm::serializer<bm::bvector<> >::buffer sbuf;
//     bvs.serialize(bv, sbuf, 0);
//     target_buf.resize(sbuf.size());
//     ::memcpy(target_buf.data(), sbuf.buf(), sbuf.size());
// }
 
// static
// void make_BLOB(vector<unsigned char>& target_buf, bm::bvector<>& bv)
// {
//     BM_DECLARE_TEMP_BLOCK(tb)
//     bm::serializer<bm::bvector<> > bvs(tb);
//     bvs.set_compression_level(4);
    
//     bv.optimize(tb, bm::bvector<>::opt_compress); // memory compression
 
//     bm::serializer<bm::bvector<> >::buffer sbuf;
//     bvs.serialize(bv, sbuf, 0);
//     target_buf.resize(sbuf.size());
//     ::memcpy(target_buf.data(), sbuf.buf(), sbuf.size());
// }
 



// void testSerialize(){
//     // Prepare a serializer class
//         //  for best performance - create serilizer once and reuse it
//         //
//         BM_DECLARE_TEMP_BLOCK(tb)
//         bm::serializer<bm::bvector<> > bvs(tb);
//         bvs.set_compression_level(4);
 
        
//         for(bv in bvs){
//              // compress bit-vectors and compute statistics
//             // (later re-used in serialization)
//             //
//             bm::bvector<>::statistics st1;
//             bv.optimize(tb, bm::bvector<>::opt_compress, &st1);
            
//             // declare serialization buffers
//             bm::serializer<bm::bvector<> >::buffer sbuf1;
            
//             // perform serialization
//             //
//             bvs.serialize(bv1, sbuf1, &st1);

//             // Serialized bvectors (sbuf1 and sbuf2) now ready to be
//             // saved to a database, file or send over a network.
//             // to simulate this we just copy content to std::vector<>
//             //
//             /* uncomment to check that it is correct*/
//             // std::vector<unsigned char> vect1;
            
//             // vect1.resize(sbuf1.size());
            
//             // ::memcpy(vect1.data(), sbuf1.buf(), sbuf1.size());
    
//         }
        
       
 
 
// }

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

void readMatrixFromFile(const std::string& filename) {
    // declare serialization buffers
    bm::serializer<bm::bvector<> >::buffer sbuf[ADVANCED_N];
        
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

        if(D==0){
            D = pos;
        }

        unsigned int* arr = new unsigned int[positions.size()];
        copy(positions.begin(),positions.end(),arr);
        for(int i=0; i<positions.size();i++){
            std::cout<<i<<" "<<arr[i]<<endl;
        }
        bv.import_sorted(arr, positions.size(), true);
        delete[] arr;
        bm::bvector<>::statistics st;
        bv.optimize(tb, bm::bvector<>::opt_compress, &st);
            
          
        // perform serialization
        bvs.serialize(bv, sbuf[N++], &st);

        //adjacencyMatrix.push_back(bv);
        //break;
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
    //add to saved_group_count_j and group_count_j for each j
    //gets freebv_A
    // for(int x=0; x<N; x++){
    //     //int cnt = od.deserialize(*vect[x].get(), sbuf[x].buf(), bm::set_COUNT_AND);
    //     //bvenumerate(bv_results[x]);
    //     //std::cout<<"cnt: "<<cnt<<endl;
    //     bvenumerate(bv_results[x], "getsfree");
    //     bm::bvector<>::statistics st;
    //     bv_results[x].optimize(tb, bm::bvector<>::opt_compress, &st);
 
    // }   


   
    int cnt = od.deserialize(bv_results[0], sbuf[0].buf(), bm::set_COUNT_OR);
    //int cnt = od.deserialize(bv_results[1], sbuf[0].buf(), bm::set_COUNT_OR);
    cout<<"heck cnt "<<cnt<<endl;
    if(cnt==0 || cnt==N){
        bm::deserialize(bv_results[0], sbuf[0].buf()); //actually deserailize

        group_count[0]=1;
        if(cnt == 0){
            
            bv_results[0].resize(ADVANCED_N);
            bv_results[0].invert();
        }
        cout<<"CT0N: "<<bv_results[0].count()<<endl;
    }else{
         group_count[0]=2;
         //needs checking
         bv_results[1].bit_or(bv_results[0]);
         bv_results[0].resize(ADVANCED_N);
         bv_results[0].invert();
        cout<<"CT: "<<bv_results[1].count()<<endl;
    }


    bm::bvector<> bvnew;
    for(int i = 1 ; i< ADVANCED_N; i++){
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
            
            // bm::bvector<> bv_result_copy;
            // bv_result_copy.bit_or(bv_results[j]) ;
            
            //
            
            
            //bvenumerate(bv_results[j],"BEFORE");
            int countAND = od.deserialize(bv_results[j], sbuf[i].buf(), bm::set_COUNT_AND); //sure its not 0
            int countSUB = od.deserialize(bv_results[j], sbuf[i].buf(), bm::set_COUNT_SUB_AB); //sure its not 0
            //bvenumerate(bv_results[j],"AFTER");
            
            //bv_result.optimize(tb2);
            // //bv_result_copy.optimize(tb2);
            //
            
            // bm::bvector<> bv_result_copy2(bv_result) ;
            // bvenumerate(bv_result_copy2);
               // bvs2.serialize(bv_result_copy2, sbuf_results[i], &sts[i]);
            

            // if(countAND==0){
            //     //bv_results[j] = bv_result_copy;
            //     //bm::bv_results[j].AND(group_count[i-1]);
            //     //group_count[i] = group_count[i-1];
            // }else if(countSUB==0){
            //     //group_count[i] = group_count[i-1];
            //     //bv_results[group_count[i]] = 
            // }else{
            //     group_count[i]++;
            // }

            if(countAND==0){
                //bv_results[j] = bv_result_copy;
                //bm::bv_results[j].AND(group_count[i-1]);
                //group_count[i] = group_count[i-1];
                //bv_results[j].resize(ADVANCED_N);
                cout<<"countSUB0:"<<i<<","<<j<< "="<<od.deserialize(bv_results[j], sbuf[i].buf(), tb, bm::set_SUB)<<endl; //sure its not 0
                
                //bm::deserialize(bvnew, sbuf[i].buf()); //actually deserailize
                //bv_results[j].and_bit_no_check(bvnew);
                cout<<"bv"<<i<<","<<j<<endl;
                bvenumerate(bv_results[ j]);

                
            }else if(countSUB==0){
                //group_count[i] = group_count[i-1];
                //bv_results[group_count[i]] = 
                //bv_results[j].resize(ADVANCED_N);

                cout<<"countAND0: "<<i<<","<<j<<"="<<od.deserialize(bv_results[j], sbuf[i].buf(), tb, bm::set_AND)<<endl; //sure its not 0


                cout<<"bv"<<i<<","<<j<<endl;
                bvenumerate(bv_results[ j]);
                //bm::deserialize(bvnew, sbuf[i].buf()); //actually deserailize
                //bv_results[j].bit_and(bvnew);
            }else{
                //bv_results[j].resize(0);

                bv_results[ group_count[i] ].bit_or(bv_results[j]);
                bvenumerate(bv_results[ group_count[i] ],"BEFORE");

                //bv_results[ group_count[i] ].resize(ADVANCED_N);

               
                cout<<"countSUBAB: "<<i<<","<<j<< "="<<od.deserialize(bv_results[ group_count[i] ], sbuf[i].buf(), tb,bm::set_SUB)<<endl; //sure its not 0
                cout<<"countAND: "<<i<<","<<j<< "="<< od.deserialize(bv_results[j], sbuf[i].buf(), tb, bm::set_AND)<<endl; //sure its not 0
                
      
                cout<<"bv"<<i<<","<<j<<endl;
            bvenumerate(bv_results[ j]);
            bvenumerate(bv_results[ group_count[i]]);

                group_count[i]++;
                
                //IDEA-OPT
            }
            //break;            
            
        }
        //break;
        //cout << i<< ": bv count = " << bv_deserialized.count() << endl;
        cout<<"GC["<<i<<"]"<<group_count[i]<<endl;
        //bvenumerate(bv_deserialized);
    }
    for(int i = 0; i< ADVANCED_N; i++){
        cout<<"GC["<<i<<"]="<<group_count[i]<<endl;
    }
}


  int main(int argc, char **argv)
  {   
    vector<string> args(argv + 1, argv + argc);
    string filename;
    for (auto i = args.begin(); i != args.end(); ++i) {
        if (*i == "-h" || *i == "--help") {
            cout << "Syntax: tool -g debug -i <DE-DUP-bitmatrix> -d <dup-bitmatrix> -c <num-colors> -m <M> -k <num-kmers> -s <spss-bound> -x <max-run>" << endl;
            return 0;
        } else if (*i == "-i") {
            filename = *++i;
        } 
        //else if (*i == "-d") {
        //     dup_bitmatrix_fname = *++i;
        // }else if (*i == "-c") {
        //     C = std::stoi(*++i);
        // }else if (*i == "-m") {
        //     M = std::stoi(*++i);
        // }
    }
    
    //readMatrixFromFile("data/out20000.impute.hap");
    //readMatrixFromFile("data/out3.txt");
    readMatrixFromFile("test.hap");
    
    
    return 0;
}

              