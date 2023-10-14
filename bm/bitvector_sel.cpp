#include "bm.h"
#include "bmserial.h"
#include "bmundef.h" 
#include "bm.h"
#include "bmaggregator.h"


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


void readMatrixHalfTime(const std::string& filename) {
    int group_count[ADVANCED_D];
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
    int group_count[ADVANCED_D];
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