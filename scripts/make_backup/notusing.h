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

// template <typename T> void print_a(T* arr, string name="v"){
//     cout<<"vector: "<<name<<": ";
//     for (int i=0; i<ADVANCED_N ; i++){
//         cout<<arr[i]<<" ";
//     }
//     cout<<endl;
// }

// void print_v(vector<unsigned int> v, string name="v"){
//     cout<<"vector: "<<name<<": ";
//     for (int i: v){
//         cout<<i<<" ";
//     }
//     cout<<endl;
// }

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