#include <vector>
#include <string>
#include "feature.hpp"
#include<thread>
#include<queue>
#include "hapmap.hpp"
#include "utils.hpp"
using namespace std;

Feature::Feature () {
}

Feature::~Feature () {

}

void Feature::cli_prepare(CLI::App * app) {
	this->subapp = app->add_subcommand("feature", "Trying a new feature.");
    // CLI::Option * opt;
	// opt = subapp->add_option("-f", input_file_mphb, "A text file with columns: Left, Right, Seq Id, Num Seq")->capture_default_str();
	// opt->check(CLI::ExistingFile);
    // opt = subapp->add_option("--vcf", input_filename_vcf, "A VCF file with one row per haplotype, \n"
    //                                                             "and one column per variant.\n"
    //                                                             "Variants should be coded 0/1")->capture_default_str();
    // opt->check(CLI::ExistingFile);
	app->get_formatter()->column_width(10);
}


void loadMPHB(string input_mphb_file = "example/final_mphb_unsorted", int numHaps=6404, int numSnps=2000){
    vector< map<int, priority_queue<pair<int, int>  > > > mphbs(numHaps);
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
        mphbs[seqId][left].push(make_pair(numSeq, right));

        //cout<<seqId<<" "<<left<<" "<<right<<" "<<numSeq<<endl;
        //verify that it is larger than stack top
    }

    for (int i = 0 ; i < numHaps; i++){
        for (int j = 0 ; j < numSnps; j++){
            if(mphbs[i].count(j)){
                priority_queue<pair<int, int> > q = mphbs[i][j];
                //if (q.size()!=0 )
                    cout<<i << " "<< j << " "<< q.size()<<": ";
                    while (!q.empty()){
                        int right = q.top().first;
                        //int right = q.front().first;
                        q.pop();
                        cout<<right<<" ";
                    }
                    
                    cout<<endl;
            }
            
        }
    }

}


void Feature::init() {  
    loadMPHB();
}

void Feature::exec() {
	this->init(); //initialize
}
