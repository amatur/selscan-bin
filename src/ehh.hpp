#include <string>
#include <iostream>
#include <vector>
#include <map>
#include <thread>
#include <omp.h>
#include<mutex>
#include "logger.hpp"


#include<queue>
using namespace std;

#include "CLI11.hpp"
#include "selscan.hpp"
#include "hapmap.hpp"

#ifndef EHH_H
#define EHH_H

class EHH: public Selscan {

public:
	EHH();
	void cli_prepare(CLI::App * subapp);
	void init();
	void exec();
	~EHH();

	std::string input_filename_hap = "data/out20000.impute.hap";
	//hapbin-1.3.0/build/ihsbin  --map selscan-bin/data/out.20k.simple.map --hap selscan-bin/data/out20000.impute.hap -s 100

//	std::string input_filename_hap = "data/out500.impute.hap";
    //std::string input_filename_hap = "test4x6.hap";
    std::string input_filename_map = "test4x6.map";
    //std::string output_filename = "selbin.ihh.out";// <outfile>.nsl[.alt].out
    
	
	std::string output_filename_ehh = "selbin.ehh.out";
	std::string output_filename_ihs = "selbin.ihs.out";

	//
	int numThread = 1;
	unsigned int locus = 0;
	double cutoff = 0.05;
	unsigned int max_extend = 20000000;
	//flags
	bool calc_all = false;
	bool keep_low_freq = false;
	double min_maf = 0.05;
	bool alt = false;

	FILE* out_fp;
	uint64_t max_gap = 200000;
	uint64_t gap_scale = 20000;
	bool openmp_enabled = false;

	


	ofstream out_ihs;
	ofstream out_ehh;
	
	std::string logger_filename = "selbin.log";
    
	
    

protected:
	//void calc_EHH2(int locus, map<int, vector<int> > & m, map<int, priority_queue<pair<int, int>  > >* mphbs, bool downstream=false);
	void calc_EHH2(int locus, unordered_map<int, vector<int> > & m, bool downstream=false);
	void calc_EHH2_v2(int locus, unordered_map<int, vector<int> > & m, bool downstream=false);

	void calc_EHH(int locus, unordered_map<int, queue<pair<int, int> > >& outmap);
	void calc_EHH(int locus);

    void calc_iHS();
    void static thread_ihs(int tid, unordered_map<int, vector<int> >& m, unordered_map<int, vector<int> >& md, EHH* ehh_obj);
	void loadMPHB(string input_mphb_file, int numHaps, int numSnps, vector< map<int, priority_queue<pair<int, int>  > > > & mphbs); //vector< map<int, priority_queue<pair<int, int>  > > > & mphbs
	void loadMPHB(string input_mphb_file, int numHaps, int numSnps, map< int, map <int, priority_queue<pair<int, int>  > > > & mphbs); //vector< map<int, priority_queue<pair<int, int>  > > > & mphbs
	
private:
	HapMap hm;
	vector< map<int, priority_queue<pair<int, int>  > > > mphbs;
    //std::vector<std::vector<unsigned int> > all_positions;
    std::vector<unsigned int> loc_map;
	unsigned int numHaps;
	unsigned int numSnps;
    double* iHH0;
    double* iHH1;
    string* log_string_per_thread;
};

#endif