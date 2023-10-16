#include <string>
#include <iostream>
#include <vector>
#include <map>
#include<thread>



using namespace std;

#include "CLI11.hpp"
#include "selscan.hpp"


#ifndef EHH_H
#define EHH_H

class EHH: public Selscan {
// private:
// 	std::string input_filename;
// 	std::string output_filename;

public:
	EHH();
	void cli_prepare(CLI::App * subapp);
	void init();
	void exec();
	~EHH();

	std::string input_filename_hap = "data/out500.impute.hap";

    //std::string input_filename_hap = "test4x6.hap";
    std::string input_filename_map = "test4x6.map";
    std::string output_filename = "outfile";// <outfile>.nsl[.alt].out
	//
	int numThread = 1;
	unsigned int locus = 0;
	double cutoff = 0.0;
	unsigned int max_extend = 200;
	//flags
	bool calc_all = false;
	bool keep_low_freq = false;
	double min_maf = 0.05;
	bool alt = false;

	FILE* out_fp;

protected:
	void calc_EHH2(int locus, map<int, vector<int> > & m, bool downstream=false);
	void calc_EHH(int locus);
	
	// void calc_EHH_downstream(int locus, map<int, vector<int> > & m);
    float calc_iHS();
    void static thread_ihs(int tid, map<int, vector<int> >& m, map<int, vector<int> >& md, EHH* ehh_obj);


private:
    std::vector<std::vector<unsigned int> > all_positions;
    std::vector<unsigned int> loc_map;
	unsigned int ADVANCED_N;
	unsigned int ADVANCED_D;

    // float* ehh0;
    // float* ehh1;
    // float* ehh0_downstream;
    // float* ehh1_downstream;

    float* iHH0;
    float* iHH1;
    // float* iHH0_downstream;
    // float* iHH1_downstream;

    string* logg;
};

#endif