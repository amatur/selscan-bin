#include <string>
#include <iostream>
#include <vector>
#include <map>
#include <thread>

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
    std::string output_filename = "selbin.ihh.out";// <outfile>.nsl[.alt].out
	//
	int numThread = 1;
	unsigned int locus = 0;
	double cutoff = 0.9;
	unsigned int max_extend = 20000000;
	//flags
	bool calc_all = false;
	bool keep_low_freq = false;
	double min_maf = 0.05;
	bool alt = false;

	FILE* out_fp;
	uint64_t max_gap = 200000;
	uint64_t gap_scale = 20000;


protected:
	void calc_EHH2(int locus, map<int, vector<int> > & m, bool downstream=false);
	void calc_EHH(int locus);
    void calc_iHS();
    void static thread_ihs(int tid, map<int, vector<int> >& m, map<int, vector<int> >& md, EHH* ehh_obj);

private:
	HapMap hm;
    //std::vector<std::vector<unsigned int> > all_positions;
    std::vector<unsigned int> loc_map;
	unsigned int numHaps;
	unsigned int numSnps;
    double* iHH0;
    double* iHH1;
    string* logg;
};

#endif