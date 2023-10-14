#include <string>
#include <iostream>
#include <vector>

#include "CLI11.hpp"
#include "selscan.hpp"
// #include "utils.hpp"

using namespace std;


#ifndef NSL_H
#define NSL_H

class NSL: public Selscan {
// private:
// 	std::string input_filename;
// 	std::string output_filename;

public:
	NSL();
	void cli_prepare(CLI::App * subapp);
	void init();
	void exec();

	void calc_nSL(int locus);
	void calc_nSL_all();
	~NSL();

    std::string input_filename_hap = "test4x6.hap";
    std::string input_filename_map = "test4x6.map";
    std::string output_filename = "outfile";// <outfile>.nsl[.alt].out

private:
	long* nsl0;
	long* nsl1;
	//
	std::vector<std::vector<unsigned int> > all_positions;
	unsigned int ADVANCED_N;
	unsigned int ADVANCED_D;
	//
	int numThread = 1;
	unsigned int locus = 0;
	double cutoff = 0.0;
	unsigned int max_extend_nsl = 200;
	//flags
	bool calc_all = false;
	bool keep_low_freq = false;
	double min_maf = 0.05;
	bool alt = false;

protected:
	void calc_nSL_upstream(int locus, map<int, vector<int> > & m);
	void calc_nSL_downstream(int locus, map<int, vector<int> > & m);

template <typename T> void print_a(T* arr, string name="v"){
    std::cout<<"vector: "<<name<<": ";
    for (int i=0; i<ADVANCED_N ; i++){
        std::cout<<arr[i]<<" ";
    }
    std::cout<<std::endl;
}

void print_v(vector<unsigned int> v, string name="v"){
    std::cout<<"vector: "<<name<<": ";
    for (int i: v){
        std::cout<<i<<" ";
    }
    std::cout<<std::endl;
}

};

#endif