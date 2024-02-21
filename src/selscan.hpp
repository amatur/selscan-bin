#include "CLI11.hpp"
#include "assert.h"
#define DEBUG true
#include<vector>
#include <fstream>

using namespace std; 

// REMINDERS
// 1. unsigned int vs. signed int


/*! \file selscan.h
    \brief Functions and utilities from selscan
*/



//public stuffs


#ifndef SEL__H__INCLUDED__
#define SEL__H__INCLUDED__


class Selscan {
public:
	CLI::App * subapp = nullptr;
	virtual void cli_prepare(CLI::App * subapp) = 0;
	virtual void exec() = 0;
	virtual ~Selscan() {};

	bool useVCF = true;
	std::string input_filename_vcf = "";
	//std::string input_filename_vcf = "/home/art/workspace/filezilla-sync/out.1000_100000_1.vcf";
	
};








// void splitList(std::vector<unsigned>& a, std::vector<unsigned>& b){

// }



// void calculateMAF(){
//    // do a remap: 
//    //bool skip 
//    // if MAF < minMaf Skip
//    // make all skip
//    // then if pos > 
//    // if < thr add to skip
//    // all that u need to skip u skip
   
// }

#endif



