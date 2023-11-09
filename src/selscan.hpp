#include "CLI11.hpp"
#include "assert.h"
#define DEBUG true
#include<vector>





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

	std::string input_filename_vcf = "/home/art/workspace/filezilla-sync/out.1000_100000_1.vcf";
	
	bool useVCF = false;
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



