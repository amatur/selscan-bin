#include <string>
#include <iostream>
#include <vector>

#include "CLI11.hpp"
#include "selscan.hpp"
#include "hapmap.hpp"
// #include "utils.hpp"

using namespace std;


#ifndef FEATURE_H
#define FEATURE_H


struct hapblock {
//   int left;
  int right;
  int numSeq;
//   int seqId;
};

class Feature: public Selscan {

public:
	Feature();
	void cli_prepare(CLI::App * subapp);
	void init();
	void exec();

	~Feature();
	string input_file_mphb;

};

#endif