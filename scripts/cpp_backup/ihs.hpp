#include <string>
#include <iostream>
#include <vector>

#include "CLI11.hpp"
#include "selscan.hpp"


#ifndef IHS_H
#define IHS_H

class NSL: public Selscan {
// private:
// 	std::string input_filename;
// 	std::string output_filename;

public:
	IHS();
	void cli_prepare(CLI::App * subapp);
	void ihs(std::string input, std::string output);
	void exec();
};

#endif