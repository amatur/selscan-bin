//need to add check if n^2 causes overflow
// main.cpp
// version 1:
// sep 26
// upstream only
// dont store it at all
// just store EHH values 


#include <thread>                            
#include <iostream>
#include <algorithm>

#include <fstream>

#include <cmath>

#include <memory>
#include <iostream>
#include <sstream>
#include <vector>
#include <set>
#include <map>


// #include "benchmark.hpp"

#include "selscan.hpp"
#include "nsl.hpp"
#include "ehh.hpp"


using namespace std;

// #define ADVANCED_N 4
// #define ADVANCED_D 6404


// #define FILENAME "data/out5000.imp#define FILENAME "data/out5000.impute.hap"
// #define ADVANCED_N 6404
// #define ADVANCED_D 5000ute.hap"
// #define ADVANCED_N 6404
// #define ADVANCED_D 5000


// #define FILENAME "/storage/home/aur1111/s/dataset/out20000.impute.hap"
// #define ADVANCED_N 6404
// #define ADVANCED_D 20000

// #define FILENAME "/storage/home/aur1111/s/msdir/ms5k6k.hap2"
// #define ADVANCED_N 6000
// #define ADVANCED_D 5000


// #define FILENAME "data/out500.impute.hap"
// #define ADVANCED_N 6404
// #define ADVANCED_D 500

// #define ADVANCED_N 5
// #define ADVANCED_D 5
// #define FILENAME "test.hap"

#define ADVANCED_N 4
#define ADVANCED_D 6
#define FILENAME "test4x6.hap"


// #define ADVANCED_N 3
// #define ADVANCED_D 6404

// #define ADVANCED_N 6404
// #define ADVANCED_D 20000

// #define ADVANCED_N 6404
// #define ADVANCED_D 5000

double start_time=0;

Selscan * parse_args(int argc, char** argv, vector<Selscan *> tools) {
	// Main command
	CLI::App app{"Selbin is a software for scan. For more details on kff format, please refer to https://github.com/"};
	app.require_subcommand(1);
    app.get_formatter()->column_width(40);
	app.get_formatter()->label("REQUIRED", "(REQUIRED)");
	app.get_formatter()->label("TEXT:FILE", "<text>");
	app.get_formatter()->label("TEXT", "<text>");
	app.get_formatter()->label("FLOAT", "<float>");
	app.get_formatter()->label("INT", "<int>");
	app.get_formatter()->label("BOOL", "<bool>");
	app.get_formatter()->label("UINT", "<int>");

	CLI::Option * help =	app.get_help_ptr();

	// Subcommands prepare
	for (Selscan * tool : tools) {
		tool->cli_prepare(&app);
	}

	// Parsing and return status if wrong
	try {
    app.parse(argc, argv);
	} catch (const CLI::ParseError &e) {
    auto val = app.exit(e);
		if (val != 0) {
			exit(val);
		}
	}

	// Help detection
	if (!help->empty()) {
		return nullptr;
	}

	// Read the command line return
	for (Selscan * tool : tools) {
		if (tool->subapp->parsed()) {
			// Help on tool triggered
			if (not tool->subapp->get_help_ptr()->empty()) {
				return nullptr;
			} else {
				return tool;
			}
		}
	}

	return nullptr;
}              

int main(int argc, char** argv) {
	
	// --- Prepare tools ---
	vector<Selscan *> tools;
	tools.push_back(new NSL());
	// tools.push_back(new IHS());
	tools.push_back(new EHH());

	// Get the one selected
	Selscan * tool = parse_args(argc, argv, tools);
   
	if (tool != nullptr)
		tool->exec();

	for (Selscan * tool : tools)
		delete tool;

	return 0;
}