/*
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#include <cstdint>
#include <vector>
#include <string>
#include <bitset>
#include <fstream>
#include <unordered_map>
#include <map>
#include "selscan.hpp"

#ifndef HAPMAP_HPP
#define HAPMAP_HPP

#include <stdexcept>
#include <iostream>

#include "gzstream.h"
using namespace std;

struct map_entry {
//string chr;
  int locId;  
  int phyPos;
  int genPos;
  bool flipped; // required for algorithm where 0 and 1 are flipped
};

class HapMap
{
public:
//    HapMap();

    bool loadHap(const char* filename, double minmaf, std::vector<std::vector<unsigned int> >& all_positions, std::vector<unsigned int>& loc_map);
    bool loadHapMap(const char* filename, double minmaf);
    bool loadVCF(const char* filename, double minmaf);
    void readMapData(string filename, int expected_loci, bool USE_PMAP);

    void print();

    std::size_t numSnps() const { return m_numSnps; }
    std::size_t numHaps() const { return m_numHaps; }    
    
    std::vector<struct map_entry> data; // hold additional info
    std::vector<std::vector<unsigned int> > all_positions; // hold the 0, 1 matrix in position form
    std::vector<std::vector<unsigned int> > all_positions_two; // hold the 0, 1 matrix in position form

//    ~HapMap();
    double get1F(int loc);
    bool flip_to_minor = false;
    bool unphased = false;


private:
    int countFields( string& s);
    double getMAF(std::vector<unsigned int>& positions, std::vector<unsigned int>& positions_two);
    double getMAF(std::vector<unsigned int>& positions);
    inline bool skipLocus(std::vector<unsigned int>& positions);
    inline bool skipLocus(unsigned int locus);

    
protected:
    std::size_t m_numSnps;
    std::size_t m_numHaps;
    double minmaf;

    // std::size_t& ADVANCED_D =  m_numSnps;
    // std::size_t& ADVANCED_N = m_numHaps;

};

#endif 