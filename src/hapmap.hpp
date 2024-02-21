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


#define BITSET_SIZE 1048576
#include <stdexcept>
#include <iostream>

#include "gzstream.h"
using namespace std;

struct map_entry {
  int locId;  
  int phyPos;
  int genPos;
  bool flipped;
};

class HapMap
{
public:
//    HapMap();
    // static std::size_t querySnpLength(const char* filename);
    // double geneticPosition(std::size_t line) const;
    // long long unsigned int physicalPosition(std::size_t line) const;
    // void loadMap(const char* mapFileName);
    // std::string lineToId(std::size_t line) const;
    // std::size_t idToLine(const std::string& id) const;
    bool loadHap(const char* filename, double minmaf, std::vector<std::vector<unsigned int> >& all_positions, std::vector<unsigned int>& loc_map);
    bool loadHapMap(const char* filename, const char* mapfile, double minmaf);
    bool loadVCF(const char* filename, double minmaf);
    void print();

    std::size_t numSnps() const { return m_numSnps; }
    std::size_t numHaps() const { return m_numHaps; }    
    
    std::vector<struct map_entry> mentries; // hold additional info
    std::vector<std::vector<unsigned int> > all_positions; // hold the 0, 1 matrix in position form
    std::vector<std::bitset<BITSET_SIZE> > all_bitsets; // hold the xor
    std::vector<std::vector<unsigned int> > all_xors; // hold the xor

    
//    ~HapMap();
    double getMAF(int loc);

private:
    int countFields( string& s);
protected:
    // std::map<std::size_t, std::string> m_idMap;
    // unsigned long long* m_physPos;
    // double* m_genPos;
    // std::map<std::size_t, std::string> m_idMap;
    std::vector<unsigned int> monomorphic;
    std::size_t m_numSnps;
    std::size_t m_numHaps;


    // std::size_t& ADVANCED_D =  m_numSnps;
    // std::size_t& ADVANCED_N = m_numHaps;

};

#endif 