/*
 * Hapbin: A fast binary implementation EHH, iHS, and XPEHH
 * Copyright (C) 2014  Colin MacLean <s0838159@sms.ed.ac.uk>
 *
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

#include "hapmap.hpp"
#include <cassert>
#include <cstdlib>

// HapMap::HapMap()
//     : m_numSnps{}
//     , m_numHaps{}
// {

// }

// HapMap::~HapMap()
// {
// }

// void HapMap::loadMap(const char* mapFilename)
// {
//     assert(m_numSnps != 0ULL); //Must loadHap first.
//     if (m_physPos)
//         delete m_physPos;
//     if (m_genPos)
//         delete m_genPos;
//     m_physPos = new unsigned long long[m_numSnps];
//     m_genPos = new double[m_numSnps];
        
//     std::ifstream file(mapFilename);

//     std::string line;
    
//     std::size_t lineNum = 0ULL;
//     while (std::getline(file, line)) 
//     {
//         if (lineNum == m_numSnps)
//         {
//             std::cerr << "WARNING: Map file has more loci than hap file!" << std::endl;
//             return;
//         }
//         std::vector<std::string> split = splitString(line, ' ');
// 	if (split.size() == 1)
//             split = splitString(line, '\t');
//         if (split.size() == 4) {
//             m_idMap[lineNum] = split[1];
//             m_genPos[lineNum] = atof(split[2].c_str());
//             m_physPos[lineNum] = strtoull(split[3].c_str(), 0, 10);
//             ++lineNum;
//         }
//     }
//     if (lineNum != m_numSnps)
//     {
//         std::cerr << "ERROR: Map file must have the same number of loci as hap file! Hap file has " << m_numSnps << " SNPs. Map file has " << lineNum << " SNPs." << std::endl;
//         if (lineNum == 0)
//             std::cerr << "Perhaps the Map file format is wrong?" << std::endl;
//         abort();
//     }
// }


bool HapMap::loadHap(const char* filename, std::vector<std::vector<unsigned int> >& all_positions)
{
    std::ifstream f(filename, std::ios::in );
    if (!f.good())
    {
        std::cerr << "ERROR: Cannot open file or file not found: " << filename << std::endl;
        return false;
    }
 
    std::string line;
    m_numSnps = 0;
    m_numHaps = 0;
    while (std::getline(f, line)) {
        if(DEBUG) std::cout<<line<<std::endl;
        std::vector<unsigned int> positions;
 
        std::istringstream rowStream(line);
        char value;
        int pos=0;
        while (rowStream >> value) {        
            if(value=='1'){
                positions.push_back(pos++);
            }else if(value=='0'){
                pos++;
            }   
        }


        if(m_numHaps==0){
            m_numHaps = pos;
        }else{
            if(pos!=m_numHaps){ //integrity check: all lines must be of same length
                std::cerr<<"ERROR: site"<<m_numSnps<<" has incorrect number of haplotypes."<<std::endl;
                exit(2);
            }
        }

        if(positions.size()==0 or positions.size()==m_numHaps){
            //This implies site is monomorphic
            std::cout<<"WARNING: monomorphic site"<< m_numSnps << std::endl;
            //this->monomorphic.push_back(m_numSnps);
        }
        ++m_numSnps;
        all_positions.push_back(positions); //check if all 0
    }
    f.close();
    return true;
}
