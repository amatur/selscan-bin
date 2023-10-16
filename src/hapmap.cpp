
#include "hapmap.hpp"
#include <cassert>
#include <cstdlib>
#include <cstdint>
#include <cstring>
#include <iostream>

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


bool HapMap::loadHapMap(const char* filename, const char* mapfile, double minmaf, std::vector<std::vector<unsigned int> >& all_positions, std::vector<struct map_entry>& mentries)
{
    //std::ifstream fmap(mapfile, std::ios::in );
    std::ifstream f(filename, std::ios::in );
    if (!f.good())
    {
        std::cerr << "ERROR: Cannot open file or file not found: " << filename << std::endl;
        return false;
    }
 
    std::string line;
    m_numSnps = 0;
    m_numHaps = 0;
    unsigned int actual_snp_id = 0;
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
        double maf = positions.size()*1.0/m_numHaps ;
        std::cout<<"Loc: "<<actual_snp_id<<"1 freq: "<<maf<<std::endl;
        if(maf < minmaf || 1-maf < minmaf){
            //skip
            std::cout<<"WARNING: skipping site" << actual_snp_id<< std::endl;

        }else{
            ++m_numSnps;
            all_positions.push_back(positions); //check if all 0
            struct map_entry mentry;
            mentry.genPos = actual_snp_id;
            mentry.phyPos = actual_snp_id;
            mentry.locId = actual_snp_id;
            mentries.push_back(mentry);

        }
        ++actual_snp_id;
    }
    
    // mentries_arr = new struct map_entry[m_numSnps];
    // for(int k = 0; k<m_numSnps; k++){
    //     mentries_arr[k] = mentries[k];
    // }
    //std::memcpy(mentries_arr, mentries.data(), sizeof(struct map_entry)*m_numSnps);
    f.close();
    return true;
}

bool HapMap::loadHap(const char* filename, double minmaf, std::vector<std::vector<unsigned int> >& all_positions)
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
            //std::cout<<"WARNING: monomorphic site"<< m_numSnps << std::endl;
            //this->monomorphic.push_back(m_numSnps);
        }
        double maf = positions.size()*1.0/m_numHaps ;
        if(maf < minmaf || 1-maf < minmaf){
            //skip
            //std::cout<<"WARNING: skipping site" << std::endl;

        }else{
            ++m_numSnps;
            all_positions.push_back(positions); //check if all 0
        }

        
        
    }
    f.close();
    return true;
}
