
#include "hapmap.hpp"
#include <cassert>
#include <cstdlib>
#include <cstdint>
#include <cstring>
#include <iostream>
using namespace std;
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




int HapMap::countFields(const string &str)
{
    string::const_iterator it;
    int result;
    int numFields = 0;
    int seenChar = 0;
    for (it = str.begin() ; it < str.end(); it++)
    {
        result = isspace(*it);
        if (result == 0 && seenChar == 0)
        {
            numFields++;
            seenChar = 1;
        }
        else if (result != 0)
        {
            seenChar = 0;
        }
    }
    return numFields;
}

bool HapMap::loadHapMapVCF(const char* filename, const char* mapfile, double minmaf, std::vector<std::vector<unsigned int> >& all_positions, std::vector<struct map_entry>& mentries)
{
    bool unphased = false;
    igzstream fin;
    std::cerr << "Opening " << filename << "...\n";
    fin.open(filename);

    if (fin.fail())
    {
        std::cerr << "ERROR: Failed to open " << filename << " for reading.\n";
        throw 0;
    }

    int numMapCols = 9;
    //int fileStart = fin.tellg();
    std::string line;
    int nloci = 0;
    int previous_nhaps = -1;
    int current_nhaps = 0;
    //Counts number of haps (cols) and number of loci (rows)
    //if any lines differ, send an error message and throw an exception
    while (getline(fin, line))
    {
        if (line[0] == '#') {
            continue;
        }
        //getline(fin,line);
        //if(fin.eof()) break;
        nloci++;
        current_nhaps = countFields(line);
        //cout << "nloci: " << current_nhaps << endl;
        if (previous_nhaps < 0)
        {
            previous_nhaps = current_nhaps;
            continue;
        }
        else if (previous_nhaps != current_nhaps)
        {
            cerr << "ERROR: line " << nloci << " of " << filename << " has " << current_nhaps
                 << " fields, but the previous line has " << previous_nhaps << " fields.\n";
            throw 0;
        }
        previous_nhaps = current_nhaps;
    }

    fin.clear(); // clear error flags
    //fin.seekg(fileStart);
    fin.close();

    fin.open(filename);

    if (fin.fail())
    {
        cerr << "ERROR: Failed to open " << filename << " for reading.\n";
        throw 0;
    }

    int nhaps = unphased ? (current_nhaps - numMapCols) : (current_nhaps - numMapCols) * 2;
    int nfields = (current_nhaps - numMapCols);
    cerr << "Loading " << nhaps << " haplotypes and " << nloci << " loci...\n";

    //  std::vector<unsigned int> positions;
 
    //     std::istringstream rowStream(line);
    //     char value;
    //     int pos=0;
    //     while (rowStream >> value) {        
    //         if(value=='1'){
    //             positions.push_back(pos++);
    //         }else if(value=='0'){
    //             pos++;
    //         }   
    //     }


    string junk;
    char allele1, allele2, separator;
    bool skipLine = false;

    for (int locus = 0; locus < nloci; locus++)
    {
        for (int i = 0; i < numMapCols; i++) {
            fin >> junk;
            if (i == 0 && junk[0] == '#') {
                skipLine = true;
                break;
            }
        }
        if (skipLine) {
            getline(fin, junk);
            skipLine = false;
            locus--;
            continue;
        }

        std::vector<unsigned int> positions; // holds positions for this locus
        
        for (int field = 0; field < nfields; field++)
        {
            fin >> junk;
            allele1 = junk[0];
            separator = junk[1];
            allele2 = junk[2];
            if ( (allele1 != '0' && allele1 != '1') || (allele2 != '0' && allele2 != '1') )
            {
                cerr << "ERROR: Alleles must be coded 0/1 only.\n";
                cerr << allele1 << " " << allele2 << endl;
                throw 0;
            }

            //if(separator != '|'){
            //    cerr << "ERROR:  Alleles must be coded 0/1 only.\n";
            //    throw 0;
            //}
            if(unphased){
            }
            else{
                // data->data[2 * field][locus] = allele1;
                // data->data[2 * field + 1][locus] = allele2;
                //TODO

                if(allele1=='1'){
                    positions.push_back(2 * field);
                }
                if(allele2=='1'){
                    positions.push_back(2 * field + 1);
                }
            }
        }

        if(positions.size()==0 or positions.size()==m_numHaps){
            //This implies site is monomorphic
            std::cout<<"WARNING: monomorphic site"<< locus << std::endl;
            //this->monomorphic.push_back(m_numSnps);
        }
        double maf = positions.size()*1.0/m_numHaps ;
        std::cout<<"Loc: "<<locus<<"1 freq: "<<maf<<std::endl;
        if(maf < minmaf || 1-maf < minmaf){
            //skip
            std::cout<<"WARNING: skipping site" << locus<< std::endl;

        }else{
            all_positions.push_back(positions); //check if all 0
            struct map_entry mentry;
            mentry.genPos = locus;
            mentry.phyPos = locus;
            mentry.locId = locus;
            mentries.push_back(mentry);
        }

        
    }

    fin.close();
    this-> m_numSnps = nloci;
    this-> m_numHaps = nhaps;

    return true;
}


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
        //if(DEBUG) std::cout<<line<<std::endl;
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
    std::cout<<"Finished loading file."<<std::endl;
    return true;
}

bool HapMap::loadHap(const char* filename, double minmaf, std::vector<std::vector<unsigned int> >& all_positions, std::vector<unsigned int>& loc_map)
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
    unsigned int act_snp_count = 0;
    while (std::getline(f, line)) {
        //if(DEBUG) std::cout<<line<<std::endl;
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
            loc_map.push_back(act_snp_count);
            ++m_numSnps;
            all_positions.push_back(positions); //check if all 0
        }
        act_snp_count++;
    }
    f.close();
    std::cout<<"Finished loading file."<<std::endl;

    return true;
}
