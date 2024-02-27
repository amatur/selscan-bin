
#include "hapmap.hpp"
#include <cassert>
#include <cstdlib>
#include <cstdint>
#include <cstring>
#include <iostream>
#include <bits/stdc++.h>

//#include "dynamic_bitset.hpp"
using namespace std;
// HapMap::HapMap()
//     : m_numSnps{}
//     , m_numHaps{}.vscode
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





void HapMap::print(){
    for(int i = 0 ; i< this->m_numSnps; i++){
        cout<<"Loc "<<i<<":";
        for (unsigned int v : this->all_positions[i]){
            cout<<v<<" ";
        }
        cout<<endl;
    }
    
}

int HapMap::countFields(string& str)
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

double HapMap::getMAF(int loc){
    return all_positions[loc].size()*1.0/m_numHaps;
}

bool HapMap::loadVCF(const char* filename, double minmaf)
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
    uint64_t nloci = 0;
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
    m_numHaps = nhaps;
    //m_numSnps = nloci; 


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


    //std::bitset<BITSET_SIZE> prev_bitset;

    
    for (int locus = 2011111; locus < nloci; locus++)
    {
        printProgress((locus+1)*1.0/nloci);
        uint64_t physpos;
        for (int i = 0; i < numMapCols; i++) {
            fin >> junk;
            if (i == 0 && junk[0] == '#') {
                skipLine = true;
                break;
            }else{
                if(i==1){ //column for physpos //0 for chr //2 for id
                    std::istringstream reader(junk);
                    reader >> physpos;
                    //cout<<"PHYSPOS:"<<junk<<" "<<physpos<<endl;
                }
            }
            
        }
        if (skipLine) {
            getline(fin, junk);
            skipLine = false;
            locus--;
            continue;
        }

        

        std::vector<unsigned int> positions; // holds positions for this locus

        
        // std::bitset<BITSET_SIZE> curr_bitset;
        // sul::dynamic_bitset<> prev_bitset_d(m_numHaps+1);
        // sul::dynamic_bitset<> curr_bitset_d(m_numHaps+1);
        // sul::dynamic_bitset<> zero_bitset(m_numHaps+1);



        bool flip_enabled = true;

        
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
                    //curr_bitset.set(2 * field);
                    //curr_bitset_d[2 * field] = 1;
                }
                if(allele2=='1'){
                    positions.push_back(2 * field + 1);
                    //curr_bitset.set(2 * field + 1);
                    //curr_bitset_d[2 * field+1] = 1;
                }
            }
        }
        //cout<<endl;
        if(positions.size()==0 or positions.size()==m_numHaps){
            //This implies site is monomorphic
            //std::cout<<"WARNING: monomorphic site"<< locus << std::endl;
            //this->monomorphic.push_back(m_numSnps);
        }
        if(positions.size()==1 or positions.size()==m_numHaps-1){
            //std::cout<<"WARNING: MAC must be > 1 : site " << locus << std::endl;
        }
        double derived_af = positions.size()*1.0/m_numHaps ;
        //cout<<maf<<" maf"<<endl;
        //std::cout<<"Loc: "<<locus<<"1 freq: "<<maf<<std::endl;
        if(derived_af < minmaf || 1-derived_af < minmaf || positions.size()==1 || positions.size()==m_numHaps-1 || positions.size()==1 || positions.size()==m_numHaps-1){
        //if(positions.size()==0 or positions.size()==m_numHaps){
            //skip
            //std::cout<<"WARNING: skipping site" << locus<< std::endl;

        }else{

            /*
            all_bitsets.push_back(curr_bitset);
            if(all_positions.size() == 0){
                //all_bitsets.push_back(curr_bitset); // for 0th loc, just push the bitset. for others do xor
                prev_bitset_d=curr_bitset_d;

                sul::dynamic_bitset b(curr_bitset_d);
                vector<unsigned int> v;

                int first = b.find_first();
                while(first<m_numHaps){
                    v.push_back(first);
                    b[first] = 0;
                    first = b.find_first();
                }
                all_xors.push_back(v);
            }else{
                
                
                cout<<"before xor:"<<prev_bitset_d.count()<<" "<<curr_bitset_d.count()<<endl;
                sul::dynamic_bitset b(prev_bitset_d ^ curr_bitset_d);
                cout<<"after xor:"<<b.count()<<" "<<curr_bitset_d.count()<<endl;

                //auto b = (  curr_bitset); // uncomment to test the speedup of xor

                vector<unsigned int> v;
                int first = prev_bitset_d.find_first();
                while(first<m_numHaps){
                    v.push_back(first);
                    prev_bitset_d[first] = 0;
                    first = prev_bitset_d.find_first();
                }
                all_xors.push_back(v);

                //cout<< m_numSnps <<":" <<all_bitsets.back().count()<<" "<<curr_bitset.count()<<endl; 
                prev_bitset_d = curr_bitset_d;
                cout<<"after xor xor:"<<prev_bitset_d.count()<<" "<<curr_bitset_d.count()<<endl;

                curr_bitset_d.clear();
                cout<<"after xor xor xor:"<<prev_bitset_d.count()<<" "<<curr_bitset_d.count()<<endl;

            }
            */


            struct map_entry mentry;
            if(flip_enabled &&  derived_af > 0.5){
                // vector<unsigned int> v;
                // //cout<<curr_bitset.count()<<" ";
                // curr_bitset.flip();
                // int first = curr_bitset._Find_first();
                // while(first<m_numHaps){
                //     v.push_back(first);
                //     curr_bitset[first] = 0;
                //     first = curr_bitset._Find_first();
                //     //cout<<first<<" ";
                // }
                // //cout<<6404-v.size()<<" ";
                // //cout<<endl;
                // all_positions.push_back(v);

                if(all_positions.size()==0){
                    all_xors.push_back(positions);
                }else{
                    vector<unsigned int> curr_xor;
                    std::set_symmetric_difference(positions.begin(), positions.end(), all_positions.back().begin(), all_positions.back().end(),
                                  std::back_inserter(curr_xor));
                    // std::set_symmetric_difference(positions.begin(), positions.end(), all_positions.back().begin(), all_positions.back().end(),
                    //               std::inserter(curr_xor, curr_xor.end()));
                                  
                    all_xors.push_back(curr_xor);
                }
                

                vector<unsigned int> copy_pos;
                int cnt = 0;
                for(int i = 0; i<m_numHaps; i++){
                    unsigned int curr = positions[cnt];
                    if(i==curr){
                        cnt++;
                    }else{
                        copy_pos.push_back(i);
                    }
                }
                all_positions.push_back(copy_pos); 
                
                
                mentry.flipped = true;
            }else{
                if(all_positions.size()==0){
                    all_xors.push_back(positions);
                }else{
                    vector<unsigned int> curr_xor;
                    std::set_symmetric_difference(positions.begin(), positions.end(), all_positions.back().begin(), all_positions.back().end(), std::back_inserter(curr_xor));
                    // std::set_symmetric_difference(positions.begin(), positions.end(), all_positions.back().begin(), all_positions.back().end(),
                    //               std::inserter(curr_xor, curr_xor.end()));
                    all_xors.push_back(curr_xor);
                }
                
                all_positions.push_back(positions); //check if all 0
                mentry.flipped = false;
            }

            
            
            // mentry.genPos = physpos;
            mentry.phyPos = physpos;
            mentry.locId = locus;
            mentries.push_back(mentry);
        }
    }

    fin.close();
    this-> m_numSnps = all_positions.size();
    cout<<"retaining "<<this-> m_numSnps <<endl;
    // this-> m_numHaps = nhaps;

    //vector<vector<int> >all_xors;
    // for(auto b : all_bitsets){
    //     vector<unsigned int> v;

    //     int first = b._Find_first();
    //     while(first!=BITSET_SIZE){
    //         v.push_back(first);
    //         b[first] = 0;
    //         first = b._Find_first();
    //     }
    //     all_xors.push_back(v);
    // }


//================================================================================================
    //comment this block for debugging
    for(int i = 0 ; i< this->m_numSnps; i++){
        //std::unordered_set<unsigned int> position_set(this->all_positions[i].begin(), this->all_positions[i].end());
        //this->all_positions.push_back(position_set);

        // cout<<"Loc (xor)"<<i<<":";
        // for (unsigned int v : all_xors[i]){
        //     cout<<v<<" ";
        // }
        // cout<<endl;
        // cout<<"Loc "<<i<<":";
        // for (unsigned int v : this->all_positions[i]){
        //     cout<<v<<" ";
        // }
        // cout<<endl;
    }
    

    return true;
}

// bool HapMap::loadHapMap(const char* hapfile, const char* mapfile, double minmaf)
// {
//     //std::ifstream fmap(mapfile, std::ios::in );
//     std::ifstream f(hapfile, std::ios::in );
//     if (!f.good())
//     {
//         std::cerr << "ERROR: Cannot open file or file not found: " << hapfile << std::endl;
//         return false;
//     }
 
//     std::string line;
//     m_numSnps = 0;
//     m_numHaps = 0;
//     unsigned int actual_snp_id = 0;
//     while (std::getline(f, line)) {
//         //if(DEBUG) std::cout<<line<<std::endl;
//         std::vector<unsigned int> positions;
 
//         std::istringstream rowStream(line);
//         char value;
//         int pos=0;
//         while (rowStream >> value) {        
//             if(value=='1'){
//                 positions.push_back(pos++);
//             }else if(value=='0'){
//                 pos++;
//             }   
//         }


//         if(m_numHaps==0){
//             m_numHaps = pos;
//         }else{
//             if(pos!=m_numHaps){ //integrity check: all lines must be of same length
//                 std::cerr<<"ERROR: site"<<m_numSnps<<" has incorrect number of haplotypes."<<std::endl;
//                 exit(2);
//             }
//         }

//         if(positions.size()==0 or positions.size()==m_numHaps){
//             //This implies site is monomorphic
//             std::cout<<"WARNING: monomorphic site"<< m_numSnps << std::endl;
//             //this->monomorphic.push_back(m_numSnps);
//         }
//         double maf = positions.size()*1.0/m_numHaps ;
//         std::cout<<"Loc: "<<actual_snp_id<<"1 freq: "<<maf<<std::endl;
//         if(maf < minmaf || 1-maf < minmaf){
//             //skip
//             std::cout<<"WARNING: skipping site" << actual_snp_id<< std::endl;

//         }else{
//             ++m_numSnps;
//             all_positions.push_back(positions); //check if all 0
//             struct map_entry mentry;
//             mentry.genPos = actual_snp_id;
//             mentry.phyPos = actual_snp_id;
//             mentry.locId = actual_snp_id;
//             mentries.push_back(mentry);

//         }
//         ++actual_snp_id;
//     }
    
//     // mentries_arr = new struct map_entry[m_numSnps];
//     // for(int k = 0; k<m_numSnps; k++){
//     //     mentries_arr[k] = mentries[k];
//     // }
//     //std::memcpy(mentries_arr, mentries.data(), sizeof(struct map_entry)*m_numSnps);
//     f.close();
//     std::cout<<"Finished loading file."<<std::endl;
//     return true;
// }



// bool HapMap::loadHap(const char* filename, double minmaf, std::vector<std::vector<unsigned int> >& all_positions, std::vector<unsigned int>& loc_map)
// {
//     std::ifstream f(filename, std::ios::in );
//     if (!f.good())
//     {
//         std::cerr << "ERROR: Cannot open file or file not found: " << filename << std::endl;
//         return false;
//     }
 
//     std::string line;
//     m_numSnps = 0;
//     m_numHaps = 0;
//     unsigned int act_snp_count = 0;
//     std::bitset<BITSET_SIZE> prev_bitset;
//     while (std::getline(f, line)) {
//         //if(DEBUG) std::cout<<line<<std::endl;
//         std::vector<unsigned int> positions;
//         std::bitset<BITSET_SIZE> current_bitset;
        

 
//         std::istringstream rowStream(line);
//         char value;
//         int pos=0;
//         while (rowStream >> value) {        
//             if(value=='1'){
//                 current_bitset.set(pos);
//                 positions.push_back(pos++);
//             }else if(value=='0'){
//                 pos++;
//             }   
//         }

//         if(m_numHaps==0){
//             m_numHaps = pos;
//         }else{
//             if(pos!=m_numHaps){ //integrity check: all lines must be of same length
//                 std::cerr<<"ERROR: site"<<m_numSnps<<" has incorrect number of haplotypes."<<std::endl;
//                 exit(2);
//             }
//         }

//         if(positions.size()==0 or positions.size()==m_numHaps){
//             //This implies site is monomorphic
//             //std::cout<<"WARNING: monomorphic site"<< m_numSnps << std::endl;
//             //this->monomorphic.push_back(m_numSnps);
//         }
//         double maf = positions.size()*1.0/m_numHaps ;
//         if(maf < minmaf || 1-maf < minmaf){
//             //skip
//             //std::cout<<"WARNING: skipping site" << std::endl;

//         }else{
//             if(m_numSnps == 0){
//                 all_bitsets.push_back(current_bitset);
//             }else{
//                 all_bitsets.push_back(prev_bitset ^ current_bitset);
//                 cout<< m_numSnps <<":" <<all_bitsets.back().count()<<" "<<current_bitset.count()<<endl; 
//                 prev_bitset = current_bitset;

//             }
//             loc_map.push_back(act_snp_count);
//             ++m_numSnps;
//             all_positions.push_back(positions); //check if all 0

            
            
            
//         }
//         act_snp_count++;
//     }
//     f.close();
//     std::cout<<"Finished loading file."<<std::endl;

//     return true;
// }
