
#include "hapmap.hpp"
#include <cassert>
#include <cstdlib>
#include <cstdint>
#include <cstring>
#include <iostream>
using namespace std;


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

double HapMap::get1F(int loc){
    if(data[loc].flipped){
        return (m_numHaps-all_positions[loc].size())*1.0/m_numHaps;

    }else{
        return all_positions[loc].size()*1.0/m_numHaps;
    }    
}

double HapMap::getMAF(std::vector<unsigned int>& positions,std::vector<unsigned int>& positions_two ){
    unsigned int freq1 = positions.size();
    unsigned int freq2 = positions_two.size();
    unsigned int freq0 = this->m_numHaps-freq1-freq2;
    unsigned int minim = std::min({freq1, freq2, freq0});
    return minim*1.0/this->m_numHaps;
}

double HapMap::getMAF(std::vector<unsigned int>& positions){
    unsigned int freq1 = positions.size();
    unsigned int freq0 = this->m_numHaps-freq1;
    unsigned int minim = std::min({freq1, freq0});
    return minim*1.0/this->m_numHaps;
}

/***
 * warning: assumes m_numHaps and positions are already correctly computed
*/
inline bool HapMap::skipLocus(std::vector<unsigned int>& positions){
    double derived_af = positions.size()/this->m_numHaps;
    return (derived_af < this->minmaf || 1-derived_af < this->minmaf || positions.size()==1 || positions.size()==m_numHaps-1 || positions.size()==1 || positions.size() == m_numHaps-1);
}

/***
 * warning: assumes m_numHaps and positions are already correctly computed
*/
inline bool HapMap::skipLocus(unsigned int locus){
    std::vector<unsigned int> positions = this->all_positions[locus];
    double derived_af = positions.size()/this->m_numHaps;
    return (derived_af < this->minmaf || 1-derived_af < this->minmaf || positions.size()==1 || positions.size()==m_numHaps-1 || positions.size()==1 || positions.size() == m_numHaps-1);
}

/**
 * Load vcf -> physPos, get positions matrix // here physPos field is populated, but genPos field is populated from map file
*/
bool HapMap::loadVCF(const char* filename, double minmaf)
{
    this->minmaf = minmaf;
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

    for (int locus = 0; locus < nloci; locus++)
    {
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


        std::vector<unsigned int> positions; // holds positions of 1's in this locus
        std::vector<unsigned int> positions_two; // holds positions of 2's in this locus
        

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
                char allele = '0';
                if (allele1 == '1' && allele2 == '1'){
                    //allele = '2';
                    positions_two.push_back(field);
                }
                else if (allele1 == '1' || allele2 == '1'){
                    //allele = '1';
                    positions.push_back(field);
                }
                //else{
                    //allele is 0 is implied 
                //}


                //data->data[field][locus] = allele;
            }
            else{
                // data->data[2 * field][locus] = allele1;
                // data->data[2 * field + 1][locus] = allele2;
                //TODO

                if(allele1=='1'){
                    positions.push_back(2 * field);
                    //cout<<2 * field<<" ";
                }
                if(allele2=='1'){
                    positions.push_back(2 * field + 1);
                    //cout<<2 * field + 1<<" ";
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
        vector<unsigned int> copy_pos;
        
        //cout<<maf<<" maf"<<endl;
        //std::cout<<"Loc: "<<locus<<"1 freq: "<<maf<<std::endl;

        if(derived_af < minmaf || 1-derived_af < minmaf || positions.size()==1 || positions.size()==m_numHaps-1 || positions.size()==1 || positions.size() == m_numHaps-1){
            //if(positions.size()==0 or positions.size()==m_numHaps){
            //skip
            //std::cout<<"WARNING: skipping site" << locus<< std::endl;

        }else{
            struct map_entry mentry;
            if(unphased){
                all_positions_two.push_back(positions_two); //check if all 0
                all_positions.push_back(positions); //check if all 0
                mentry.flipped = false;
            }else{
                if( flip_to_minor && derived_af < 0.5){
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
                    all_positions.push_back(positions); //check if all 0
                    mentry.flipped = false;
                }
            }
            
            // mentry.genPos = physpos;
            mentry.phyPos = physpos;
            mentry.locId = locus;
            data.push_back(mentry);
        }
    }

    fin.close();
    this-> m_numSnps = all_positions.size();
    cout<<"retaining "<<this-> m_numSnps <<endl;
    // this-> m_numHaps = nhaps;
    this->unphased = unphased;

    return true;
}


/***
 * impute hap format with map information imputed
*/
bool HapMap::loadHapMap(const char* hapfile, double minmaf)
{
    //std::ifstream fmap(mapfile, std::ios::in );
    std::ifstream f(hapfile, std::ios::in );
    if (!f.good())
    {
        std::cerr << "ERROR: Cannot open file or file not found: " << hapfile << std::endl;
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

        // if(positions.size()==0 or positions.size()==m_numHaps){
        //     //This implies site is monomorphic
        //     std::cout<<"WARNING: monomorphic site"<< m_numSnps << std::endl;
        // }

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
            data.push_back(mentry);
        }
        ++actual_snp_id;
    }
    
    f.close();
    std::cout<<"Finished loading file."<<std::endl;
    return true;
}


/**
 * IMPUTE hap format (matches hapbin format)
*/
bool HapMap::loadHap(const char* filename, double minmaf, std::vector<std::vector<unsigned int> >& all_positions, std::vector<unsigned int>& loc_map)
{
    this->minmaf = minmaf;
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

/***
 * untested
*/
//reads in map data and also does basic checks on integrity of format
//returns a populated MapData structure if successful
//throws an exception otherwise
void  HapMap::readMapData(string filename, int expected_loci, bool USE_PMAP)
{
    igzstream fin;
    cerr << "Opening " << filename << "...\n";
    fin.open(filename.c_str());

    if (fin.fail())
    {
        cerr << "ERROR: Failed to open " << filename << " for reading.\n";
        throw 0;
    }

    string line;
    int nloci = 0;
    int num_cols = 4;
    int current_cols = 0;
    while (getline(fin, line))
    {
        nloci++;
        current_cols = countFields(line);
        if (current_cols != num_cols)
        {
            cerr << "ERROR: line " << nloci << " of " << filename << " has " << current_cols
                 << ", but expected " << num_cols << ".\n";
            throw 0;
        }
    }

    if (nloci != expected_loci)
    {
        cerr << "ERROR: Expected " << expected_loci << " loci in map file but found " << nloci << ".\n";
        throw 0;
    }

    fin.clear(); // clear error flags
    fin.close();
    fin.open(filename.c_str());

    if (fin.fail())
    {
        cerr << "ERROR: Failed to open " << filename << " for reading.\n";
        throw 0;
    }

    cerr << "Loading map data for " << nloci << " loci\n";

    double Mb = 1000000.0;
    
    string skip_chr;
    int skip_int;
    
    int ith_snp = 0; 
    for (int locus = 0; locus < this->m_numSnps; locus++) //numSnps == nloci
    {
        if (skipLocus(locus)){
            fin >> skip_chr;
            fin >> skip_int;
            fin >> skip_int;
            fin >> skip_int;
        }else{
            fin >> skip_chr; //fin >> this->data[locus]->chr;
            fin >> this->data[ith_snp].locId; //data->locusName[locus];
            fin >> this->data[ith_snp].genPos;
            fin >> this->data[ith_snp].phyPos;
            if (USE_PMAP) this->data[ith_snp].genPos = double(this->data[ith_snp].phyPos)/Mb;
            ith_snp++;
        }
    }

    fin.close();
}