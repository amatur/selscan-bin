
#include<string>
using namespace std;
#ifndef UTILS_H
#define UTILS_H

// std::string GetDirectory (const std::string& path)
// {
//     size_t found = path.find_last_of("/\\");
//     string ret = (path.substr(0, found));
//     if(ret==""){
//         return "./";
//     }
//     return ret;
// }

#define PBSTR "||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||"
#define PBWIDTH 60

void static printProgress(double percentage) {
    int val = (int) (percentage * 100);
    int lpad = (int) (percentage * PBWIDTH);
    int rpad = PBWIDTH - lpad;
    printf("\r%3d%% [%.*s%*s]", val, lpad, PBSTR, rpad, "");
    fflush(stdout);
}


#endif