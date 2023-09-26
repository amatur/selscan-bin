// main.cpp

                            
#include <iostream>
#include <algorithm>
#include "bm.h"

using namespace std;

void Print(unsigned n)
{
    cout << n << endl;;
}

int main(void)
{
    bm::bvector<>   bv;    

    bv[10] = true;
    bv[100] = true;
    bv[10000] = true;

    bm::bvector<>::enumerator en = bv.first();
    bm::bvector<>::enumerator en_end = bv.end();

    while (en < en_end)
    {
        cout << *en << endl;
        ++en;  // pre-increment - fastest way to increment enumerator
    }

    en = bv.first();

    // This is not the fastest way to do the job, because for_each 
    // often will try to calculate difference between iterators,
    // which is expensive for enumerators.
    // But it can be useful for some STL loyal applications. 

    std::for_each(en, en_end, Print);

    return 0;
}

              