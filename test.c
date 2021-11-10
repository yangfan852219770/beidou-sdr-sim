// B1I 信号模拟
#include<stdio.h>
#include<stdbool.h>

#include "beidou_sim.h"

int main(int argc, char const *argv[])
{
    ephemeris eph[EPHEM_ARRAY_SIZE][MAX_SAT];

    ionoutc_t ionoutc;
    int neph,ieph;
    char navfile[MAX_CHAR] = "2021.1.3.txt";
    neph = readRinexNavAll(1, eph, &ionoutc, navfile);

    return 0;
}
