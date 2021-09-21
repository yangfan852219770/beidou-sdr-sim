// B1I 信号模拟
#include<stdio.h>
#include<stdbool.h>

#include "beidou_sim.h"
void UTC2beidou_time(const UTC_time *U_time, beidou_time *bd_time){
    int day[12] = {0,31,59,90,120,151,181,212,243,273,304,334};
    int start_year, start_month, start_day;
    int d_days, leap_days, d_years;
    
    start_year = 2006, start_month = 1, start_day = 1;
    leap_days = 0;
    // 相差的年份
    d_years = U_time->year - start_year;

    // 计算闰月的累计天数
    while(start_year < U_time->year){
        if((start_year % 4 == 0 && start_year % 100 != 0) || (start_year % 400 == 0))
            ++leap_days;
        ++start_year;
    }
    
    // 相差的天数
    d_days =  d_years * 365 + day[U_time->month - 1] + U_time->day + leap_days - 1;
    // 计算北斗周数
    bd_time->week = (d_days / 7) % 8192;
    // 计算北斗周数归零次数,8192周后归零
    bd_time->back0time = bd_time->week / 8192;
    // 计算北斗秒数
    bd_time->second = (double)(d_days % 7)* SECONDS_IN_DAY + U_time->hour * SECONDS_IN_HOUR
                        + U_time->minute * SECONDS_IN_MINUTE + U_time->second;
    return;
}
// Neumann-Hoffman码
int NH_code[] = {
    0, 0, 0, 0, 0, 1, 0, 0, 1, 1,
    0, 1, 0, 1, 0, 0, 1, 1, 1, 0
};



int main(int argc, char const *argv[])
{
    ephemeris eph[EPHEM_ARRAY_SIZE][MAX_SAT];

    ionoutc_t ionoutc;
    int neph,ieph;
    char navfile[MAX_CHAR] = "2021.1.3.txt";
    neph = readRinexNavAll(eph, &ionoutc, navfile);

    return 0;
}
