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

    /*
    int a1[11] = {0,0,1,0,0,1,1,1,1,1,1},a2[15];
    BCH_code_gen(a1, 11, a2);
    printf("Msg:");
    for(int i = 0; i < 15; ++i)
        printf("%d ", a2[i]);
    putchar('\n');
    // parse_BCH(a2,15);
    return 0;
    // 时间转换测试
     */
    UTC_time u_time;
    beidou_time bd_time;
    u_time.year = 2021, u_time.month = 1, u_time.day = 1;
    u_time.hour = 0, u_time.minute = 0, u_time.second = 0;
    UTC2beidou_time(&u_time, &bd_time);
    printf("周数为:%d,归零次数为:%d, 秒数为:%.2f\n", bd_time.week, bd_time.back0time, bd_time.second);

     return 0;
}
