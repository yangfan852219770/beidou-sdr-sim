// B1I 信号模拟
#include<stdio.h>
#include<stdbool.h>

#include "beidou_sim.h"

// Neumann-Hoffman码
int NH_code[] = {
    0, 0, 0, 0, 0, 1, 0, 0, 1, 1,
    0, 1, 0, 1, 0, 0, 1, 1, 1, 0
};

int main(int argc, char const *argv[])
{

    int a1[11] = {0,1,1,1,1,0,1,0,0,1,0},a2[15];
    BCH_code_gen(a1, 11, a2);
    printf("Msg:");
    for(int i = 0; i < 15; ++i)
        printf("%d ", a2[i]);
    putchar('\n');
    parse_BCH(a2,15);
    return 0;

    // 时间转换测试
    /*
    UTC_time u_time;
    beidou_time bd_time;
    u_time.year = 2021, u_time.month = 4, u_time.day = 30;
    u_time.hour = 13, u_time.minute = 10, u_time.second = 10;
    UTC2beidou_time(&u_time, &bd_time);
    printf("周数为:%d,归零次数为:%d, 秒数为:%.2f\n", bd_time.week, bd_time.back0time, bd_time.second);
    */

     return 0;
}
