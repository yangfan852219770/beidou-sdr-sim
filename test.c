// B1I 信号模拟
#include<stdio.h>
#include<stdbool.h>

#include "beidou_sim.h"

#define BCH_CORRECT_CODE_LEN (4)
// Neumann-Hoffman码
int NH_code[] = {
    0, 0, 0, 0, 0, 1, 0, 0, 1, 1,
    0, 1, 0, 1, 0, 0, 1, 1, 1, 0
};
void BCH_code_gen(int *n1, int n1_length, int *nav_msg){
    // BCH_code 默认为0
    int BCH_code[BCH_CORRECT_CODE_LEN] = {0};
    int temp_d0, temp_d1;
    int i, j;

    for(i = 1; i < n1_length; ++i){
        temp_d0 = BCH_code[0];
        temp_d1 = BCH_code[1];
        BCH_code[0] = n1[i] ^ BCH_code[3];
        BCH_code[1] = temp_d0 ^ BCH_code[0];
        BCH_code[3] = BCH_code[2];
        BCH_code[2] = temp_d1;
    }

    /*
     * 第3位为纠错码的第一位，BCH数组倒着输出
    printf("4bit BCH code:");
    for(i = BCH_CORRECT_CODE_LEN - 1; i >= 0; --i)
        printf("%d ", BCH_code[i]);
    putchar('\n');
    */

    
    for(i = 0; i < n1_length; ++i)
        nav_msg[i] = n1[i];
    // 将BCH_code与n1合成导航信息,第3位为纠错码的第一位
    for( j = BCH_CORRECT_CODE_LEN - 1; j >=0; --j)
        nav_msg[i++] = BCH_code[j];
    return;
}

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
