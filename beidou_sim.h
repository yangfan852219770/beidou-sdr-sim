#ifndef BEIDOU_SIM_H
#define BEIDOU_SIM_H
#include<stdbool.h>

// RINEX 文件中卫星的最大数量
#define MAX_SAT (60)
#define PI 3.1415926535898

// PRN码长度
#define PRN_SEQ_LEN (2046)
// 线性反馈移位寄存器长度
#define LFSR_LEN (11)
//BCH纠错码长度，4位
#define BCH_CORRECT_CODE_LEN (4)
// 导航电文一个字的长度，30bit
#define WORD_LEN (30)
// 导航电文一个子的一半长度，15bit
// #define HALF_WORD_LEN (15)
// 子帧个数
#define SUBFRAME_NUM (5)
// 字个数，暂时模仿10个
#define WORD_NUM (10)

#define SECONDS_IN_WEEK (604800.0)
#define SECONDS_IN_HALF_WEEK (302400.0)
#define SECONDS_IN_DAY (86400.0)
#define SECONDS_IN_HOUR (3600.0)
#define SECONDS_IN_MINUTE (60.0)

// 载波频率
#define CARR_FREQ (1561.098e6)
// PRN码速率
#define PRN_CODE_FREQ (2.046e6)
// 载波频率转PRN速率
#define CARR_TO_PRN (1.0/763.0)
// 模拟的通道数量
#define MAX_CHAN_SIM (5)
// NH码长度
#define NH_CODE_LEN (20)

#define R2D 57.2957795131

// UTC 时间
typedef struct{
    int year;
    int month;
    int day;
    int hour;
    int minute;
    double second;
} UTC_time;

// 北斗时间格式
typedef struct{
    // 北斗周数归零次数
    int back0time;
    // 北斗周数，从2006年1月1日起，8191后重新归0
    int week;
    // 周内秒计数
    double second;

} beidou_time;

// 星历
typedef struct
{
    int flag;
    beidou_time toc; // 时钟时间
    beidou_time toe; // 星历中的时间
} ephemeris;

typedef struct
{
    int flag;
} ionoutc_t;

typedef struct 
{
    int prn_num; // 卫星编号
    int prn_code[PRN_SEQ_LEN]; // prn码

    double f_carr;	/*< Carrier frequency */
    double f_code;	/*< Code frequency */
    unsigned int carr_phase; /*< Carrier phase */
    int carr_phasestep;	/*< Carrier phasestep */
    double code_phase; /*< Code phase */
    int icode;	/*!< initial code */

    int prn_code_bit;
    beidou_time bd_t0; // 北斗时间
    /**
     * TODO 目前只模仿第一个子帧
     */
    unsigned long subframe[1][WORD_NUM];
    // TODO 暂时不用
    unsigned int subfranme_num; // 子帧号
    // TODO 暂时不用
    unsigned long subframe_word[3];
    // 目前只模仿第一帧
    int subframe_word_bits[WORD_NUM][WORD_LEN];
    // 导航信息数据
    int data_bit;
    // 第几位导航数据
    int ibit;
    // 第几个字
    int iword;
    // NH码
    unsigned long NH_code;
    // NH码位数据
    int NH_code_bit;
    // 第几位NH码
    int iNH_code;
} beidou_channel;

/** Generate PRN code
 * prn_code: A array
 * prn_number: PRN number(Or the number of satellite)
 */
void prn_code_gen(int *prn_code, int prn_number);

void nav_word_gen(unsigned long source_word, bool first_word_flag, int *nav_msg);

/**
 * BCH(15,11,1)纠错码，共15位，11位数据信息，4位纠错码，纠错能力为1位
 * n1: Original navagation message, 11 bit or 15 bit
 * n1_length: 
 * nav_msg: After coded by BCH, 15 bit
 */
void BCH_code_gen(int *n1, int n1_length, int *nav_msg);
/**
 * 将两个15bit导航信息合成30bit(一个字)，串/并变换
 * n1: 导航电文信息1
 * n2：导航电文信息2
 * nav_msg: 将n1和n2合并后的导航信息(包含纠错码) 30bit
 */
void nav_msg_merge(int *n1, int *n2, int *nav_msg);
// test
void parse_BCH(int *nav_msg, int length);

/**
 * 将UTC时间转为北斗系统时间格式
 * U_time： UTC时间
 * bd_time： 格式化后的北斗时间
 */
void UTC2beidou_time(const UTC_time *U_time, beidou_time *bd_time);

/**
 * 星历转导航电文子帧
 * !!! 目前只模仿第一帧
 * @param eph
 * @param ion
 * @param sbf
 */
void eph2sbf(const ephemeris eph, const ionoutc_t ion, unsigned long sbf[][WORD_NUM]);

/**
 * 分配信道
 * @param chan
 * @param eph
 * @return
 */
int allocate_channel(beidou_channel *chan, const ephemeris eph, ionoutc_t ionoutc, beidou_time bdt);

/**
 * 生成导航信息
 * @param bd_time
 * @param chan
 * @param init
 * @return
 */
int nav_msg_gen(beidou_time bd_time, beidou_channel *chan, int init);

void init(beidou_channel *chan);

/**
 * 信号模拟
 * @param arg
 * @return
 */
void *beidou_task(void *arg);
#endif 