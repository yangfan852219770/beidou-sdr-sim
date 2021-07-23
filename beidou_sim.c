#include<stdio.h>
#include<stdlib.h>
#include<stdbool.h>
#include<string.h>

#include "beidou_sim.h"
#include "lime_beidou.h"
#include "util.h"

int sinTable512[] = {
        2,   5,   8,  11,  14,  17,  20,  23,  26,  29,  32,  35,  38,  41,  44,  47,
        50,  53,  56,  59,  62,  65,  68,  71,  74,  77,  80,  83,  86,  89,  91,  94,
        97, 100, 103, 105, 108, 111, 114, 116, 119, 122, 125, 127, 130, 132, 135, 138,
        140, 143, 145, 148, 150, 153, 155, 157, 160, 162, 164, 167, 169, 171, 173, 176,
        178, 180, 182, 184, 186, 188, 190, 192, 194, 196, 198, 200, 202, 204, 205, 207,
        209, 210, 212, 214, 215, 217, 218, 220, 221, 223, 224, 225, 227, 228, 229, 230,
        232, 233, 234, 235, 236, 237, 238, 239, 240, 241, 241, 242, 243, 244, 244, 245,
        245, 246, 247, 247, 248, 248, 248, 249, 249, 249, 249, 250, 250, 250, 250, 250,
        250, 250, 250, 250, 250, 249, 249, 249, 249, 248, 248, 248, 247, 247, 246, 245,
        245, 244, 244, 243, 242, 241, 241, 240, 239, 238, 237, 236, 235, 234, 233, 232,
        230, 229, 228, 227, 225, 224, 223, 221, 220, 218, 217, 215, 214, 212, 210, 209,
        207, 205, 204, 202, 200, 198, 196, 194, 192, 190, 188, 186, 184, 182, 180, 178,
        176, 173, 171, 169, 167, 164, 162, 160, 157, 155, 153, 150, 148, 145, 143, 140,
        138, 135, 132, 130, 127, 125, 122, 119, 116, 114, 111, 108, 105, 103, 100,  97,
        94,  91,  89,  86,  83,  80,  77,  74,  71,  68,  65,  62,  59,  56,  53,  50,
        47,  44,  41,  38,  35,  32,  29,  26,  23,  20,  17,  14,  11,   8,   5,   2,
        -2,  -5,  -8, -11, -14, -17, -20, -23, -26, -29, -32, -35, -38, -41, -44, -47,
        -50, -53, -56, -59, -62, -65, -68, -71, -74, -77, -80, -83, -86, -89, -91, -94,
        -97,-100,-103,-105,-108,-111,-114,-116,-119,-122,-125,-127,-130,-132,-135,-138,
        -140,-143,-145,-148,-150,-153,-155,-157,-160,-162,-164,-167,-169,-171,-173,-176,
        -178,-180,-182,-184,-186,-188,-190,-192,-194,-196,-198,-200,-202,-204,-205,-207,
        -209,-210,-212,-214,-215,-217,-218,-220,-221,-223,-224,-225,-227,-228,-229,-230,
        -232,-233,-234,-235,-236,-237,-238,-239,-240,-241,-241,-242,-243,-244,-244,-245,
        -245,-246,-247,-247,-248,-248,-248,-249,-249,-249,-249,-250,-250,-250,-250,-250,
        -250,-250,-250,-250,-250,-249,-249,-249,-249,-248,-248,-248,-247,-247,-246,-245,
        -245,-244,-244,-243,-242,-241,-241,-240,-239,-238,-237,-236,-235,-234,-233,-232,
        -230,-229,-228,-227,-225,-224,-223,-221,-220,-218,-217,-215,-214,-212,-210,-209,
        -207,-205,-204,-202,-200,-198,-196,-194,-192,-190,-188,-186,-184,-182,-180,-178,
        -176,-173,-171,-169,-167,-164,-162,-160,-157,-155,-153,-150,-148,-145,-143,-140,
        -138,-135,-132,-130,-127,-125,-122,-119,-116,-114,-111,-108,-105,-103,-100, -97,
        -94, -91, -89, -86, -83, -80, -77, -74, -71, -68, -65, -62, -59, -56, -53, -50,
        -47, -44, -41, -38, -35, -32, -29, -26, -23, -20, -17, -14, -11,  -8,  -5,  -2
};

int cosTable512[] = {
        250, 250, 250, 250, 250, 249, 249, 249, 249, 248, 248, 248, 247, 247, 246, 245,
        245, 244, 244, 243, 242, 241, 241, 240, 239, 238, 237, 236, 235, 234, 233, 232,
        230, 229, 228, 227, 225, 224, 223, 221, 220, 218, 217, 215, 214, 212, 210, 209,
        207, 205, 204, 202, 200, 198, 196, 194, 192, 190, 188, 186, 184, 182, 180, 178,
        176, 173, 171, 169, 167, 164, 162, 160, 157, 155, 153, 150, 148, 145, 143, 140,
        138, 135, 132, 130, 127, 125, 122, 119, 116, 114, 111, 108, 105, 103, 100,  97,
        94,  91,  89,  86,  83,  80,  77,  74,  71,  68,  65,  62,  59,  56,  53,  50,
        47,  44,  41,  38,  35,  32,  29,  26,  23,  20,  17,  14,  11,   8,   5,   2,
        -2,  -5,  -8, -11, -14, -17, -20, -23, -26, -29, -32, -35, -38, -41, -44, -47,
        -50, -53, -56, -59, -62, -65, -68, -71, -74, -77, -80, -83, -86, -89, -91, -94,
        -97,-100,-103,-105,-108,-111,-114,-116,-119,-122,-125,-127,-130,-132,-135,-138,
        -140,-143,-145,-148,-150,-153,-155,-157,-160,-162,-164,-167,-169,-171,-173,-176,
        -178,-180,-182,-184,-186,-188,-190,-192,-194,-196,-198,-200,-202,-204,-205,-207,
        -209,-210,-212,-214,-215,-217,-218,-220,-221,-223,-224,-225,-227,-228,-229,-230,
        -232,-233,-234,-235,-236,-237,-238,-239,-240,-241,-241,-242,-243,-244,-244,-245,
        -245,-246,-247,-247,-248,-248,-248,-249,-249,-249,-249,-250,-250,-250,-250,-250,
        -250,-250,-250,-250,-250,-249,-249,-249,-249,-248,-248,-248,-247,-247,-246,-245,
        -245,-244,-244,-243,-242,-241,-241,-240,-239,-238,-237,-236,-235,-234,-233,-232,
        -230,-229,-228,-227,-225,-224,-223,-221,-220,-218,-217,-215,-214,-212,-210,-209,
        -207,-205,-204,-202,-200,-198,-196,-194,-192,-190,-188,-186,-184,-182,-180,-178,
        -176,-173,-171,-169,-167,-164,-162,-160,-157,-155,-153,-150,-148,-145,-143,-140,
        -138,-135,-132,-130,-127,-125,-122,-119,-116,-114,-111,-108,-105,-103,-100, -97,
        -94, -91, -89, -86, -83, -80, -77, -74, -71, -68, -65, -62, -59, -56, -53, -50,
        -47, -44, -41, -38, -35, -32, -29, -26, -23, -20, -17, -14, -11,  -8,  -5,  -2,
        2,   5,   8,  11,  14,  17,  20,  23,  26,  29,  32,  35,  38,  41,  44,  47,
        50,  53,  56,  59,  62,  65,  68,  71,  74,  77,  80,  83,  86,  89,  91,  94,
        97, 100, 103, 105, 108, 111, 114, 116, 119, 122, 125, 127, 130, 132, 135, 138,
        140, 143, 145, 148, 150, 153, 155, 157, 160, 162, 164, 167, 169, 171, 173, 176,
        178, 180, 182, 184, 186, 188, 190, 192, 194, 196, 198, 200, 202, 204, 205, 207,
        209, 210, 212, 214, 215, 217, 218, 220, 221, 223, 224, 225, 227, 228, 229, 230,
        232, 233, 234, 235, 236, 237, 238, 239, 240, 241, 241, 242, 243, 244, 244, 245,
        245, 246, 247, 247, 248, 248, 248, 249, 249, 249, 249, 250, 250, 250, 250, 250
};

double ant_pat_db[37] = {
        0.00,  0.00,  0.22,  0.44,  0.67,  1.11,  1.56,  2.00,  2.44,  2.89,  3.56,  4.22,
        4.89,  5.56,  6.22,  6.89,  7.56,  8.22,  8.89,  9.78, 10.67, 11.56, 12.44, 13.33,
        14.44, 15.56, 16.67, 17.78, 18.89, 20.00, 21.33, 22.67, 24.00, 25.56, 27.33, 29.33,
        31.56
};

// 生成测距码
void prn_code_gen(int *prn_code, int prn_number)
{
    int g1[PRN_SEQ_LEN], g2[PRN_SEQ_LEN];
    int x1[LFSR_LEN], x2[LFSR_LEN];
    int i,j,k;
    int c1, c2;
    int *prn_selector;
    int prn_selector_size;
    // 目前定义60行，最多4列，第一列代表长度
    int g2_phase[][4] = {
        {2,0,2},{2,0,3},{2,0,4},{2,0,5},{2,0,7},{2,0,8},{2,0,9},{2,0,10},
        {2,1,6},
        {2,2,3},{2,2,4},{2,2,5},{2,2,7},{2,2,8},{2,2,9},{2,2,10},
        {2,3,4},{2,3,5},{2,3,7},{2,3,8},{2,3,9},{2,3,10},
        {2,4,5},{2,4,7},{2,4,8},{2,4,9},{2,4,10},
        {2,5,7},{2,5,8},{2,5,9},{2,5,10},
        {2,7,8},{2,7,9},{2,7,10},
        {2,8,9},{2,8,10},
        {2,9,10},
        {3,0,1,6},
        {3,0,2,3},{3,0,2,5},{3,0,2,7},{3,0,2,9},{3,0,2,10},
        {3,0,3,4},{3,0,3,8},
        {3,0,4,5},{3,0,4,7},{3,0,4,9},{3,0,4,10},
        {3,0,5,8},
        {3,0,7,8},
        {3,0,8,9},{3,0,8,10},
        {3,1,2,6},{3,1,4,6},{3,1,6,8},
        {3,2,3,4},{3,2,3,8},{3,2,4,5},{3,2,4,7}
    };
    
    if(prn_number < 1 || prn_number > 60)
        return;
    
    prn_selector = g2_phase[prn_number - 1];
    prn_selector_size = prn_selector[0];
    // 初始化x1,x2相位
    // 初始相位为01010101010
    for(i = 0; i < LFSR_LEN; ++i){
        if(i%2 == 0)
            x1[i] = x2[i] = 0;
        else
            x1[i] = x2[i] = 1;
    }

    for(i = 0; i < PRN_SEQ_LEN; ++i){
        g1[i] = x1[LFSR_LEN - 1];
        // 初值为0，不会改变异或结果
        g2[i] = 0;
        for(k = 1; k <= prn_selector_size; ++k)
            g2[i] ^= x2[prn_selector[k]];
            
        /**
        for(int k = 0; k < LFSR_LEN; ++k)
            printf("%d ",x1[k]);
        printf("|%d",g1[i]);
        printf("|%d", ( g1[i] + g2[i] ) % 2);
        printf("|%d| ", g2[i]);
        for(int k = 0; k < LFSR_LEN; ++k)
            printf("%d ",x2[k]);
        putchar('\n');
        */

        c1 = x1[0] ^ x1[6] ^ x1[7] ^ x1[8] ^ x1[9] ^ x1[10];
        c2 = x2[0] ^ x2[1] ^ x2[2] ^ x2[3] ^ x2[4] ^ x2[7] ^ x2[8] ^ x2[10];
        for(j = LFSR_LEN -1; j >0; --j){
            x1[j] = x1[j-1];
            x2[j] = x2[j-1];
        }
        x1[0] = c1;
        x2[0] = c2; 

    }

    for(i = 0; i < PRN_SEQ_LEN; ++i)
        prn_code[i] = g1[i] ^ g2[i];
    return;
}

void nav_word_gen(unsigned long source_word, bool first_word_flag, int *nav_msg){
    // 将一个字拆为两个数组
    int *wrd1, wrd2[11];
    int nav_msg1[WORD_LEN/2], nav_msg2[WORD_LEN/2];
    // 11位掩码
    unsigned long mask_11 = 0x7FF;
    unsigned long wrd_hex;

    // 一个帧中的第一个字，第一个字的前15bits不做校验
    if(first_word_flag){
        // 移除后四位校验位
        wrd_hex = source_word >> 4;
        wrd1 = (int *) malloc(sizeof (int) * 15);
        // 转化为二进制形式
        hex2binary(wrd_hex >> 11, wrd1, 15);
        // 导航数据
        for(int i = 0; i < 15; ++i)
            nav_msg1[i] = wrd1[i];
    }
    else{
        // 移除后八位校验位
        wrd_hex = source_word >> 8;
        wrd1 = (int *) malloc(sizeof (int) * 11);
        // 转化为二进制形式
        hex2binary(wrd_hex >> 11, wrd1, 11);
        // 生成前11位的校验位
        BCH_code_gen(wrd1, 11, nav_msg1);
    }
    // 处理后11位
    hex2binary((wrd_hex & mask_11), wrd2, 11);
    BCH_code_gen(wrd2, 11, nav_msg2);

    // 合成30位导航电文
    // 第一个字不做交织编码
    if(first_word_flag){
        for(int i = 0; i < 15; ++i){
            nav_msg[i] = nav_msg1[i];
            nav_msg[i+15] = nav_msg2[i];
        }
    }
    else
        nav_msg_merge(nav_msg1, nav_msg2, nav_msg);

    free(wrd1);
}

/**
 * BCH(15,11,1)纠错码，共15位，11位数据信息，4位纠错码，纠错能力为1位
 * n1: Original navagation message, 11 bit or 15 bit
 * n1_length: 
 * nav_msg: After coded by BCH, 15 bit
 */
void BCH_code_gen(int *n1, int n1_length, int *nav_msg){
    // BCH_code 默认为0
    int BCH_code[BCH_CORRECT_CODE_LEN] = {0};
    int temp_d0, temp_d1;
    int i, j;

    for(i = 0; i < n1_length; ++i){
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

/**
 * 将两个15bit导航信息合成30bit(一个字)，串/并变换
 * n1: 导航电文信息1
 * n1_length: 15bit
 * n2：导航电文信息2
 * n2_length: 15bit
 * nav_msg: 将n1和n2合并后的导航信息(包含纠错码) 30bit
 * */
void nav_msg_merge(int *n1, int *n2, int *nav_msg){
    int i, j, k;
    j = 0, k = 0;
    // n1和n2 经串/并变换
    for(i = 0; i < WORD_LEN; ++i){
        // 偶数
        if( (i & 1) == 0)
            nav_msg[i] = n1[j++];
        else
            nav_msg[i] = n2[k++];
    }
    return;
}

void parse_BCH(int *nav_msg, int length){
    int BCH_code[4] = {0};
    int temp_D0, temp_D3;
    for(int i = 0; i < length; ++i){
        temp_D0 = BCH_code[0];
        temp_D3 = BCH_code[3];
        BCH_code[0] = nav_msg[i] ^ BCH_code[3];
        BCH_code[3] = BCH_code[2];
        BCH_code[2] = BCH_code[1];
        BCH_code[1] = temp_D3 ^ temp_D0;
    }
    printf("The parsing of BCH code: ");
    for(int i = 3; i >= 0; --i)
        printf("%d ", BCH_code[i]);
    putchar('\n');
}

/**
 * 将UTC时间转为北斗系统时间格式
 * U_time： UTC时间
 * bd_time： 格式化后的北斗时间
 * TODO 合法性校验
 */
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

void eph2sbf_D1(const ephemeris eph, const ionoutc_t ion, unsigned long sbf[][WORD_NUM]){

    //TODO 目前只模仿第一帧

    // 帧同步码
    unsigned long pre = 0x712UL;
    // 保留字段
    unsigned long rev = 0UL;
    /**
     * 设置模拟的固定时间，模拟的时间为:2021.1.1 0:0:0
     * second 19 bits
     * week 10 bits
     * mask 12 bits, 每位都为1，取second的后12位
     */
    unsigned long second = 0x69780UL;
    unsigned long week = 0x30EUL;
    //TODO 暂定为一天内的秒数
    unsigned long toc = 0UL;
    // TODO TGD1 暂存整数，82, 存补码形式
    unsigned long TGD1 = 0x82UL;
    // TODO TGD1 暂存整数，-15, 存补码形式
    unsigned long TGD2 = 0x3F1UL;

    ////////////////////////////////////////////////////////////
    // 第一帧
    ////////////////////////////////////////////////////////////

    /**
     * 第一个字
     * Pre | Rev | FraId | SOW(前8bit)
     * Pre(11bit)，帧同步码，默认为 11100010010
     * Rev(4bit)，保留字段，未参与计算
     * FraId(3bit)，子帧计数，第一帧为 1
     * SOW(前8bit)，周内秒计数
     */
    sbf[0][0] = pre << 19 | rev << 15 | 0x1UL << 12 | (second >> 12) << 4;
    /**
     * 第二个字
     * SOW | SatH1 | AODC | URAI
     * SOW(后12bit)，周内秒计数
     * SatH1(1bit)，卫星健康表示，0表示可用，1表示不可用
     * AODC(5bit)，时钟数据龄期，暂定为 27
     * URAI(4bit)，用户距离精度指数，暂定为 2
     */
    sbf[0][1] = ( second & 0xFFFUL) << 18 | 0 << 17 | 0x1BUL << 12 | 0x2UL << 8;
    /**
     * 第三个字
     * WN | toc
     * WN(13bit)，整周计数
     * toc(前9bit)，时钟时间
     */
    sbf[0][2] = week << 17 | toc << 8;
    /**
     * 第四个字
     * toc | TGD1 | TGD2
     * toc(后8bit)，时钟时间
     * TGD1(10bit) 星上设备时延差
     * TGD2(前4bit) 星上设备时延差
     */
    sbf[0][3] = (toc & 0xFF) << 22 | TGD1 << 12 | (TGD2 >> 6) << 8;
    /**
     * 第五个字
     * TGD2 | a0 | a1
     * TGD2(后6bit) 星上设备时延差
     */
    sbf[0][4] = (TGD2 & 0x3FUL) << 24;
    /**
     * 第六个字
     */
    sbf[0][5] = pre << 19 | rev << 15 | 0x1UL << 12 | (second >> 12) << 4;
     /**
      * 第七个字
      */
    sbf[0][6] = ( second & 0xFFFUL) << 18 | 0x1UL << 17 | 0x1BUL << 12 | 0x2UL << 8;
    /**
      * 第八个字
      */
    sbf[0][7] = week << 17 | toc << 8;
    /**
      * 第九个字
      */
    sbf[0][8] = (toc & 0xFF) << 22 | TGD1 << 12 | (TGD2 >> 6) << 8;
    /**
      * 第十个字
      */
    sbf[0][9] = (TGD2 & 0x3FUL) << 24;
}

void eph2sbf_D2(const ephemeris eph, const ionoutc_t ion, unsigned long sbf[][WORD_NUM]){

    //TODO 目前只模仿第一帧

    // 帧同步码
    unsigned long pre = 0x712;

    /**
     * 设置模拟的固定时间，模拟的时间为:2021.1.1 0:0:0
     * week 10 bits
     */
    unsigned long second = 0x9780;
    unsigned long week = 0x20E;
    //TODO 暂定为一天内的秒数
    unsigned long toc = 0x15;
    // TODO TGD1 暂存整数，82, 存补码形式
    unsigned long TGD1 = 0x52UL;
    // TODO TGD1 暂存整数，-15, 存补码形式
    unsigned long TGD2 = 0x3F1UL;

    ////////////////////////////////////////////////////////////
    // 第一帧
    ////////////////////////////////////////////////////////////

    /**
     * 第一个字
     * Pre | Rev | FraId | SOW(前8bit)
     * Pre(11bit)，帧同步码，默认为 11100010010
     * Rev(4bit)，保留字段，未参与计算
     * FraId(3bit)，子帧计数，第一帧为 1
     * SOW(前8bit)，周内秒计数
     */
    sbf[0][0] = pre << 19 | 0 << 15 | 0x1UL << 12 | (second >> 12) << 4;
     /** 第二个字
      * SOW | Pnum1 | SatH1 | AODC
      * SOW(后12bit)，周内秒计数
      * Pnum1(4bit), 页面号
      * SatH1(1bit)，卫星健康表示，0表示可用，1表示不可用
      * AODC(5bit)，时钟数据龄期，暂定为 24
      */
    sbf[0][1] = ( second & 0xFFFUL) << 18 | 0x1UL << 14 | 0 << 13 | 0x18UL << 8;
     /** 第三个字
      * URAI | WN | toc
      * URAI(4bit)，用户距离精度指数，暂定为 2
      * WN(13bit)，整周计数
      * toc(前5bit)，时钟时间
      */
    sbf[0][2] = 0x2UL << 26 | week << 13 | (toc >> 12) << 8 ;
    /**
    * 第四个字
    * toc | TGD1
    * toc(后12bit)，时钟时间
    * TGD1(10bit) 星上设备时延差
    */
    sbf[0][3] = (toc & 0xFFF) << 18 | TGD1 << 8;
    /**
     * 第五个字
     * TGD2 | Rev
     * TGD2(后6bit) 星上设备时延差
     * Rev 保留字
     */
    sbf[0][4] = TGD2 << 20 | 0< 8;
    /**
     * 第六个字
     */
    sbf[0][5] = pre << 19 | 0 << 15 | 0x1UL << 12 | (second >> 12) << 4;
    /**
     * 第七个字
     */
    sbf[0][6] = ( second & 0xFFFUL) << 18 | 0x2UL << 14;
    /**
      * 第八个字
      */
     sbf[0][7] = 0;
    /**
      * 第九个字
      */
    sbf[0][8] = 0;
    /**
      * 第十个字
      */
    sbf[0][9] = 0;
}

int allocate_channel(beidou_channel *chan, const ephemeris eph, ionoutc_t ionoutc, beidou_time bdt){
    // 目前只模仿prn号为1-5的卫星
    int start = 1, end = 5;
    int j = 0;
    for(int i = start; i <= end; ++i){
        chan[j].prn_num = i;
        // 生成测距码
        prn_code_gen(chan[j].prn_code, i);
        // 从星历生成子帧数据
        eph2sbf_D2(eph, ionoutc, chan[j].subframe);
        // 生成导航数据
        nav_msg_gen(bdt, &chan[j], 1);
        chan[j].carr_phase = (unsigned int)25050620;
        ++j;
    }
    return j;
}

int nav_msg_gen(beidou_time bd_time, beidou_channel *chan, int init){

    // 第一帧第一个字单独处理
    nav_word_gen(chan->subframe[0][0], true, chan->subframe_word_bits[0]);
    nav_word_gen(chan->subframe[0][5], true, chan->subframe_word_bits[5]);

    // 处理剩余的字
    int start = 1;
    for(int i = start; i < WORD_NUM; ++i)
        if(5 != i)
            nav_word_gen(chan->subframe[0][i], false, chan->subframe_word_bits[i]);
}

void init(beidou_channel *chan){

    // TODO 暂时为静态
    chan->f_carr = 108.158;
    chan->f_code = 2046000.14;
    chan->code_phase = 10.34;

    chan->iword = 7;
    chan->ibit = 26;

    chan->icode = 0;

    chan->prn_code_bit = chan->prn_code[(int)chan->code_phase] * 2 - 1;
    chan->data_bit = chan->subframe_word_bits[chan->iword][chan->ibit] * 2  - 1;

    // 20bits NH码
    chan->NH_code = 0x4D4EUL;
    // 从高位开始为第一位
    chan->iNH_code = chan->icode;
    chan->NH_code_bit = (int)((chan->NH_code >> (19 - chan->iNH_code)) & 0x1UL) * 2 - 1;

    return;
}

void *beidou_task(void *arg){
    sim_t *s = (sim_t *)arg;

    //TODO eph目前不用
    ephemeris eph[1][1];
    beidou_time bdt;
    // 模拟4个信道
    beidou_channel chan[MAX_CHAN_SIM];

    int ip, qp;
    int iTable;
    short *iq_buff = NULL;
    int iq_buff_size;
    int isamp;

    // TODO ionoutc目前不用
    ionoutc_t ionoutc;
    // 模拟信号的时间，模拟300s
    int duration = 30000;
    int idura;

    int i, j;

    // 采样周期
    double delt;

    int gain[MAX_CHAN_SIM];
    double timer = 0;
    iTable = 0;

    iq_buff_size = NUM_IQ_SAMPLES;
    iq_buff = calloc(2*iq_buff_size, 2);
    if (iq_buff==NULL)
    {
        printf("ERROR: 分配 I/Q buff失败.\n");
        goto exit;
    }
    delt = 1.0/(double)TX_SAMPLERATE;

    ////////////////////////////////////////////////////////////
    // 初始化发送信道
    ////////////////////////////////////////////////////////////

    // 清除信道信息
    for(i = 0; i < MAX_CHAN_SIM; ++i)
        chan[i].prn_num = 0;

    // 分配信道
    allocate_channel(chan, eph[0][0], ionoutc, bdt);

    ////////////////////////////////////////////////////////////
    // 生成基带信号
    ////////////////////////////////////////////////////////////

    //初始化信道数据
    for( i = 0; i < MAX_CHAN_SIM; ++i){
        if(chan[i].prn_num > 0){
            init(&chan[i]);
            chan[i].carr_phasestep = (int)(512 * 65536.0 * chan[i].f_carr * delt);
            gain[i] = 64 + (i + 1) * 8;
        }
    }

    for(idura = 1; idura < duration; ++idura){
        //TODO 天线增益、路径损失

        // 生成发射数据
        for(isamp = 0; isamp < iq_buff_size; ++isamp){
            int i_acc = 0;
            int q_acc = 0;

            for(i = 0; i < MAX_CHAN_SIM; ++i){
                if(chan[i].prn_num > 0){
                    iTable = (chan[i].carr_phase >> 16) & 511;

                    ip =  chan[i].data_bit * chan[i].prn_code_bit * cosTable512[iTable];
                    qp =  chan[i].data_bit * chan[i].prn_code_bit * sinTable512[iTable];

                    i_acc += ip;
                    q_acc += qp;

                    ////////////////////////////////////////////////////////////
                    // PRN码，NH码，导航信息 处理
                    ////////////////////////////////////////////////////////////
                    // 更新PRN码相位
                    chan[i].code_phase += chan[i].f_code * delt;

                    if(chan[i].code_phase >= PRN_SEQ_LEN){
                        chan[i].code_phase -= PRN_SEQ_LEN;
                        ++chan[i].icode;
                        // 2 prn_code = 1 nav bit
                        if(chan[i].icode >= 2){
                            chan[i].icode = 0;
                            ++chan[i].ibit;
                            // 30 bits = 1 word
                            if(chan[i].ibit >= WORD_LEN){
                                chan[i].ibit = 0;
                                //
                                printf("%d\n",chan[i].iword);

                                ++chan[i].iword;

                                // TODO 10word 等于1帧
                                if(chan[i].iword >= WORD_NUM){
                                    chan[i].iword = 0;
                                }
                            }
                            chan[i].data_bit = chan[i].subframe_word_bits[chan[i].iword][chan[i].ibit] * 2 -1;
                        }

                    }
                    // PRN数据处理
                    chan[i].prn_code_bit = chan[i].prn_code[(int)chan[i].code_phase] * 2 -1;
                    chan[i].carr_phase += chan[i].carr_phasestep;
                }
            }
            // Scaled by 2^7
           // i_acc = (i_acc+64)>>7;
            //q_acc = (q_acc+64)>>7;

            // 存储I/Q buff
            iq_buff[isamp*2] = i_acc;
            iq_buff[isamp*2+1] = q_acc;
        }
        ////////////////////////////////////////////////////////////
        // 写入发射缓存
        ///////////////////////////////////////////////////////////

        if(!s->beidou.ready){
            printf("Bei Dou signal is ready.\n");
            s->beidou.ready = 1;
            pthread_cond_signal(&(s->beidou.initialization_done));
        }


        // 等待FIFO缓存
        pthread_mutex_lock(&(s->beidou.lock));
        while (!is_fifo_write_ready(s))
            pthread_cond_wait(&(s->fifo_write_ready), &(s->beidou.lock));
        pthread_mutex_unlock(&(s->beidou.lock));
        // 向FIFO缓存写入
        memcpy(&(s->fifo[s->head * 2]), iq_buff, NUM_IQ_SAMPLES * 2 * sizeof(short));

        s->head += (long)NUM_IQ_SAMPLES;
        if (s->head >= FIFO_LENGTH)
            s->head -= FIFO_LENGTH;
        pthread_cond_signal(&(s->fifo_read_ready));
        // 分配发射信道
        allocate_channel(chan, eph[0][0], ionoutc, bdt);
        timer +=0.1;
    }
    s->finished = true;
    free(iq_buff);

    exit:
        return (NULL);
}