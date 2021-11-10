#ifndef BEIDOU_SIM_H
#define BEIDOU_SIM_H
#include<stdbool.h>

/*! \brief Maximum length of a line in a text file (RINEX, motion) */
#define MAX_CHAR (100)
// TODO RINEX 文件中GEO卫星的最大数量
#define MAX_SAT (60)
#define MAX_SAT_EPH (60)
#define EPHEM_ARRAY_SIZE (13) // for daily GPS broadcast ephemers file (brdc)

#define GM_EARTH 3.986004418e14
#define OMEGA_EARTH 7.2921150e-5
#define PI 3.1415926535898

#define WGS84_RADIUS	6378137.0
#define WGS84_ECCENTRICITY 0.0818191908426

//星历参数缩放因子
#define POW2_M5  0.03125
#define POW2_M6  0.015625
#define POW2_M11 4.8828125e-4

#define POW2_M19 1.907348632812500e-6
#define POW2_M21 4.76837158e-7
#define POW2_M23 1.1920929e-7
#define POW2_M24 5.960464477539063e-008
#define POW2_M27 7.450580596923828e-009
#define POW2_M29 1.862645149230957e-9
#define POW2_M30 9.313225746154785e-010
#define POW2_M31 4.656612873077393e-10
#define POW2_M33 1.164153218269348e-10
#define POW2_M43 1.136868377216160e-13
#define POW2_M50 8.881784197001252e-016
#define POW2_M55 2.775557561562891e-17
#define POW2_M66 1.3552527e-20

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
// 字个数，10个
#define WORD_NUM (10)

#define  N_WORD (WORD_NUM *(SUBFRAME_NUM+1))

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
#define MAX_CHAN_SIM (12)

/*! \brief Maximum number of user motion points */
#define USER_MOTION_SIZE (864000) // for 24 hours at 10Hz

// NH码长度
#define NH_CODE_LEN (20)

#define R2D 57.2957795131
#define SPEED_OF_LIGHT 2.99792458e8
#define LAMBDA_B1 0.1920394863102765

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

typedef struct
{
	beidou_time g;
	double range; // pseudorange
	double rate;
	double d; // geometric distance
	double azel[2];
	double iono_delay;
} range_t;

// 星历
typedef struct
{
    int vflag;
    UTC_time t;
    beidou_time toc; /*!< Time of Clock */
    beidou_time toe; /*!< Time of Ephemeris */
    int aodc;	/*!< 时钟数据期龄 */
    int aode;	/*!< 星历数据期龄 */
    double deltan;	/*!< Delta-N (radians/sec) */
    double cuc;	/*!< Cuc (radians) */
    double cus;	/*!< Cus (radians) */
    double cic;	/*!< Correction to inclination cos (radians) */
    double cis;	/*!< Correction to inclination sin (radians) */
    double crc;	/*!< Correction to radius cos (meters) */
    double crs;	/*!< Correction to radius sin (meters) */
    double ecc;	/*!< e Eccentricity */
    double sqrta;	/*!< sqrt(A) (sqrt(m)) */
    double m0;	/*!< Mean anamoly (radians) */
    double omg0;	/*!< Longitude of the ascending node (radians) */
    double inc0;	/*!< Inclination (radians) */
    double omega; //GPS中的aop
    double omgdot;	/*!< Omega dot (radians/s) */
    double idot;	/*!< IDOT (radians/s) */
    double af0;	/*!< Clock offset (seconds) */
    double af1;	/*!< rate (sec/sec) */
    double af2;	/*!< acceleration (sec/sec^2) */
 
    // Working variables follow
    double n; 	/*!< Mean motion (Average angular velocity) */
    double sq1e2;	/*!< sqrt(1-e^2) */
    double A;	/*!< Semi-major axis */
    double omgkdot; /*!< OmegaDot-OmegaEdot */
    
    //北斗添加
    double sv_acc; // 卫星精度
    int sath1; // 卫星可用性0可用
    double tgd1;
    double tgd2;
} ephemeris;

typedef struct
{
    int enable;
    int vflag;
    double alpha0,alpha1,alpha2,alpha3;
    double beta0,beta1,beta2,beta3;
    double A0,A1;
    int dtls,tot,wnt;
    int dtlsf,dn,wnlsf;
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
    double azel[2];
	range_t rho0;

    int prn_code_bit;
    beidou_time g0; // 北斗时间
    /**
     * TODO 目前只模仿第一个子帧
     */
    unsigned long subframe[SUBFRAME_NUM][WORD_NUM];
    // TODO 暂时不用
    unsigned int subfranme_num; // 子帧号
    // TODO 暂时不用
    unsigned long subframe_word[3];
    // 目前只模仿第一帧
    int subframe_word_bits[N_WORD][WORD_LEN];
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
 * 星历转D1导航电文子帧
 * !!! 目前只模仿第一帧
 * @param eph
 * @param ion
 * @param sbf
 */
void eph2sbf_D1(const ephemeris eph, const ionoutc_t ion, unsigned long sbf[][WORD_NUM]);

/**
 * 星历转D2导航电文子帧
 * @param eph
 * @param ion
 * @param sbf
 */
void eph2sbf_D2(const ephemeris eph, const ionoutc_t ion, unsigned long sbf[][WORD_NUM]);

/**
 * 分配信道
 * @param chan
 * @param eph
 * @return
 */
int allocate_channel_D2(beidou_channel *chan, const ephemeris *eph, ionoutc_t ionoutc, beidou_time grx, double *xyz, double elvMask);

int allocate_channel_D1(beidou_channel *chan, const ephemeris *eph, ionoutc_t ionoutc, beidou_time grx, double *xyz, double elvMask);

/**
 * 生成导航信息
 * @param bd_time
 * @param chan
 * @param init
 * @return
 */
int nav_msg_gen_D2(beidou_time bd_time, beidou_channel *chan, int init);

int nav_msg_gen_D1(beidou_time bd_time, beidou_channel *chan, int init);

void init(beidou_channel *chan, int prn_number);

/**
 * 信号模拟
 * @param arg
 * @return
 */
void *beidou_task(void *arg);

/**
 * 读取星历内容
 * @param flag 1为D1电文，2为D2
 * @param eph
 * @param ionoutc
 * @param fname
 * @return
 */
int readRinexNavAll(int flag, ephemeris eph[][MAX_SAT], ionoutc_t *ionoutc, const char *fname);


/**
 * 计算两个时间差的秒数
 * @param bd1
 * @param bd0
 * @return
 */
double subBeidouTime(beidou_time bd1, beidou_time bd0);

int checkSatVisibility(ephemeris eph, beidou_time g, double *xyz, double elvMask, double *azel);

double normVect(const double *x);

void xyz2llh(const double *xyz, double *llh);

void llh2xyz(const double *llh, double *xyz);

void ltcmat(const double *llh, double t[3][3]);

void satpos(ephemeris eph, beidou_time g, double *pos, double *vel, double *clk);

void subVect(double *y, const double *x1, const double *x2);

void ecef2neu(const double *xyz, double t[3][3], double *neu);

void neu2azel(double *azel, const double *neu);

double dotProd(const double *x1, const double *x2);

double ionosphericDelay(const ionoutc_t *ionoutc, beidou_time g, double *llh, double *azel);

/**
 * 将星历数据转换为导航电文所要求的补码
 * @param source
 * @param bits
 * @return
 */
unsigned long convert_nav_msg_complement(long source, int bits);

#endif 