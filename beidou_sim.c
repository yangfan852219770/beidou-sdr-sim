#include<stdio.h>
#include<stdlib.h>
#include<stdbool.h>
#include<string.h>
#include<math.h>

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

/*
double ant_pat_db[37] = {
        0.00,  0.00,  0.22,  0.44,  0.67,  1.11,  1.56,  2.00,  2.44,  2.89,  3.56,  4.22,
        4.89,  5.56,  6.22,  6.89,  7.56,  8.22,  8.89,  9.78, 10.67, 11.56, 12.44, 13.33,
        14.44, 15.56, 16.67, 17.78, 18.89, 20.00, 21.33, 22.67, 24.00, 25.56, 27.33, 29.33,
        31.56
};
*/
int allocatedSat[MAX_SAT];

/*! \brief Subtract two vectors of double
 *  \param[out] y Result of subtraction
 *  \param[in] x1 Minuend of subtracion
 *  \param[in] x2 Subtrahend of subtracion
 */
void subVect(double *y, const double *x1, const double *x2)
{
    y[0] = x1[0]-x2[0];
    y[1] = x1[1]-x2[1];
    y[2] = x1[2]-x2[2];

    return;
}

/*! \brief Compute Norm of Vector
 *  \param[in] x Input vector
 *  \returns Length (Norm) of the input vector
 */
double normVect(const double *x)
{
    return(sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]));
}

/*! \brief Convert Earth-centered Earth-fixed (ECEF) into Lat/Long/Heighth
 *  \param[in] xyz Input Array of X, Y and Z ECEF coordinates
 *  \param[out] llh Output Array of Latitude, Longitude and Height
 */
void xyz2llh(const double *xyz, double *llh)
{
    double a,eps,e,e2;
    double x,y,z;
    double rho2,dz,zdz,nh,slat,n,dz_new;

    a = WGS84_RADIUS;
    e = WGS84_ECCENTRICITY;

    eps = 1.0e-3;
    e2 = e*e;

    if (normVect(xyz)<eps)
    {
        // Invalid ECEF vector
        llh[0] = 0.0;
        llh[1] = 0.0;
        llh[2] = -a;

        return;
    }

    x = xyz[0];
    y = xyz[1];
    z = xyz[2];

    rho2 = x*x + y*y;
    dz = e2*z;

    while (1)
    {
        zdz = z + dz;
        nh = sqrt(rho2 + zdz*zdz);
        slat = zdz / nh;
        n = a / sqrt(1.0-e2*slat*slat);
        dz_new = n*e2*slat;

        if (fabs(dz-dz_new) < eps)
            break;

        dz = dz_new;
    }

    llh[0] = atan2(zdz, sqrt(rho2));
    llh[1] = atan2(y, x);
    llh[2] = nh - n;

    return;
}

/*! \brief Convert Lat/Long/Height into Earth-centered Earth-fixed (ECEF)
 *  \param[in] llh Input Array of Latitude, Longitude and Height
 *  \param[out] xyz Output Array of X, Y and Z ECEF coordinates
 */
void llh2xyz(const double *llh, double *xyz)
{
    double n;
    double a;
    double e;
    double e2;
    double clat;
    double slat;
    double clon;
    double slon;
    double d,nph;
    double tmp;

    a = WGS84_RADIUS;
    e = WGS84_ECCENTRICITY;
    e2 = e*e;

    clat = cos(llh[0]);
    slat = sin(llh[0]);
    clon = cos(llh[1]);
    slon = sin(llh[1]);
    d = e*slat;

    n = a/sqrt(1.0-d*d);
    nph = n + llh[2];

    tmp = nph*clat;
    xyz[0] = tmp*clon;
    xyz[1] = tmp*slon;
    xyz[2] = ((1.0-e2)*n + llh[2])*slat;

    return;
}

/*! \brief Compute the intermediate matrix for LLH to ECEF
 *  \param[in] llh Input position in Latitude-Longitude-Height format
 *  \param[out] t Three-by-Three output matrix
 */
void ltcmat(const double *llh, double t[3][3])
{
    double slat, clat;
    double slon, clon;

    slat = sin(llh[0]);
    clat = cos(llh[0]);
    slon = sin(llh[1]);
    clon = cos(llh[1]);

    t[0][0] = -slat*clon;
    t[0][1] = -slat*slon;
    t[0][2] = clat;
    t[1][0] = -slon;
    t[1][1] = clon;
    t[1][2] = 0.0;
    t[2][0] = clat*clon;
    t[2][1] = clat*slon;
    t[2][2] = slat;

    return;
}

/*! \brief Convert Earth-centered Earth-Fixed to ?
 *  \param[in] xyz Input position as vector in ECEF format
 *  \param[in] t Intermediate matrix computed by \ref ltcmat
 *  \param[out] neu Output position as North-East-Up format
 */
void ecef2neu(const double *xyz, double t[3][3], double *neu)
{
    neu[0] = t[0][0]*xyz[0] + t[0][1]*xyz[1] + t[0][2]*xyz[2];
    neu[1] = t[1][0]*xyz[0] + t[1][1]*xyz[1] + t[1][2]*xyz[2];
    neu[2] = t[2][0]*xyz[0] + t[2][1]*xyz[1] + t[2][2]*xyz[2];

    return;
}

/*! \brief Convert North-Eeast-Up to Azimuth + Elevation
 *  \param[in] neu Input position in North-East-Up format
 *  \param[out] azel Output array of azimuth + elevation as double
 */
void neu2azel(double *azel, const double *neu)
{
    double ne;

    azel[0] = atan2(neu[1],neu[0]);
    if (azel[0]<0.0)
        azel[0] += (2.0*PI);

    ne = sqrt(neu[0]*neu[0] + neu[1]*neu[1]);
    azel[1] = atan2(neu[2], ne);

    return;
}

/*! \brief Compute Satellite position, velocity and clock at given time
 *  \param[in] eph Ephemeris data of the satellite
 *  \param[in] g GPS time at which position is to be computed
 *  \param[out] pos Computed position (vector)
 *  \param[out] vel Computed velociy (vector)
 *  \param[clk] clk Computed clock
 */
void satpos(ephemeris eph, beidou_time g, double *pos, double *vel, double *clk)
{
    // Computing Satellite Velocity using the Broadcast Ephemeris
    // http://www.ngs.noaa.gov/gps-toolbox/bc_velo.htm

    double tk;
    double mk;
    double ek;
    double ekold;
    double ekdot;
    double cek,sek;
    double pk;
    double pkdot;
    double c2pk,s2pk;
    double uk;
    double ukdot;
    double cuk,suk;
    double ok;
    double sok,cok;
    double ik;
    double ikdot;
    double sik,cik;
    double rk;
    double rkdot;
    double xpk,ypk;
    double xpkdot,ypkdot;

    double relativistic, OneMinusecosE, tmp;

    tk = g.second - eph.toe.second;

    if(tk>SECONDS_IN_HALF_WEEK)
        tk -= SECONDS_IN_WEEK;
    else if(tk<-SECONDS_IN_HALF_WEEK)
        tk += SECONDS_IN_WEEK;

    mk = eph.m0 + eph.n*tk;
    ek = mk;
    ekold = ek + 1.0;

    OneMinusecosE = 0; // Suppress the uninitialized warning.
    while(fabs(ek-ekold)>1.0E-14)
    {
        ekold = ek;
        OneMinusecosE = 1.0-eph.ecc*cos(ekold);
        ek = ek + (mk-ekold+eph.ecc*sin(ekold))/OneMinusecosE;
    }

    sek = sin(ek);
    cek = cos(ek);

    ekdot = eph.n/OneMinusecosE;

    relativistic = -4.442807633E-10*eph.ecc*eph.sqrta*sek;

    pk = atan2(eph.sq1e2*sek,cek-eph.ecc) + eph.omega;
    pkdot = eph.sq1e2*ekdot/OneMinusecosE;

    s2pk = sin(2.0*pk);
    c2pk = cos(2.0*pk);

    uk = pk + eph.cus*s2pk + eph.cuc*c2pk;
    suk = sin(uk);
    cuk = cos(uk);
    ukdot = pkdot*(1.0 + 2.0*(eph.cus*c2pk - eph.cuc*s2pk));

    rk = eph.A*OneMinusecosE + eph.crc*c2pk + eph.crs*s2pk;
    rkdot = eph.A*eph.ecc*sek*ekdot + 2.0*pkdot*(eph.crs*c2pk - eph.crc*s2pk);

    ik = eph.inc0 + eph.idot*tk + eph.cic*c2pk + eph.cis*s2pk;
    sik = sin(ik);
    cik = cos(ik);
    ikdot = eph.idot + 2.0*pkdot*(eph.cis*c2pk - eph.cic*s2pk);

    xpk = rk*cuk;
    ypk = rk*suk;
    xpkdot = rkdot*cuk - ypk*ukdot;
    ypkdot = rkdot*suk + xpk*ukdot;

    ok = eph.omg0 + tk*eph.omgkdot - OMEGA_EARTH*eph.toe.second;
    sok = sin(ok);
    cok = cos(ok);

    pos[0] = xpk*cok - ypk*cik*sok;
    pos[1] = xpk*sok + ypk*cik*cok;
    pos[2] = ypk*sik;

    tmp = ypkdot*cik - ypk*sik*ikdot;

    vel[0] = -eph.omgkdot*pos[1] + xpkdot*cok - tmp*sok;
    vel[1] = eph.omgkdot*pos[0] + xpkdot*sok + tmp*cok;
    vel[2] = ypk*cik*ikdot + ypkdot*sik;

    // Satellite clock correction
    tk = g.second - eph.toc.second;

    if(tk>SECONDS_IN_HALF_WEEK)
        tk -= SECONDS_IN_WEEK;
    else if(tk<-SECONDS_IN_HALF_WEEK)
        tk += SECONDS_IN_WEEK;

    clk[0] = eph.af0 + tk*(eph.af1 + tk*eph.af2) + relativistic - eph.tgd1;
    clk[1] = eph.af1 + 2.0*tk*eph.af2;

    return;
}

double subBeidouTime(beidou_time bd1, beidou_time bd0)
{
    double dt;

    dt = bd1.second - bd0.second;
    dt += (double)(bd1.week - bd0.week) * SECONDS_IN_WEEK;

    return(dt);
}

beidou_time inBeidouTime(beidou_time g0, double dt)
{
    beidou_time g1;

    g1.week = g0.week;
    g1.second = g0.second + dt;

    g1.second = round(g1.second*1000.0)/1000.0; // Avoid rounding error

    while (g1.second>=SECONDS_IN_WEEK)
    {
        g1.second -= SECONDS_IN_WEEK;
        g1.week++;
    }

    while (g1.second<0.0)
    {
        g1.second += SECONDS_IN_WEEK;
        g1.week--;
    }

    return(g1);
}

// 读取星历
int readRinexNavAll(ephemeris eph[][MAX_SAT], ionoutc_t *ionoutc, const char *fname)
{
    FILE *fp;
    int ieph;

    int sv;
    char str[MAX_CHAR];
    char tmp[20];

    UTC_time t;
    beidou_time bd;
    beidou_time bd0;
    double dt;

    // 读取头文件标识
    int flags = 0x0;

    if (NULL==(fp=fopen(fname, "rt")))
        return(-1);

    //清除有效标识
    for (ieph=0; ieph<EPHEM_ARRAY_SIZE; ieph++)
        for (sv=0; sv<MAX_SAT; sv++)
            eph[ieph][sv].vflag = 0;

    //TODO 读取头文件内容，目前只读一组即可
    while(1){
        if (NULL==fgets(str, MAX_CHAR, fp))
            break;
        // 前有60个空格，到达头文件尾
        if (strncmp(str+60, "END OF HEADER", 13)==0)
            break;
            // 读取alpha
        else if(strncmp(str, "BDSA", 4)==0){
            // alpha0
            strncpy(tmp, str+5, 12);
            // 0 将结束符添加到字符串的末尾 '\0'与0的ASCII码的值都是0
            tmp[12] = 0;
            /*atof()函数将字符转换为double
			可跳过前面的空白字符
			遇到非字符或'\0'结束
			*/
            ionoutc->alpha0 = atof(tmp);

            // alpha1
            strncpy(tmp, str+17, 12);
            tmp[12] = 0;

            ionoutc->alpha1 = atof(tmp);

            // alpha2
            strncpy(tmp, str+29, 12);
            tmp[12] = 0;
            ionoutc->alpha2 = atof(tmp);

            // alpha3
            strncpy(tmp, str+41, 12);
            tmp[12] = 0;
            ionoutc->alpha3 = atof(tmp);

            flags |= 0x1;
        }
            // 读取beta
        else if(strncmp(str, "BDSB", 4)==0){
            //beta0
            strncpy(tmp, str+5, 12);
            tmp[12] = 0;
            ionoutc->beta0 = atof(tmp);

            //beta1
            strncpy(tmp, str+17, 12);
            tmp[12] = 0;
            ionoutc->beta1 = atof(tmp);

            //beta2
            strncpy(tmp, str+29, 12);
            ionoutc->beta2 = atof(tmp);

            //beta3
            strncpy(tmp, str+41, 12);
            tmp[12] = 0;
            ionoutc->beta3 = atof(tmp);

            flags |= 0x1<<1;
        }
        // TODO 先有一组数据即可


        //TODO 头文件中BDS1、BDS2 、BDS3的行
    }

    ionoutc->vflag = false;
    // 读完电离层参数
    if (flags==0x3)
        ionoutc->vflag = true;

    bd0.week = -1;
    ieph = 0;

    // 读取导航电文
    // TODO 当前只读取GEO卫星1-5
    while(1){
        if (NULL==fgets(str, MAX_CHAR, fp))
            break;
        // PRN号
        strncpy(tmp, str+1, 2);
        tmp[2] = 0;
        // atoi()函数，将字符串转为int
        sv = atoi(tmp)-1;
        // TODO 只读GEO卫星1-5
        if(sv >= 5){
            // 跳过后面的7行
            fgets(str, MAX_CHAR, fp);
            fgets(str, MAX_CHAR, fp);
            fgets(str, MAX_CHAR, fp);
            fgets(str, MAX_CHAR, fp);
            fgets(str, MAX_CHAR, fp);
            fgets(str, MAX_CHAR, fp);
            fgets(str, MAX_CHAR, fp);
            continue;
        }
        // 星历的时间信息
        strncpy(tmp, str+4, 4);
        tmp[4] = 0;
        t.year = atoi(tmp);

        strncpy(tmp, str+9, 2);
        tmp[2] = 0;
        t.month = atoi(tmp);

        strncpy(tmp, str+12, 2);
        tmp[2] = 0;
        t.day = atoi(tmp);

        strncpy(tmp, str+15, 2);
        tmp[2] = 0;
        t.hour = atoi(tmp);

        strncpy(tmp, str+18, 2);
        tmp[2] = 0;
        t.minute = atoi(tmp);

        strncpy(tmp, str+21, 2);
        tmp[2] = 0;
        t.second = atoi(tmp);

        //将UTC时间转换为北斗时间
        UTC2beidou_time(&t, &bd);

        if(bd0.week == -1)
            bd0 = bd;

        //检查当前时钟
        dt = subBeidouTime(bd, bd0);

        if (dt>SECONDS_IN_HOUR)
        {
            bd0 = bd;
            ieph++; // a new set of ephemerides

            if (ieph>=EPHEM_ARRAY_SIZE)
                break;
        }
        // 设置星历时间
        eph[ieph][sv].t = t;

        //设置卫星时钟
        eph[ieph][sv].toc = bd;

        // af0,钟差
        strncpy(tmp, str+23, 19);
        tmp[19] = 0;
        eph[ieph][sv].af0 = atof(tmp);

        //钟漂
        strncpy(tmp, str+42, 19);
        tmp[19] = 0;
        eph[ieph][sv].af1 = atof(tmp);

        //钟漂的速率
        strncpy(tmp, str+61, 19);
        tmp[19] = 0;
        eph[ieph][sv].af2 = atof(tmp);

        if (NULL==fgets(str, MAX_CHAR, fp))
            break;

        // 轨道参数1
        // 星历数据期龄
        strncpy(tmp, str+4, 19);
        tmp[19] = 0;
        eph[ieph][sv].aode = atof(tmp);

        // 轨道半径的正弦调和改正项的振幅
        strncpy(tmp, str+23, 19);
        tmp[19] = 0;
        eph[ieph][sv].crs = atof(tmp);

        // 卫星平均运动速率与计算值之差
        strncpy(tmp, str+42, 19);
        tmp[19] = 0;
        eph[ieph][sv].deltan = atof(tmp);

        // 参考时间的平近点角
        strncpy(tmp, str+61, 19);
        tmp[19] = 0;
        eph[ieph][sv].m0 = atof(tmp);

        // 轨道参数2
        if (NULL==fgets(str, MAX_CHAR, fp))
            break;

        // 纬度幅角的余弦调和改正项的振幅
        strncpy(tmp, str+4, 19);
        tmp[19] = 0;
        eph[ieph][sv].cuc = atof(tmp);

        // 偏心率
        strncpy(tmp, str+23, 19);
        tmp[19] = 0;
        eph[ieph][sv].ecc = atof(tmp);

        // 纬度幅角的正弦调和改正项的振幅
        strncpy(tmp, str+42, 19);
        tmp[19] = 0;
        eph[ieph][sv].cus = atof(tmp);

        // 长半轴的平方根
        strncpy(tmp, str+61, 19);
        tmp[19] = 0;
        eph[ieph][sv].sqrta = atof(tmp);

        // 轨道参数3
        if (NULL==fgets(str, MAX_CHAR, fp))
            break;

        // 星历参考时刻
        strncpy(tmp, str+4, 19);
        tmp[19] = 0;
        eph[ieph][sv].toe.second = atof(tmp);

        // 轨道倾角的余弦调和改正项的振幅
        strncpy(tmp, str+23, 19);
        tmp[19] = 0;
        eph[ieph][sv].cic = atof(tmp);

        // 升交点经度
        strncpy(tmp, str+42, 19);
        tmp[19] = 0;
        eph[ieph][sv].omg0 = atof(tmp);

        // 轨道倾角的正弦调和改正项的振幅
        strncpy(tmp, str+61, 19);
        tmp[19] = 0;
        eph[ieph][sv].cis = atof(tmp);

        // 轨道参数4

        if (NULL==fgets(str, MAX_CHAR, fp))
            break;

        // i0 参考时间的轨道倾角
        strncpy(tmp, str+4, 19);
        tmp[19] = 0;
        eph[ieph][sv].inc0 = atof(tmp);

        // 轨道半径的余弦调和改正项的振幅
        strncpy(tmp, str+23, 19);
        tmp[19] = 0;
        eph[ieph][sv].crc = atof(tmp);

        // 近地点幅角
        strncpy(tmp, str+42, 19);
        tmp[19] = 0;
        eph[ieph][sv].omega = atof(tmp);

        // 升交点赤经变化率
        strncpy(tmp, str+61, 19);
        tmp[19] = 0;
        eph[ieph][sv].omgdot = atof(tmp);

        // 轨道参数5
        if (NULL==fgets(str, MAX_CHAR, fp))
            break;

        // IDOT 轨道倾角变化率
        strncpy(tmp, str+4, 19);
        tmp[19] = 0;
        eph[ieph][sv].idot = atof(tmp);

        // BDT 周数
        strncpy(tmp, str+42, 19);
        tmp[19] = 0;
        eph[ieph][sv].toe.week = atoi(tmp);

        // 轨道参数6
        if (NULL==fgets(str, MAX_CHAR, fp))
            break;

        // 轨道精度
        strncpy(tmp, str+4, 19);
        tmp[19] = 0;
        eph[ieph][sv].sv_acc = atof(tmp);

        // 卫星可用标识，0为可用
        strncpy(tmp, str+23, 19);
        tmp[19] = 0;
        eph[ieph][sv].sath1 = atoi(tmp);

        // 设备时延差
        strncpy(tmp, str+42, 19);
        tmp[19] = 0;
        eph[ieph][sv].tgd1 = atof(tmp);

        strncpy(tmp, str+61, 19);
        tmp[19] = 0;
        eph[ieph][sv].tgd2 = atof(tmp);

        // 轨道参数7
        if (NULL==fgets(str, MAX_CHAR, fp))
            break;

        // 时钟数据期龄
        strncpy(tmp, str+23, 19);
        tmp[19] = 0;
        eph[ieph][sv].aodc = atoi(tmp);

        eph[ieph][sv].vflag = 1;

        eph[ieph][sv].A = eph[ieph][sv].sqrta * eph[ieph][sv].sqrta;
        eph[ieph][sv].n = sqrt(GM_EARTH/(eph[ieph][sv].A*eph[ieph][sv].A*eph[ieph][sv].A)) + eph[ieph][sv].deltan;
        eph[ieph][sv].sq1e2 = sqrt(1.0 - eph[ieph][sv].ecc*eph[ieph][sv].ecc);
        eph[ieph][sv].omgkdot = eph[ieph][sv].omgdot - OMEGA_EARTH;
    }

    fclose(fp);

    if (bd0.week>=0)
        ieph += 1;

    return(ieph);
}

double ionosphericDelay(const ionoutc_t *ionoutc, beidou_time g, double *llh, double *azel)
{
    double iono_delay = 0.0;
    double E,phi_u,lam_u,F;

    if (ionoutc->enable==false)
        return (0.0); // No ionospheric delay

    E = azel[1]/PI;
    phi_u = llh[0]/PI;
    lam_u = llh[1]/PI;

    // Obliquity factor
    F = 1.0 + 16.0*pow((0.53 - E),3.0);

    if (ionoutc->vflag==false)
        iono_delay = F*5.0e-9*SPEED_OF_LIGHT;
    else
    {
        double t,psi,phi_i,lam_i,phi_m,phi_m2,phi_m3;
        double AMP,PER,X,X2,X4;

        // Earth's central angle between the user position and the earth projection of
        // ionospheric intersection point (semi-circles)
        psi = 0.0137/(E + 0.11) - 0.022;

        // Geodetic latitude of the earth projection of the ionospheric intersection point
        // (semi-circles)
        phi_i = phi_u + psi*cos(azel[0]);
        if(phi_i>0.416)
            phi_i = 0.416;
        else if(phi_i<-0.416)
            phi_i = -0.416;

        // Geodetic longitude of the earth projection of the ionospheric intersection point
        // (semi-circles)
        lam_i = lam_u + psi*sin(azel[0])/cos(phi_i*PI);

        // Geomagnetic latitude of the earth projection of the ionospheric intersection
        // point (mean ionospheric height assumed 350 km) (semi-circles)
        phi_m = phi_i + 0.064*cos((lam_i - 1.617)*PI);
        phi_m2 = phi_m*phi_m;
        phi_m3 = phi_m2*phi_m;

        AMP = ionoutc->alpha0 + ionoutc->alpha1*phi_m
              + ionoutc->alpha2*phi_m2 + ionoutc->alpha3*phi_m3;
        if (AMP<0.0)
            AMP = 0.0;

        PER = ionoutc->beta0 + ionoutc->beta1*phi_m
              + ionoutc->beta2*phi_m2 + ionoutc->beta3*phi_m3;
        if (PER<72000.0)
            PER = 72000.0;

        // Local time (sec)
        t = SECONDS_IN_DAY/2.0*lam_i + g.second;
        while(t>=SECONDS_IN_DAY)
            t -= SECONDS_IN_DAY;
        while(t<0)
            t += SECONDS_IN_DAY;

        // Phase (radians)
        X = 2.0*PI*(t - 50400.0)/PER;

        if(fabs(X)<1.57)
        {
            X2 = X*X;
            X4 = X2*X2;
            iono_delay = F*(5.0e-9 + AMP*(1.0 - X2/2.0 + X4/24.0))*SPEED_OF_LIGHT;
        }
        else
            iono_delay = F*5.0e-9*SPEED_OF_LIGHT;
    }

    return (iono_delay);
}

/*! \brief Compute range between a satellite and the receiver
 *  \param[out] rho The computed range
 *  \param[in] eph Ephemeris data of the satellite
 *  \param[in] g GPS time at time of receiving the signal
 *  \param[in] xyz position of the receiver
 */
void computeRange(range_t *rho, ephemeris eph, ionoutc_t *ionoutc, beidou_time g, double xyz[])
{
    double pos[3],vel[3],clk[2];
    double los[3];
    double tau;
    double range,rate;
    double xrot,yrot;

    double llh[3],neu[3];
    double tmat[3][3];

    // SV position at time of the pseudorange observation.
    satpos(eph, g, pos, vel, clk);

    // Receiver to satellite vector and light-time.
    subVect(los, pos, xyz);
    tau = normVect(los)/SPEED_OF_LIGHT;

    // Extrapolate the satellite position backwards to the transmission time.
    pos[0] -= vel[0]*tau;
    pos[1] -= vel[1]*tau;
    pos[2] -= vel[2]*tau;

    // Earth rotation correction. The change in velocity can be neglected.
    xrot = pos[0] + pos[1]*OMEGA_EARTH*tau;
    yrot = pos[1] - pos[0]*OMEGA_EARTH*tau;
    pos[0] = xrot;
    pos[1] = yrot;

    // New observer to satellite vector and satellite range.
    subVect(los, pos, xyz);
    range = normVect(los);
    rho->d = range;

    // Pseudorange.
    rho->range = range - SPEED_OF_LIGHT*clk[0];

    // Relative velocity of SV and receiver.
    rate = dotProd(vel, los)/range;

    // Pseudorange rate.
    rho->rate = rate; // - SPEED_OF_LIGHT*clk[1];

    // Time of application.
    rho->g = g;

    // Azimuth and elevation angles.
    xyz2llh(xyz, llh);
    ltcmat(llh, tmat);
    ecef2neu(los, tmat, neu);
    neu2azel(rho->azel, neu);

    // Add ionospheric delay
    rho->iono_delay = ionosphericDelay(ionoutc, g, llh, rho->azel);
    rho->range += rho->iono_delay;

    return;
}

/*! \brief Compute the code phase for a given channel (satellite)
 *  \param chan Channel on which we operate (is updated)
 *  \param[in] rho1 Current range, after \a dt has expired
 *  \param[in dt delta-t (time difference) in seconds
 */
void computeCodePhase(beidou_channel *chan, range_t rho1, double dt)
{
    double ms;
    int ims;
    double rhorate;

    // Pseudorange rate.
    rhorate = (rho1.range - chan->rho0.range)/dt;

    // Carrier and code frequency.
    chan->f_carr = -rhorate/LAMBDA_B1;
    chan->f_code = PRN_CODE_FREQ + chan->f_carr*CARR_TO_PRN;

    // Initial code phase and data bit counters.
    ms = ((subBeidouTime(chan->rho0.g,chan->g0)+0.6) - chan->rho0.range/SPEED_OF_LIGHT)*1000.0;

    ims = (int)ms;
    chan->code_phase = (ms-(double)ims)*PRN_SEQ_LEN; // in chip

    chan->iword = ims/60; // 1 word = 30 bits = 60 ms

    ims -= chan->iword*60;

    chan->ibit = ims/2; // 1 bit = 2 code = 2 ms
    ims -= chan->ibit*2;

    chan->icode = ims; // 1 code = 1 ms

    chan->prn_code_bit = chan->prn_code[(int)chan->code_phase] * 2 - 1;
    chan->data_bit = chan->subframe_word_bits[chan->iword][chan->ibit] * 2  - 1;

    // Save current pseudorange
    chan->rho0 = rho1;

    return;
}

/*! \brief Compute dot-product of two vectors
 *  \param[in] x1 First multiplicand
 *  \param[in] x2 Second multiplicand
 *  \returns Dot-product of both multiplicands
 */
double dotProd(const double *x1, const double *x2)
{
    return(x1[0]*x2[0]+x1[1]*x2[1]+x1[2]*x2[2]);
}

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

    // 一个帧中的第一个字，第一个字的前15bits不做校验，只对后11bits做校验
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

unsigned long convert_nav_msg_complement(long source, int bits){
    unsigned long mask;
    unsigned long result;
    // 符号位
    unsigned long sign_bit;

    if(0 == source)
        return 0;

    // 正数不用转换，直接操作符号位即可。0为正，1为负。
    if(source > 0){
        result = (0 << (bits-1)) & source;
    }
    else{
        sign_bit = 1 << (bits-1);

        // 数据位
        mask = 1;
        for(int i = 0; i < bits - 2; ++i)
            mask = (mask << 1) | 0x1UL;

        // 结果的数据位
        result = source & mask;
        // 符号位
        result = sign_bit | result;
    }

    //printf("转换后的结果为:%ld\n",result);

    return result;
}

void eph2sbf_D2(const ephemeris eph, const ionoutc_t ionoutc, unsigned long sbf[][WORD_NUM]){

    unsigned long wn;
    unsigned long toe;
    unsigned long toc;
    unsigned long aode;
    unsigned long aodc;
    long deltan;
    long cuc;
    long cus;
    long cic;
    long cis;
    long crc;
    long crs;
    unsigned long ecc;
    unsigned long sqrta;
    long m0;
    long omg0;
    long inc0;
    long omega;
    long omgdot;
    long idot;
    long af0;
    long af1;
    long af2;
    long tgd1,tgd2;


    unsigned long ura = 0UL;

    unsigned long wna;
    unsigned long toa;

    signed long alpha0,alpha1,alpha2,alpha3;
    signed long beta0,beta1,beta2,beta3;
    signed long A0,A1;
    signed long dtls,dtlsf;
    unsigned long tot,wnt,wnlsf,dn;

    // 转换后的补码
    unsigned long deltan_c;
    unsigned long cuc_c;
    unsigned long cus_c;
    unsigned long cic_c;
    unsigned long cis_c;
    unsigned long crc_c;
    unsigned long crs_c;

    unsigned long m0_c;
    unsigned long omg0_c;
    unsigned long inc0_c;
    unsigned long omega_c;
    unsigned long omgdot_c;
    unsigned long idot_c;
    unsigned long af0_c;
    unsigned long af1_c;
    unsigned long af2_c;
    unsigned long tgd1_c,tgd2_c;

    unsigned long alpha0_c,alpha1_c,alpha2_c,alpha3_c;
    unsigned long beta0_c,beta1_c,beta2_c,beta3_c;

    unsigned long second;
    unsigned long week;

    // 帧同步码
    unsigned long pre = 0x712;

    // TODO BDT的周内秒和周数
    second = (unsigned long)(eph.toc.second);
    week = (unsigned long)(eph.toc.week);

    // 有缩放因子，所以需要除
    // TODO modify
    toe = (unsigned long)(eph.toe.second/8.0);
    toc = (unsigned long)(eph.toc.second/8.0);

    aode = (unsigned long)(eph.aode);
    aodc = (unsigned long)(eph.aodc);
    // 星历中deltan的单位是rad/s，导航电文是pi/s
    deltan = (long)(eph.deltan/POW2_M43/PI);
    cuc = (long)(eph.cuc/POW2_M31);
    cus = (long)(eph.cus/POW2_M31);
    cic = (long)(eph.cic/POW2_M31);
    cis = (long)(eph.cis/POW2_M31);
    crc = (long)(eph.crc/POW2_M6);
    crs = (long)(eph.crs/POW2_M6);
    ecc = (unsigned long)(eph.ecc/POW2_M33);
    sqrta = (unsigned long)(eph.sqrta/POW2_M19);
    m0 = (long)(eph.m0/POW2_M31/PI);
    omg0 = (long)(eph.omg0/POW2_M31/PI);
    inc0 = (long)(eph.inc0/POW2_M31/PI);
    omega = (long)(eph.omega/POW2_M31/PI);
    omgdot = (long)(eph.omgdot/POW2_M43/PI);
    idot = (long)(eph.idot/POW2_M43/PI);
    af0 = (long)(eph.af0/POW2_M33);
    af1 = (long)(eph.af1/POW2_M50);
    af2 = (long)(eph.af2/POW2_M66);
    tgd1 = (long)(eph.tgd1/0.1);
    tgd2 = (long)(eph.tgd2/0.1);

    // 电离层改正数
    alpha0 = (signed long)round(ionoutc.alpha0/POW2_M30);
    alpha1 = (signed long)round(ionoutc.alpha1/POW2_M27);
    alpha2 = (signed long)round(ionoutc.alpha2/POW2_M24);
    alpha3 = (signed long)round(ionoutc.alpha3/POW2_M24);
    beta0 = (signed long)round(ionoutc.beta0/2048.0);
    beta1 = (signed long)round(ionoutc.beta1/16384.0);
    beta2 = (signed long)round(ionoutc.beta2/65536.0);
    beta3 = (signed long)round(ionoutc.beta3/65536.0);

    // 将数据转换为导航电文所要求的补码格式
    af0_c = convert_nav_msg_complement(af0,24);
    af1_c = convert_nav_msg_complement(af1, 22);
    af2_c = convert_nav_msg_complement(af2, 11);
    tgd1_c = convert_nav_msg_complement(tgd1, 10);
    tgd2_c = convert_nav_msg_complement(tgd2, 10);
    omega_c = convert_nav_msg_complement(omega, 32);
    deltan_c = convert_nav_msg_complement(deltan,16);
    m0_c = convert_nav_msg_complement(m0, 32);
    omg0_c = convert_nav_msg_complement(omg0, 32);
    omgdot_c = convert_nav_msg_complement(omgdot, 24);
    inc0_c = convert_nav_msg_complement(inc0, 32);
    idot_c = convert_nav_msg_complement(idot, 14);
    cuc_c = convert_nav_msg_complement(cuc, 18);
    cus_c = convert_nav_msg_complement(cus, 18);
    crc_c = convert_nav_msg_complement(crc, 18);
    crs_c = convert_nav_msg_complement(crs, 18);
    cic_c = convert_nav_msg_complement(cic, 18);
    cis_c = convert_nav_msg_complement(cis, 18);

    alpha0_c = convert_nav_msg_complement(alpha0, 8);
    alpha1_c = convert_nav_msg_complement(alpha1, 8);
    alpha2_c = convert_nav_msg_complement(alpha2, 8);
    alpha3_c = convert_nav_msg_complement(alpha3, 8);
    beta0_c = convert_nav_msg_complement(beta0, 8);
    beta1_c = convert_nav_msg_complement(beta1, 8);
    beta2_c = convert_nav_msg_complement(beta0, 8);
    beta3_c = convert_nav_msg_complement(beta3, 8);

    week = week & 0x1FFFUL;
    second = second & 0xFFFFFUL;
    aodc = aodc & 0x1FUL;
    toc = toc & 0x1FFFFUL;
    toe = toe & 0x1FFFFUL;
    aode = aode & 0x1FUL;

    ////////////////////////////////////////////////////////////
    // 第一帧 第一字
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
     * SatH1(1bit)，卫星健康表示，0表示可用，1表示不可用p
     * AODC(5bit)，时钟数据龄期，暂定为 0
     */
    sbf[0][1] = ( second & 0xFFF) << 18 | 0x1UL << 14 | 0 << 13 | aodc << 8;
    /** 第三个字
     * URAI | WN | toc
     * URAI(4bit)，用户距离精度指数，暂定为 8
     * WN(13bit)，整周计数
     * toc(前5bit)，时钟时间
     */
    sbf[0][2] = 0x8UL << 26 | week << 13 | (toc >> 12) << 8 ;
    /**
    * 第四个字
    * toc | TGD1
    * toc(后12bit)，时钟时间
    * TGD1(10bit) 星上设备时延差
    */
    sbf[0][3] = (toc & 0xFFF) << 18 | tgd1 << 8;
    /**
     * 第五个字
     * TGD2 | Rev
     * TGD2(后6bit) 星上设备时延差
     * Rev 保留字
     */
    sbf[0][4] = tgd2 << 20 | 0 << 8;

    ////////////////////////////////////////////////////////////
    // 第一帧 第二字
    ////////////////////////////////////////////////////////////
    sbf[0][5] = pre << 19 | 0 << 15 | 0x1UL << 12 | (second >> 12) << 4;
    sbf[0][6] = ( second & 0xFFF) << 18 | 0x2UL << 14 | (alpha0_c >> 2) << 8;
    sbf[0][7] = (alpha0_c & 0x2UL) << 28 | alpha1_c << 20 | alpha2_c << 12 | (alpha3_c >> 4) << 8;
    sbf[0][8] = (alpha3_c & 0x4UL) << 26 | beta0_c << 18 | beta1_c << 10 | (beta2_c >> 6) << 8;
    sbf[0][9] = (beta2_c & 0x6UL) << 24 | beta3_c << 16;

    ////////////////////////////////////////////////////////////
    // 第一帧 第三字
    ////////////////////////////////////////////////////////////
    sbf[1][0] = pre << 19 | 0 << 15 | 0x1UL << 12 | (second >> 12) << 4;
    sbf[1][1] = ( second & 0xFFF) << 18 | 0x3UL << 14 | 0 << 8;
    sbf[1][2] = 0 << 8;
    sbf[1][3] = 0 << 20 | (af0_c >> 12) << 8;
    sbf[1][4] = (af0_c & 0xFFFUL) << 18 | (af1_c >> 18) << 14;

    ////////////////////////////////////////////////////////////
    // 第一帧 第四字
    ////////////////////////////////////////////////////////////
    sbf[1][5] = pre << 19 | 0 << 15 | 0x1UL << 12 | (second >> 12) << 4;
    sbf[1][6] = ( second & 0xFFFUL) << 18 | 0x4UL << 14 | ((af1_c & 0x3F000) << 8);
    sbf[1][7] = (af1_c & 0xFFFUL)<< 18 | (af2_c >> 1) << 8;
    sbf[1][8] = (af2_c & 0x1UL) << 29 | aode << 24 | deltan_c << 8;
    sbf[1][9] = (cuc_c >> 4) << 16 | 0 << 8;

    ////////////////////////////////////////////////////////////
    // 第一帧 第五字
    ////////////////////////////////////////////////////////////
    sbf[2][0] = pre << 19 | 0 << 15 | 0x1UL << 12 | (second >> 12) << 4;
    sbf[2][1] = ( second & 0xFFFUL) << 18 | 0x5UL << 14 | (cuc_c & 0x4UL) << 10 | (m0_c >> 30) << 8;
    sbf[2][2] = (m0_c & 0x3FFFFF00UL) << 8;
    sbf[2][3] =(m0_c & 0xFFUL) << 22 | (cus_c >> 4) << 8;
    sbf[2][4] = (cus_c & 0xF) << 26 | (ecc >> 22) << 8;

    ////////////////////////////////////////////////////////////
    // 第一帧 第六字
    ////////////////////////////////////////////////////////////
    sbf[2][5] = pre << 19 | 0 << 15 | 0x1UL << 12 | (second >> 12) << 4;
    sbf[2][6] = ( second & 0xFFFUL) << 18 | 0x6UL << 14 | (ecc & 0x3F0000UL) << 8;
    sbf[2][7] = (ecc & 0xFFFFUL) << 14 | (sqrta >> 22) << 8;
    sbf[2][8] = (sqrta & 0x3FFFFF0UL) << 8;
    sbf[2][9] = (sqrta & 0xFUL) << 26 | (cic_c >> 8) << 8;

    ////////////////////////////////////////////////////////////
    // 第一帧 第七字
    ////////////////////////////////////////////////////////////
    sbf[3][0] = pre << 19 | 0 << 15 | 0x1UL << 12 | (second >> 12) << 4;
    sbf[3][1] = ( second & 0xFFFUL) << 18 | 0x7UL << 14 | (cic_c & 0xFCUL) << 8;
    sbf[3][2] = (cic_c & 0x2UL) << 28 | cis_c << 10 | (toe >> 15) << 8;
    sbf[3][3] = (toe & 0x7FFFUL) << 15 | (inc0_c >> 25) << 8;
    sbf[3][4] = (inc0_c & 0x1FFF800UL) << 16 | 0 << 8;

    ////////////////////////////////////////////////////////////
    // 第一帧 第八字
    ////////////////////////////////////////////////////////////
    sbf[3][5] = pre << 19 | 0 << 15 | 0x1UL << 12 | (second >> 12) << 4;
    sbf[3][6] = ( second & 0xFFFUL) << 18 | 0x8UL << 14 | (inc0_c & 0x7E0UL) << 8;
    sbf[3][7] = (inc0_c & 0x1FUL) << 25 | (crc_c >> 1) << 8;
    sbf[3][8] = (crc_c & 0x1UL) << 29 | crs_c << 11 | (omgdot_c >> 21) << 8;
    sbf[3][9] = (omgdot_c & 0x1FFFE0UL ) << 14 | 0 << 8;

    ////////////////////////////////////////////////////////////
    // 第一帧 第九字
    ////////////////////////////////////////////////////////////
    sbf[4][0] = pre << 19 | 0 << 15 | 0x1UL << 12 | (second >> 12) << 4;
    sbf[4][1] = ( second & 0xFFFUL) << 18 | 0x9UL << 14 | (omgdot_c & 0x1FUL) << 9 | (omg0_c >> 31) << 8;
    sbf[4][2] = (omg0_c & 0x7FFFFE00UL) << 8;
    sbf[4][3] = (omg0_c & 0x1FFUL) << 21 | (omega_c >> 19) << 8;
    sbf[4][4] = (omega_c & 0x7FFE0UL) << 16 | 0 << 8;

    ////////////////////////////////////////////////////////////
    // 第一帧 第十字
    ////////////////////////////////////////////////////////////
    sbf[4][5] = pre << 19 | 0 << 15 | 0x1UL << 12 | (second >> 12) << 4;
    sbf[4][6] = ( second & 0xFFFUL) << 18 | 0xAUL << 14 | (omega_c & 0x1FUL) << 9 | (idot_c >> 13) << 8;
    sbf[4][7] = (idot_c & 0x1FFFUL) << 17 | 0 << 8;
    sbf[4][8] = 0UL;
    sbf[4][9] = 0UL;
}

int checkSatVisibility(ephemeris eph, beidou_time g, double *xyz, double elvMask, double *azel)
{
    double llh[3],neu[3];
    double pos[3],vel[3],clk[3],los[3];
    double tmat[3][3];

    if (eph.vflag != 1)
        return (-1); // Invalid

    xyz2llh(xyz,llh);
    ltcmat(llh, tmat);

    satpos(eph, g, pos, vel, clk);
    subVect(los, pos, xyz);
    ecef2neu(los, tmat, neu);
    neu2azel(azel, neu);

    //return 1;

    if (azel[1]*R2D > elvMask)
        return (1); // Visible
    // else
    return (0); // Invisible
}

int allocate_channel(beidou_channel *chan, const ephemeris *eph, ionoutc_t ionoutc, beidou_time grx, double *xyz, double elvMask)
{

    int nsat=0;
    int i,sv;
    double azel[2];

    range_t rho;
    double ref[3]={0.0};
    double r_ref,r_xyz;
    double phase_ini;

    for (sv=0; sv<MAX_SAT; sv++)
    {
        if(checkSatVisibility(eph[sv], grx, xyz, 0.0, azel)==1)
        {
            nsat++; // Number of visible satellites

            if (allocatedSat[sv]==-1) // Visible but not allocated
            {
                // Allocated new satellite
                for (i=0; i<MAX_CHAN_SIM; i++)
                {
                    if (chan[i].prn_num==0)
                    {
                        // Initialize channel
                        chan[i].prn_num = sv+1;
                        chan[i].azel[0] = azel[0];
                        chan[i].azel[1] = azel[1];

                        // C/A code generation
                        prn_code_gen(chan[i].prn_code, chan[i].prn_num);
                        //TODO Generate subframe

                        eph2sbf_D2(eph[sv], ionoutc, chan[i].subframe);

                        //TODO Generate navigation message
                        nav_msg_gen(grx, &chan[i], 1);
                        // Initialize pseudorange
                        computeRange(&rho, eph[sv], &ionoutc, grx, xyz);
                        chan[i].rho0 = rho;

                        // Initialize carrier phase
                        r_xyz = rho.range;

                        computeRange(&rho, eph[sv], &ionoutc, grx, ref);
                        r_ref = rho.range;

                        phase_ini = (2.0*r_ref - r_xyz)/LAMBDA_B1;
                        phase_ini -= floor(phase_ini);
                        chan[i].carr_phase = (unsigned int)(512 * 65536.0 * phase_ini);
                        // Done.
                        break;
                    }
                }

                // Set satellite allocation channel
                if (i<MAX_CHAN_SIM)
                    allocatedSat[sv] = i;
            }
        }
        else if (allocatedSat[sv]>=0) // Not visible but allocated
        {
            // Clear channel
            chan[allocatedSat[sv]].prn_num = 0;

            // Clear satellite allocation flag
            allocatedSat[sv] = -1;
        }
    }
    return(nsat);
}

int nav_msg_gen(beidou_time bd_time, beidou_channel *chan, int init){
    beidou_time g0;
    unsigned long week,sow;

    g0.week = bd_time.week;
    // TODO +0.05
    g0.second = (double)(((unsigned long)(bd_time.second+0.05))/3UL) * 3.0; // Align with the full frame length = 3 sec
    chan->g0 = g0; // Data bit reference time

    for(int i = 0; i < 5; ++i){
        for(int j = 0; j < 10; ++j){
            if(j % 5 != 0)
                nav_word_gen(chan->subframe[i][j], false, chan->subframe_word_bits[(i*10+j)]);
            else
                nav_word_gen(chan->subframe[i][j], true, chan->subframe_word_bits[(i*10+j)]);
        }
    }
    /*
    for(int i = 1; i < 6; ++i){
        for(int k = 0; k < 10; ++k){
            for(int j = 0; j < 30; ++j)
                chan->subframe_word_bits[i*10+k][j] = chan->subframe_word_bits[k][j];
        }
    }
    */

    for(int i = 0; i < 10; ++i)
        for(int j = 0; j < 30; ++j)
            chan->subframe_word_bits[50+i][j] = chan->subframe_word_bits[i][j];

    return 0;
}

void init(beidou_channel *chan, int prn_number){

    // TODO 暂时为静态

    chan->f_carr = 108.158;
    chan->f_code = 2046000.14;
    chan->code_phase = 10.34;

    chan->iword = 9;
    chan->ibit = 23;
    chan->icode = 1;

    chan->prn_code_bit = chan->prn_code[(int)chan->code_phase] * 2 - 1;
    chan->data_bit = chan->subframe_word_bits[chan->iword][chan->ibit] * 2  - 1;

    return;
}

void *beidou_task(void *arg)
{
    sim_t *s = (sim_t *)arg;

    int sv;

    ephemeris eph[EPHEM_ARRAY_SIZE][MAX_SAT];
    int neph,ieph;
    beidou_time g0;

    double llh[3];

    beidou_time bdt;
    // 模拟4个信道
    beidou_channel chan[MAX_CHAN_SIM];
    double elvmask = 0.0; // in degree

    int ip, qp;
    int iTable;
    short *iq_buff = NULL;
    int iq_buff_size;
    beidou_time grx;
    double delt;
    int isamp;

    int iumd;
    int numd;
    double **xyz;

    double ant_pat[37];

    UTC_time t0,tmin,tmax;
    beidou_time gmin,gmax;
    double dt;
    int igrx;

    int iduration;

    // TODO ionoutc目前不用
    ionoutc_t ionoutc;

    double tmat[3][3];

    int i,j;

    int gain[MAX_CHAN_SIM];
    double timer = 0;
    char navfile[MAX_CHAR] = "2021.1.3.txt";
    iTable = 0;

    iduration = 2000;

    //经纬高
    llh[0] = 40;
    llh[1] = 30;
    llh[2] = 100;

    iq_buff_size = NUM_IQ_SAMPLES;
    iq_buff = calloc(2*iq_buff_size, 2);
    if (iq_buff==NULL)
    {
        printf("ERROR: 分配 I/Q buff失败.\n");
        goto exit;
    }

    delt = 1.0/(double)TX_SAMPLERATE;

    // 电离层误差
    ionoutc.enable = true;

    ////////////////////////////////////////////////////////////
    // Receiver position
    ////////////////////////////////////////////////////////////
    // Allocate user motion array
    xyz = (double **)malloc(USER_MOTION_SIZE * sizeof(double**));

    if (xyz==NULL)
    {
        printf("ERROR: Faild to allocate user motion array.\n");
        goto exit;
    }

    for (i=0; i<USER_MOTION_SIZE; i++)
    {
        xyz[i] = (double *)malloc(3 * sizeof(double));

        if (xyz[i]==NULL)
        {
            for (j=i-1; j>=0; j--)
                free(xyz[i]);

            printf("ERROR: Faild to allocate user motion array.\n");
            goto exit;
        }
    }

    // Static geodetic coordinates input mode: "-l"
    printf("Using static location mode.\n");
    llh2xyz(llh,xyz[0]); // Convert llh to xyz

    numd = iduration;

    for (iumd=1; iumd<numd; iumd++)
    {
        xyz[iumd][0] = xyz[0][0];
        xyz[iumd][1] = xyz[0][1];
        xyz[iumd][2] = xyz[0][2];
    }

    // Initialize the local tangential matrix for interactive mode
    ltcmat(llh, tmat);

    printf("xyz = %11.1f, %11.1f, %11.1f\n", xyz[0][0], xyz[0][1], xyz[0][2]);
    printf("llh = %11.6f, %11.6f, %11.1f\n", llh[0]*R2D, llh[1]*R2D, llh[2]);

    neph = readRinexNavAll(eph, &ionoutc, navfile);

    if (neph==0)
    {
        printf("ERROR: No ephemeris available.\n");
        goto exit;
    }

    if (ionoutc.vflag==true)
    {
        printf("  %12.3e %12.3e %12.3e %12.3e\n",
               ionoutc.alpha0, ionoutc.alpha1, ionoutc.alpha2, ionoutc.alpha3);
        printf("  %12.3e %12.3e %12.3e %12.3e\n",
               ionoutc.beta0, ionoutc.beta1, ionoutc.beta2, ionoutc.beta3);
        printf("   %19.11e %19.11e  %9d %9d\n",
               ionoutc.A0, ionoutc.A1, ionoutc.tot, ionoutc.wnt);
        printf("%6d\n", ionoutc.dtls);
    }

    gmin.second = 0.0;
    gmax.second = 0.0;

    for (sv=0; sv<MAX_SAT; sv++)
    {
        if (eph[0][sv].vflag==1)
        {
            gmin = eph[0][sv].toc;
            tmin = eph[0][sv].t;
            break;
        }
    }

    for (sv=0; sv<MAX_SAT; sv++)
    {
        if (eph[neph-1][sv].vflag == 1)
        {
            gmax = eph[neph-1][sv].toc;
            tmax = eph[neph-1][sv].t;
            break;
        }
    }

    g0 = gmin;
    t0 = tmin;

    printf("tmin = %4d/%02d/%02d,%02d:%02d:%02.0f (%d:%.0f)\n",
           tmin.year, tmin.month, tmin.day, tmin.hour, tmin.minute, tmin.second,
           gmin.week, gmin.second);
    printf("tmax = %4d/%02d/%02d,%02d:%02d:%02.0f (%d:%.0f)\n",
           tmax.year, tmax.month, tmax.day, tmax.hour, tmax.minute, tmax.second,
           gmax.week, gmax.second);

    printf("Start time = %4d/%02d/%02d,%02d:%02d:%02.0f (%d:%.0f)\n",
           t0.year, t0.month, t0.day, t0.hour, t0.minute, t0.second, g0.week, g0.second);
    printf("Duration = %.1f [sec]\n", ((double)numd)/10.0);

    // Select the current set of ephemerides
    ieph = -1;
    for (i=0; i<neph; i++)
    {
        for (sv=0; sv<MAX_SAT; sv++)
        {
            if (eph[i][sv].vflag == 1)
            {
                dt = subBeidouTime(g0, eph[i][sv].toc);
                if (dt>=-SECONDS_IN_HOUR && dt<SECONDS_IN_HOUR)
                {
                    ieph = i;
                    break;
                }
            }
        }

        if (ieph>=0) // ieph has been set
            break;
    }

    if (ieph == -1)
    {
        printf("ERROR: No current set of ephemerides has been found.\n");
        goto exit;
    }

    ////////////////////////////////////////////////////////////
    // 初始化发送信道
    ////////////////////////////////////////////////////////////

    // 清除信道信息
    for(i = 0; i < MAX_CHAN_SIM; ++i)
        chan[i].prn_num = 0;

    // Clear satellite allocation flag
    for (sv=0; sv<MAX_SAT; sv++)
        allocatedSat[sv] = -1;

    // Initial reception time
    grx = inBeidouTime(g0, 0.0);

    // 分配信道
    allocate_channel(chan, eph[ieph], ionoutc, grx, xyz[0], elvmask);

    printf("模拟的卫星号:");
    for(i = 0; i < MAX_CHAN_SIM; ++i)
        if(chan[i].prn_num > 0)
            printf("%d,",chan[i].prn_num);
    putchar('\n');

    for(i=0; i<MAX_CHAN_SIM; i++)
    {
        if (chan[i].prn_num>0)
            printf("%02d %6.1f %5.1f %11.1f %5.1f\n", chan[i].prn_num,
                   chan[i].azel[0]*R2D, chan[i].azel[1]*R2D, chan[i].rho0.d, chan[i].rho0.iono_delay);
    }

    // Update receiver time
    grx = inBeidouTime(grx, 0.1);

    for(iumd=1; iumd<numd; iumd++){
        //TODO 天线增益、路径损失

        for (i=0; i<MAX_CHAN_SIM; i++)
        {
            if (chan[i].prn_num>0)
            {
                // Refresh code phase and data bit counters
                range_t rho;
                sv = chan[i].prn_num-1;

                // Current pseudorange
                computeRange(&rho, eph[ieph][sv], &ionoutc, grx, xyz[iumd]);
                chan[i].azel[0] = rho.azel[0];
                chan[i].azel[1] = rho.azel[1];

                // Update code phase and data bit counters
                computeCodePhase(&chan[i], rho, 0.1);
                chan[i].carr_phasestep = (int)(512 * 65536.0 * chan[i].f_carr * delt);
            }
        }
        // 生成发射数据
        for(isamp = 0; isamp < iq_buff_size; ++isamp){
            int i_acc = 0;
            int q_acc = 0;

            for(i = 0; i < MAX_CHAN_SIM; ++i){
                if(chan[i].prn_num > 0){
                    iTable = (chan[i].carr_phase >> 16) & 511;

                    ip =  chan[i].data_bit * chan[i].prn_code_bit * cosTable512[iTable] ;
                    qp =  chan[i].data_bit * chan[i].prn_code_bit * sinTable512[iTable] ;

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
                            //printf("prn=%d,iword=%d,ibit=%d\n",chan[i].prn_num,chan[i].iword, chan[i].ibit);
                            chan[i].icode = 0;
                            ++chan[i].ibit;
                            // 30 bits = 1 word
                            if(chan[i].ibit >= WORD_LEN){
                                chan[i].ibit = 0;
                                //printf("prn=%d, iword= %d\n",chan[i].prn_num, chan[i].iword);
                                ++chan[i].iword;
                                // TODO D2 5个word 等于1帧
                                if(chan[i].iword >= N_WORD){
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
            //i_acc = (i_acc+64)>>7;
            //q_acc = (q_acc+64)>>7;

            // 存储I/Q buff
            iq_buff[isamp*2] = (short)i_acc;
            iq_buff[isamp*2+1] = (short)q_acc;
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

        // TODO +0.05 reason
        igrx = (int)(grx.second*10.0+0.05);

        //TODO SOW+3s
        if (igrx%30==0) // Every 3 seconds
        {
            // Update navigation message
            for (i=0; i<MAX_CHAN_SIM; i++)
            {
                if (chan[i].prn_num>0)
                    nav_msg_gen(grx, &chan[i], 0);
            }

            // Refresh ephemeris and subframes
            // Quick and dirty fix. Need more elegant way.
            for (sv=0; sv<MAX_SAT; sv++)
            {
                if (eph[ieph+1][sv].vflag==1)
                {
                    dt = subBeidouTime(eph[ieph+1][sv].toc, grx);
                    if (dt<SECONDS_IN_HOUR)
                    {
                        ieph++;

                        for (i=0; i<MAX_CHAN_SIM; i++)
                        {
                            // Generate new subframes if allocated
                            if (chan[i].prn_num!=0)
                                eph2sbf_D2(eph[ieph][chan[i].prn_num-1],ionoutc,chan[i].subframe);
                        }
                    }

                    break;
                }
            }
            // Update channel allocation
            allocate_channel(chan, eph[ieph], ionoutc, grx, xyz[iumd], elvmask);
        }
        // Update receiver time
        grx = inBeidouTime(grx, 0.1);

        // Update time counter
        printf("\rTime into run = %4.1f", subBeidouTime(grx, g0));
        fflush(stdout);
    }
    s->finished = true;
    free(iq_buff);

    exit:
    return (NULL);
}