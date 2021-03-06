// ***********************************************************************
//  PBCG                version:  v1.0                date: 09/25/2010
//  ----------------------------------------------------------------------
//  FileName : PBCG.CPP
//  FilePath : e:\VPPlat\examples\PBCG
//  Purpose  : 预条件双共轭梯度法求解Ax=b
//  Author   : Chen Wenfei
//  Modified by
//  Created  : 2010/09/25 16:38
//  Copyright: (c) Chen Wenfei
//  ----------------------------------------------------------------------
//  Copyright (C) 2010 - All Rights Reserved.
// ***********************************************************************
// Note:
// 怎样能构成"好"的收敛依赖于应用.
// 因此程序PBCG::linbcg提供了4种选择, 可以对输入标记itol进行设置.
// ①如果设itol=1, 则当值|A·x-b|/|b|小于输入量tol时, 迭代停止.
// ②如果设itol=2, 则所需的准则为|Ã-1·(A·x-b)|/|Ã-1·b|<tol.
// ③如果设itol=3, 则程序使用其自身对x误差的估计, 并要求其大小除以x的值<tol.
// ④设置itol=4与itol=3情况相同, 只是用误差的最大量(绝对值)和x的最大分量
// 代替了向量的大小(即用L∞范数代替了L²范数)
// ★需要在实践中试着寻找到最适合自己的收敛准则.
// ***********************************************************************

//#include "stdafx.h"
#include "PBCG.h"
#include <iostream>
#include <cmath>
#include <cassert>

//////////////////////////////////////////////////////////////////////////
//internal functions
//------------------------------------------------------------------------
static void nrerror(const char error_text[])
/* Numerical Recipes standard error handler */
{
    fprintf(stderr, "Numerical Recipes run-time error...\n");
    fprintf(stderr, "%s\n", error_text);
    fprintf(stderr, "...now exiting to system...\n");
    throw "Please processing...\n";
    //exit(1);
}

//////////////////////////////////////////////////////////////////////////
// ctor & dtor of class PBCG
//------------------------------------------------------------------------
PBCG::PBCG(const double *sa, const unsigned long *ija)
: _sa(sa), _ija(ija)
{}

PBCG::~PBCG(void)
{
    _sa = 0, _ija = 0;
}
//////////////////////////////////////////////////////////////////////////
//interface functions
void PBCG::linbcg(double x[],             //解向量(out)
                  const double b[],       //方程右式向量(in)
                  const unsigned long n,  //方阵维数
                  int &iter,              //实际迭代次数
                  double &err,            //实际迭代精度
                  int itol,               //收敛方式(1,2,3,4)
                  double tol,             //所需迭代精度
                  int itmax) const        //所需迭代次数
{
    assert(x != 0 && "x is null in PBCG::linbcg()!");
    assert(b != 0 && "b is null in PBCG::linbcg()!");

    double ak, bk, bnrm, dxnrm, xnrm, zm1nrm, znrm;
    const double EPS = 1.0e-14;

    double *p,*pp,*r,*rr,*z,*zz;
    p  = (double*)calloc(n, sizeof(double));
    if (0 == p) nrerror("linbcg: allocate memory failed for p!");
    pp = (double*)calloc(n, sizeof(double));
    if (0 == pp) nrerror("linbcg: allocate memory failed for pp!");
    r  = (double*)calloc(n, sizeof(double));
    if (0 == r) nrerror("linbcg: allocate memory failed for r!");
    rr = (double*)calloc(n, sizeof(double));
    if (0 == rr) nrerror("linbcg: allocate memory failed for rr!");
    z  = (double*)calloc(n, sizeof(double));
    if (0 == z) nrerror("linbcg: allocate memory failed for z!");
    zz = (double*)calloc(n, sizeof(double));
    if (0 == zz) nrerror("linbcg: allocate memory failed for zz!");

    // 计算初始余量
    iter = 0;

    //atimes的输入为x[0..n-1], 输出为r[0..n-1], 最后的0表示使用矩阵本身(而不是其转置)
    atimes(r, x, n, 0);
    for (unsigned long i = 0; i < n; i++)
    {
        r[i]  = b[i] - r[i];
        rr[i] = r[i];
    }

    //去掉注释符号即可得到该算法的"最小余量"变形
    // atimes(_sa, _ija, r, rr, n, 0);

    if (1 == itol)
    {
        bnrm = snrm(b, n, itol);
        //asolve的输入是r[0..n-1], 输出是z[0..n-1],
        //最后的0表示使用矩阵Ã(而不是其转置)
        asolve(z, r, n, 0);
    }
    else if (2 == itol)
    {
        asolve(z, b, n, 0);
        bnrm=snrm(z, n, itol);
        asolve(z, r, n, 0);
    }
    else if (3 == itol || 4 == itol)
    {
        asolve(z, b, n, 0);
        bnrm=snrm(z,itol,n);
        asolve(z, r, n, 0);
        znrm=snrm(z, n, itol);
    }
    else
    {
        nrerror("linbcg: illegal itol!");
    }

    double bkden = 1.0;
    while (iter < itmax)
    {   //主循环
        ++iter;
        asolve(zz, rr, n, 1);   //最后的1表示使用转置矩阵Ã'

        double bknum = 0.0;
        for (unsigned long j = 0; j < n; j++) bknum += z[j] * rr[j];

        //calculate coefficient bk and direction vectors p and pp
        if (1 == iter)
        {
            for (unsigned long j = 0; j < n; j++)
            {
                p[j]  = z[j];
                pp[j] = zz[j];
            }
        }
        else
        {
            bk=bknum/bkden;
            for (unsigned long j = 0; j < n; j++)
            {
                p[j]  = bk * p[j]  + z[j];
                pp[j] = bk * pp[j] + zz[j];
            }
        }

        //计算系数ak, 新的迭代量x, 以及新的余量r和rr
        bkden = bknum;
        atimes(z, p, n, 0);

        double akden = 0.0;
        for (unsigned long j = 0; j < n; j++)   akden += z[j]*pp[j];
        ak = bknum / akden;
        atimes(zz, pp, n, 1);
        for (unsigned long j = 0; j < n; j++)
        {
            x[j]  += ak * p[j];
            r[j]  -= ak * z[j];
            rr[j] -= ak * zz[j];
        }

        //求解Ã·z=r, 并检查停止准则
        asolve(z, r, n, 0);
        if (1 == itol)
            err=snrm(r, n, itol) / bnrm;
        else if (2 == itol)
            err=snrm(z, n, itol)/bnrm;
        else if (3 == itol || 4 == itol)
        {
            zm1nrm = znrm;
            znrm = snrm(z, n, itol);

//             printf("iter=%4d err=%12.6f\n",  iter,  err);

            if (fabs(zm1nrm-znrm) > EPS * znrm)
            {
                dxnrm=fabs(ak) * snrm(p, n, itol);
                err = znrm / fabs(zm1nrm - znrm) * dxnrm;
            }
            else
            {   //误差可能不精确, 故再次循环
                err = znrm / bnrm;
                continue;
            }

            xnrm = snrm(x, n, itol);

//             printf("iter=%4d err=%12.6f\n",  iter,  err);

            if (err <= 0.5 * xnrm)
                err /= xnrm;
            else
            {   //误差可能不精确, 故再次循环
                err = znrm / bnrm;
                continue;
            }
        }

//         printf("iter=%4d err=%12.6f\n",  iter,  err);

        if (err <= tol) break;
    }

    //释放资源
    free(p),  p = 0;
    free(pp), pp = 0;
    free(r),  r = 0;
    free(rr), rr = 0;
    free(z),  z = 0;
    free(zz), zz = 0;
}

//////////////////////////////////////////////////////////////////////////
//often used computation of matrix
//------------------------------------------------------------------------
void PBCG::sprsax(double b[],                   //方程右式向量(out)
                  const double x[],             //解向量(in)
                  const unsigned long n) const  //向量维度n
{   //将行索引稀疏矩阵存储为数组_sa和_ija乘以一个向量x[0..n-1], 得到向量b[0..n-1]
    assert(x != 0 && "x is null in PBCG::sprsax()!");
    assert(b != 0 && "b is null in PBCG::sprsax()!");

    if (_ija[0] != n + 1)
        nrerror("sprsax: mismatched vector and matrix");
    for (unsigned long i = 0; i < n; i++)
    {
        b[i] = _sa[i] * x[i];
        for (unsigned long k = _ija[i]; k < _ija[i + 1]; k++)
            b[i] += _sa[k] * x[_ija[k]];
    }
}
//------------------------------------------------------------------------
void PBCG::sprstx(double b[],                   //方程右式向量(out)
                  const double x[],             //解向量(in)
                  const unsigned long n) const  //向量维度n
{   //将存储成数组_sa和_ija的行索引稀疏矩阵的转置右乘向量x[0..n-1], 得到向量b[0..n-1]
    assert(x != 0 && "x is null in PBCG::sprstx()!");
    assert(b != 0 && "b is null in PBCG::sprstx()!");

    if (_ija[0] != n + 1)
        nrerror("sprstx: mismatched vector and matrix");

    for (unsigned long i = 0; i < n; ++i)   //从对角项开始
        b[i] = _sa[i] *x[i];

    for (unsigned long i = 0; i < n; i++)  //对非对角项循环
    {
        for (unsigned long k = _ija[i]; k < _ija[i + 1]; k++)
        {
            unsigned long j = _ija[k];
            b[j] += _sa[k] * x[i];
        }
    }
}

//////////////////////////////////////////////////////////////////////////
void PBCG::asolve(double x[],               //解向量(out)
                  const double b[],         //方程右式向量(in)
                  const unsigned long n,    //向量长度
                  const int itrnsp) const   //标记(是否转置)
{   //预条件矩阵Ã是A的对角线部分, 存于_sa的前n个元素.
    //因为该矩阵的转置对角元不变, 故没有用标记itrnsp.
    assert(x != 0 && "x is null in PBCG::asolve()!");
    assert(b != 0 && "b is null in PBCG::asolve()!");

    for (unsigned long i = 0; i < n; i++)
        x[i] = (_sa[i] != 0.0 ? b[i] / _sa[i]: b[i]);
}
//------------------------------------------------------------------------
void PBCG::atimes(double b[],               //方程右式向量(out)
                  const double x[],         //解向量(in)
                  const unsigned long n,    //向量长度
                  const int itrnsp) const   //标记(0:Ax=b; 1:A'x=b)
{   //计算A或其转置A'与一个向量的乘积
    assert(x != 0 && "x is null in PBCG::atimes()!");
    assert(b != 0 && "b is null in PBCG::atimes()!");

    if (itrnsp) sprstx(b, x, n);
    else        sprsax(b, x, n);
}
//------------------------------------------------------------------------
double PBCG::snrm(const double sx[], const unsigned long n, const int itol) const
{   //根据itol的值, 计算向量sx[0..n-1]的两个范数之一
    assert(sx != 0 && "sx is null in PBCG::snrm()!");

    if (itol <= 3)
    {   //向量大小的范数
        double ans = 0.0;
        for (unsigned long i = 0; i < n; i++)
            ans += sx[i] * sx[i];
        return sqrt(ans);
    }
    else
    {
        unsigned long isamax = 0;
        for (unsigned long i = 0; i < n; i++)
        {   //最大元的范数
            if (fabs(sx[i]) > fabs(sx[isamax]))
                isamax = i;
        }
        return fabs(sx[isamax]);
    }
}

//////////////////////////////////////////////////////////////////////////
//global functions
//------------------------------------------------------------------------
void sprsin(double *const sa, unsigned long *const ija, //行索引稀疏矩阵(out)
            const unsigned long nmax,                   //向量长度一般为2*n*n+1
            const double *const a,                      //稀疏系数方阵A(in)
            const unsigned long n,                      //方阵的维数n
            const double thresh)                        //阈值
{   //将方阵a[0..n-1][0..n-1]转换成行索引稀疏存储模式.
    //仅保留数值大于阈值的元素, 输出2个维数nmax(作为输入参数)的线性数组:
    //sa[0..]包含数组值, ija[0..]为索引号.输出sa和ija中的元素数目都是ija[ija[0]-1]
    assert(sa != 0 && "sa is null in sprsin()!");
    assert(ija != 0 && "ija is null in sprsin()!");
    assert(a != 0 && "a is null in sprsin()!");

    for (unsigned long i = 0; i < n; i++)   //存储对角线元
        sa[i] = a[i * n + i];
    ija[0] = n + 1;                         //根据规则, ija[0] = n + 1

    unsigned long k = n;
    for (unsigned long i = 0; i < n; i++)
    {   //行循环
        for (unsigned long j = 0; j < n; j++)
        {   //列循环
            if (fabs(a[i * n + j]) > thresh && i != j)
            {
                if (++k > nmax)
                    nrerror("sprsin: sa and ija too small");
                sa[k] = a[i * n + j];  //存储非对角元及其列号
                ija[k] = j;
            }
        }
        ija[i + 1] = k + 1;   //每行完成后,将指针指向下一个
    }
}
//------------------------------------------------------------------------
void sprsout(double *const a,               //稀疏系数方阵A(out)
             const unsigned long n,         //方阵维数n
             const double *const sa,
             const unsigned long *const ija)//行索引稀疏矩阵(in)
{   //将行索引稀疏存储模式转换成方阵a[0..n-1][0..n-1].
    //sa[0..]包含数组值, ija[0..]为索引号.输出sa和ija中的元素数目都是ija[ija[0]-1]
    assert(sa != 0 && "sa is null in sprsout()!");
    assert(ija != 0 && "ija is null in sprsout()!");
    assert(a != 0 && "a is null in sprsout()!");

    for (unsigned long i = 0; i < n; i++)
    {
        a[i * n + i]=sa[i];
        for (unsigned long j = ija[i]; j < ija[i + 1]; j++)
        {
            a[i * n + ija[j]]=sa[j];
        }
    }
}
//------------------------------------------------------------------------
void CreateLapMat(double *const sa,
                  unsigned long *const ija, //行索引稀疏矩阵(out)
                  const int iWidth,         //源图像的宽度
                  const int iHeight)        //源图像的高度
{   //生成laplace矩阵
    assert(sa != 0 && "sa is null in CreateLapMat()!");
    assert(ija != 0 && "ija is null in CreateLapMat()!");

    unsigned long n = iWidth * iHeight;     //矩阵维数

    for (unsigned long i = 0; i < n; i++)   //存储对角线元
        sa[i] = -4.0;
    ija[0] = n + 1;                         //根据规则, ija[0] = n + 1

    unsigned long k = n;
    for (int i = 0; i < iHeight; i++)
    {   //对图像的行循环
        for (int j = 0; j < iWidth; j++)
        {   //对图像的列循环
            //获取4邻域的坐标索引
            int iRowIndex = i * iWidth + j; //等价于稀疏矩阵的行号

            if ((i - 1) >= 0)
            {   //Upper N(p)
                ++k;
                sa[k] = 1.0;  //存储非对角元及其列值
                ija[k] = (i - 1) * iWidth + j;
            }

            if ((j - 1) >= 0)
            {   //Left N(p)
                ++k;
                sa[k] = 1.0;  //存储非对角元及其列值
                ija[k] = i * iWidth + (j - 1);
            }

            if ((j + 1) < iWidth)
            {   //Right N(p)
                ++k;
                sa[k] = 1.0;  //存储非对角元及其列值
                ija[k] = i * iWidth + (j + 1);
            }

            if ((i + 1) < iHeight)
            {   //Down N(p)
                ++k;
                sa[k] = 1.0;  //存储非对角元及其列值
                ija[k] = (i + 1) * iWidth + j;
            }

            ija[iRowIndex + 1] = k + 1; //每行完成后,将指针指向下一个
        }
    }

    return;
}