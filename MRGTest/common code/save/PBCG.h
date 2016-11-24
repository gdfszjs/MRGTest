//////////////////////////////////////////////////////////////////////
// ***************************************************************
//  PBCG          version:  v1.0          date: 09/20/2010
//  -------------------------------------------------------------
//  Name     : PBCG.H
//  Purpose  : Ԥ����˫�����ݶȷ����Ax=b
//  Author   : Chen Wenfei
//  -------------------------------------------------------------
//  Copyright (C) 2010 - All Rights Reserved.
// ***************************************************************
//  ���÷���:
//     const int       ITOL = 1, ITMAX = 75;
//     const double    TOL  = 1.0e-9;
//     double          err  = 0.0;
//     int             iter = 0;
//     PBCG pbcgSolver(sa_d, ija_d);
//     pbcgSolver.linbcg(x, b, NP, iter, err, ITOL, TOL, ITMAX);
// ***************************************************************
#pragma once
#ifndef _INC_PBCG_H
#define _INC_PBCG_H

//////////////////////////////////////////////////////////////////////////
//macro define
//------------------------------------------------------------------------
//������ʽ(1,2,3,4)
#ifndef DEF_ITOL_DEFAULT
#define DEF_ITOL_DEFAULT 1
#endif
//------------------------------------------------------------------------
//��������
#ifndef DEF_TOL_DEFAULT
#define DEF_TOL_DEFAULT 1.0e-9
#endif
//------------------------------------------------------------------------
//����������
#ifndef DEF_ITMAX_DEFAULT
#define DEF_ITMAX_DEFAULT 75
#endif

//////////////////////////////////////////////////////////////////////////
// MACRO that forbids to use the copy ctor & assign statement
//------------------------------------------------------------------------
#ifndef DISALLOW_COPY_AND_ASSIGN
#define DISALLOW_COPY_AND_ASSIGN(ClassName)   \
    ClassName(const ClassName&);              \
    ClassName& operator=(const ClassName&)
#endif

//////////////////////////////////////////////////////////////////////////
class PBCG
{
public://ctor & dtor
    PBCG(const double *sa, const unsigned long *ija);
    ~PBCG(void);

public://interface functions
    void linbcg(double x[],                     //������(out)
        const double b[],                       //������ʽ����(in)
        const unsigned long n,                  //����ά��
        int &iter,                              //ʵ�ʵ�������
        double &err,                            //ʵ�ʵ�������
        int itol = DEF_ITOL_DEFAULT,            //������ʽ(1,2,3,4)
        double tol = DEF_TOL_DEFAULT,           //�����������
        int itmax = DEF_ITMAX_DEFAULT) const;   //�����������

public:
    void sprstx(double b[], const double x[], const unsigned long n) const;
    void sprsax(double b[], const double x[], const unsigned long n) const;

private:
    void atimes(double b[], const double x[], const unsigned long n, const int itrnsp) const;
    void asolve(double x[], const double b[], const unsigned long n, const int itrnsp) const;
    double snrm(const double sx[], const unsigned long n, const int itol = 1) const;

private://������ϡ�����
    const double          *_sa;
    const unsigned long   *_ija;

private://forbid to use the copy ctor & assign statement
    DISALLOW_COPY_AND_ASSIGN(PBCG);
};

//////////////////////////////////////////////////////////////////////////
//global functions
//------------------------------------------------------------------------
void sprsin(double *const sa,
            unsigned long *const ija,       //������ϡ�����(out)
            const unsigned long nmax,       //��������һ��Ϊ2*n*n+1
            const double *const a,          //ϡ��ϵ������A(in)
            const unsigned long n,          //�����ά��n
            const double thresh = 0.0);     //��ֵ
//------------------------------------------------------------------------
void sprsout(double *const a,               //ϡ��ϵ������A(out)
             const unsigned long n,         //����ά��n
             const double *const sa,
             const unsigned long *const ija);//������ϡ�����(in)
//------------------------------------------------------------------------
void CreateLapMat(double *const sa,
                  unsigned long *const ija, //������ϡ�����(out)
                  const int iWidth,         //Դͼ��Ŀ��
                  const int iHeight);       //Դͼ��ĸ߶�

#endif // !_INC_PBCG_H
