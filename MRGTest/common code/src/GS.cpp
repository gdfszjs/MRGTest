#include "stdafx.h"
#include "GS.h"
#include "math.h"
GS::GS(void)
{
}
GS::GS(double *&sa_,int *&ija_,double * &x_, double *&b_,const int n_)
:sa(sa_),ija(ija_),x(x_),b(b_),n(n_)
{
}
GS::~GS(void)
{
}
void GS::Do_GS(int max,double e)
{
	for(int k=0;k<max;k++)
	{
		int temp_e=0.0;
	for(int i=0;i<n;i++)
	{
		int sum=b[i];
		int temp=x[i];
		for(int j=ija[i];j<ija[i+1];j++)
			sum-=sa[j]*x[ija[j]];
		x[i]=sum/sa[i];
		if(fabs(x[i]-temp)>temp_e)
			temp_e=fabs(x[i]-temp);
	}
	if(temp_e<e)
		{
			fprintf(stderr,"�ﵽ�������� ��������Ϊ��%d :\n",k);
		return ;
		}
	}
	fprintf(stderr,"��ȫ���� ��������Ϊ��%d \n",max);
}
