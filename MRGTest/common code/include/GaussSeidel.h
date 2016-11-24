//2008/2/22
//Nie yong wei 
//Modified on 2008/3/22 

/***************************
@�ø�˹-���¶����������Է����顣ע�⣬ֻ���������ӵ��С��ϵ����������Է����顣
@���룺ϵ�������õ�NIE_DATASTRUCTURE�����Matrix2D���洢��
	   �Ҷ˾���ͬ�ϣ�
@������������ͬ�ϣ�
****************************/


#ifndef _GAUSS_SEIDEL_H_
#define _GAUSS_SEIDEL_H_

#include "datastructure.h"
#include <cmath>
using namespace NIE_DATASTRUCTURE;

typedef Matrix2D<double> MAT_2D;
typedef const Matrix2D<double> MAT_C_2D;


namespace NIE_GAUSSSEIDEL{

class CGaussSeidel
{
private:
	double e;
public:
	CGaussSeidel(double _e);
	void gs(MAT_2D & res,MAT_C_2D & coemat,MAT_C_2D & rmat);
};

}

NIE_GAUSSSEIDEL::CGaussSeidel::CGaussSeidel(double _e) : e(_e) { }

void NIE_GAUSSSEIDEL::CGaussSeidel::gs(MAT_2D & res,MAT_C_2D & coemat,MAT_C_2D & rmat)
{
	int n = coemat.getnrows();
	double _e,t;
	while(1)
	{
		_e = 0.0;
		for(int i = 0;i<n;i++)
		{
			t = res[i][0];
			double tmp = 0.0;
			for(int j = 0;j<n;j++)
				if(j != i)
					tmp += coemat[i][j] * res[j][0];
			res[i][0] = (rmat[i][0] - tmp) / coemat[i][i];

			if(fabs(t - res[i][0]) > _e)
				_e = fabs(t-res[i][0]);
		}
		if(_e<e)
			break;
	}
}


 
#endif   /* _GAUSS_SEIDEL_H_ */