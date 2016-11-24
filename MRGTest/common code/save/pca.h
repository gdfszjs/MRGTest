/*
	   ��ɳɷַ���(Principal Components Analysis)
	   2008/4/8 by Yongwei Nie
	   
*/
#pragma 
#include "datastructure.h"
#include <vector>
#include <cmath>

typedef NIE_DATASTRUCTURE::Matrix2D<double> MAT_D_2D;
typedef const NIE_DATASTRUCTURE::Matrix2D<double> MAT_C_D_2D;

namespace NIE_PCA
{

template<class T>
inline const T SQR(const T a) {return a*a;}

inline float SIGN(const double &a, const float &b)
{return b >= 0 ? (a >= 0 ? a : -a) : (a >= 0 ? -a : a);}



class PCA
{
public:
	/*
		������ֵ����������
	*/
	
	//��ת
	void rot(MAT_D_2D &a, const double s, const double tau, const int i,const int j, const int k, const int l);
	//����ʵ�Գƾ���a����������ֵ������������d����a������ֵ��v����������a�Ĺ淶��������
	//nrot����������ſ˱���ת����
	void jacobi(MAT_D_2D &a, std::vector<double> &d, MAT_D_2D &v, int &nrot);
	//��jacobi���������ֵ���������У�������������v������Ӧ�ĵ�����ֱ�Ӳ��뷨��
	void eigsrt(std::vector<double> &d, MAT_D_2D &v);
	double pythag(const double a, const double b);
	void tred2(MAT_D_2D &a, std::vector<double> &d, std::vector<double> &e);
	void tqli(std::vector<double> &d, std::vector<double> &e, MAT_D_2D &z);
	
	/*
		Э�������
	*/
	
	//�������s����ÿ���ϵľ�ֵm���s��Э�������a
	void covariancemat(MAT_D_2D & a,MAT_C_D_2D & s,const std::vector<double> & m);

	/* 
		����һ�����
	*/

	//����ת��
	void transpose(MAT_D_2D & u , MAT_C_D_2D & c);
	//����˷�
	void multi(MAT_D_2D & a , MAT_C_D_2D & u ,MAT_C_D_2D & c);
	//���󿽱�
	void copy(MAT_D_2D & u,MAT_C_D_2D & c);

	/*
		PCA
    */

	//����ԭ����s��ʣ��ά��dim������jacobi��ת�����ͽ������
	MAT_D_2D* pca(MAT_D_2D & s,int dim ,int & nrot); 
	
};
}


/*
  
  #############################################
  #############################################
  #############################################
  ######                                 ######
  ######   I M P L E M E N T A T I O N   ######
  ######                                 ######
  #############################################
  #############################################
  #############################################
  
*/


void NIE_PCA::PCA::rot(MAT_D_2D &a, const double s, const double tau, const int i,const int j, const int k, const int l)
{
	double g,h;

	g=a[i][j];
	h=a[k][l];
	a[i][j]=g-s*(h+g*tau);
	a[k][l]=h+s*(g-h*tau);
}


void NIE_PCA::PCA::jacobi(MAT_D_2D &a, std::vector<double> &d, MAT_D_2D &v, int &nrot)
{
	int i,j,ip,iq;
	double tresh,theta,tau,t,sm,s,h,g,c;

	int n=(int)d.size();
	std::vector<double> b(n),z(n);
	for (ip=0;ip<n;ip++) {
		for (iq=0;iq<n;iq++) v[ip][iq]=0.0;
		v[ip][ip]=1.0;
	}
	for (ip=0;ip<n;ip++) {
		b[ip]=d[ip]=a[ip][ip];
		z[ip]=0.0;
	}
	nrot=0;
	for (i=1;i<=50;i++) {
		sm=0.0;
		for (ip=0;ip<n-1;ip++) {
			for (iq=ip+1;iq<n;iq++)
				sm += fabs(a[ip][iq]);
		}
		if (sm == 0.0)
			return;
		if (i < 4)
			tresh=0.2*sm/(n*n);
		else
			tresh=0.0;
		for (ip=0;ip<n-1;ip++) {
			for (iq=ip+1;iq<n;iq++) {
				g=100.0*fabs(a[ip][iq]);
				if (i > 4 && (fabs(d[ip])+g) == fabs(d[ip])
					&& (fabs(d[iq])+g) == fabs(d[iq]))
						a[ip][iq]=0.0;
				else if (fabs(a[ip][iq]) > tresh) {
					h=d[iq]-d[ip];
					if ((fabs(h)+g) == fabs(h))
						t=(a[ip][iq])/h;
					else {
						theta=0.5*h/(a[ip][iq]);
						t=1.0/(fabs(theta)+sqrt(1.0+theta*theta));
						if (theta < 0.0) t = -t;
					}
					c=1.0/sqrt(1+t*t);
					s=t*c;
					tau=s/(1.0+c);
					h=t*a[ip][iq];
					z[ip] -= h;
					z[iq] += h;
					d[ip] -= h;
					d[iq] += h;
					a[ip][iq]=0.0;
					for (j=0;j<ip;j++)
						rot(a,s,tau,j,ip,j,iq);
					for (j=ip+1;j<iq;j++)
						rot(a,s,tau,ip,j,j,iq);
					for (j=iq+1;j<n;j++)
						rot(a,s,tau,ip,j,iq,j);
					for (j=0;j<n;j++)
						rot(v,s,tau,j,ip,j,iq);
					++nrot;
				}
			}
		}
		for (ip=0;ip<n;ip++) {
			b[ip] += z[ip];
			d[ip]=b[ip];
			z[ip]=0.0;
		}
	}
}


void NIE_PCA::PCA::eigsrt(std::vector<double> &d, MAT_D_2D &v)
{
	int i,j,k;
	double p;

	int n=(int)d.size();
	for (i=0;i<n-1;i++) {
		p=d[k=i];
		for (j=i;j<n;j++)
			if (d[j] >= p) p=d[k=j];
		if (k != i) {
			d[k]=d[i];
			d[i]=p;
			for (j=0;j<n;j++) {
				p=v[j][i];
				v[j][i]=v[j][k];
				v[j][k]=p;
			}
		}
	}
}

double NIE_PCA::PCA:: pythag(const double a, const double b)
{
	double absa,absb;

	absa=fabs(a);
	absb=fabs(b);
	if (absa > absb) return absa*sqrt(1.0+SQR(absb/absa));
	else return (absb == 0.0 ? 0.0 : absb*sqrt(1.0+SQR(absa/absb)));
}

void NIE_PCA::PCA::tred2(MAT_D_2D &a, std::vector<double> &d, std::vector<double> &e)
{
	int l,k,j,i;
	double scale,hh,h,g,f;

	int n=(int)d.size();
	for (i=n-1;i>0;i--) {
		l=i-1;
		h=scale=0.0;
		if (l > 0) {
			for (k=0;k<l+1;k++)
				scale += fabs(a[i][k]);
			if (scale == 0.0)
				e[i]=a[i][l];
			else {
				for (k=0;k<l+1;k++) {
					a[i][k] /= scale;
					h += a[i][k]*a[i][k];
				}
				f=a[i][l];
				g=(f >= 0.0 ? -sqrt(h) : sqrt(h));
				e[i]=scale*g;
				h -= f*g;
				a[i][l]=f-g;
				f=0.0;
				for (j=0;j<l+1;j++) {
				// Next statement can be omitted if eigenvectors not wanted
					a[j][i]=a[i][j]/h;
					g=0.0;
					for (k=0;k<j+1;k++)
						g += a[j][k]*a[i][k];
					for (k=j+1;k<l+1;k++)
						g += a[k][j]*a[i][k];
					e[j]=g/h;
					f += e[j]*a[i][j];
				}
				hh=f/(h+h);
				for (j=0;j<l+1;j++) {
					f=a[i][j];
					e[j]=g=e[j]-hh*f;
					for (k=0;k<j+1;k++)
						a[j][k] -= (f*e[k]+g*a[i][k]);
				}
			}
		} else
			e[i]=a[i][l];
		d[i]=h;
	}
	// Next statement can be omitted if eigenvectors not wanted
	d[0]=0.0;
	e[0]=0.0;
	// Contents of this loop can be omitted if eigenvectors not
	//	wanted except for statement d[i]=a[i][i];
	for (i=0;i<n;i++) {
		l=i;
		if (d[i] != 0.0) {
			for (j=0;j<l;j++) {
				g=0.0;
				for (k=0;k<l;k++)
					g += a[i][k]*a[k][j];
				for (k=0;k<l;k++)
					a[k][j] -= g*a[k][i];
			}
		}
		d[i]=a[i][i];
		a[i][i]=1.0;
		for (j=0;j<l;j++) a[j][i]=a[i][j]=0.0;
	}
}


void NIE_PCA::PCA::tqli(std::vector<double> &d, std::vector<double> &e, MAT_D_2D &z)
{
	int m,l,iter,i,k;
	double s,r,p,g,f,dd,c,b;

	int n=(int)d.size();
	for (i=1;i<n;i++) e[i-1]=e[i];
	e[n-1]=0.0;
	for (l=0;l<n;l++) {
		iter=0;
		do {
			for (m=l;m<n-1;m++) {
				dd=fabs(d[m])+fabs(d[m+1]);
				if (fabs(e[m])+dd == dd) break;
			}
			if (m != l) {
				g=(d[l+1]-d[l])/(2.0*e[l]);
				r=pythag(g,1.0);
				g=d[m]-d[l]+e[l]/(g+SIGN(r,g));
				s=c=1.0;
				p=0.0;
				for (i=m-1;i>=l;i--) {
					f=s*e[i];
					b=c*e[i];
					e[i+1]=(r=pythag(f,g));
					if (r == 0.0) {
						d[i+1] -= p;
						e[m]=0.0;
						break;
					}
					s=f/r;
					c=g/r;
					g=d[i+1]-p;
					r=(d[i]-g)*s+2.0*c*b;
					d[i+1]=g+(p=s*r);
					g=c*r-b;
					// Next loop can be omitted if eigenvectors not wanted
					for (k=0;k<n;k++) {
						f=z[k][i+1];
						z[k][i+1]=s*z[k][i]+c*f;
						z[k][i]=c*z[k][i]-s*f;
					}
				}
				if (r == 0.0 && i >= l) continue;
				d[l] -= p;
				e[l]=g;
				e[m]=0.0;
			}
		} while (m != l);
	}
}



void NIE_PCA::PCA::covariancemat(MAT_D_2D &a, MAT_C_D_2D & s,const std::vector<double> &m)
{
	int nn = a.getncols();
	int sn = s.getncols();
	for(int i = 0;i<nn;i++)
		for(int j = 0;j<=i;j++)
		{
			double nu = 0.0;
			for(int k = 0;k<sn;k++)
				nu+=(s[i][k] - m[i])*(s[j][k] - m[j]);
			a[i][j] = nu/(sn-1);
		}
	for(int i = 0;i<nn;i++)
		for(int j = i;j<nn;j++)
			a[i][j] = a[j][i];
}

void NIE_PCA::PCA::copy(MAT_D_2D & u,MAT_C_D_2D &c)
{
	int rr = u.getnrows();
	int cc = u.getncols();
	for(int i = 0;i<rr;i++)
		for(int j = 0;j<cc;j++)
			u[i][j] = c[i][j];
}

void NIE_PCA::PCA::transpose(MAT_D_2D & u,MAT_C_D_2D & c)
{
	int rr = u.getnrows();
	int cc = u.getncols();
	for(int i = 0;i<rr;i++)
		for(int j = 0;j<cc;j++)
			u[i][j] = c[j][i];
}

void NIE_PCA::PCA::multi(MAT_D_2D & a ,MAT_C_D_2D & u,MAT_C_D_2D & c)
{
	int rr = a.getnrows();
	int cc = a.getncols();
	int nn = u.getncols();
	for(int i = 0;i<rr;i++)
		for(int j = 0;j<cc;j++)
		{
			double tmp = 0.0;
			for(int k = 0;k<nn;k++)
				tmp += u[i][k]*c[k][j];
			a[i][j] = tmp;
		}
}

MAT_D_2D* NIE_PCA::PCA::pca(MAT_D_2D &s,int dim,int & nrot)
{
	int rr = s.getnrows();
	int cc = s.getncols();
	std::vector<double> m(rr);
	//����ÿ�еľ�ֵ
	for(int i = 0;i<rr;i++)
	{
		double nu = 0.0;
		for(int j = 0;j<cc;j++)
			nu+=s[i][j];
		m[i] = nu/cc;
	}
	//ÿ��ֵ����ȥ��Ӧ�ľ�ֵ
	for(int i = 0;i<rr;i++)
	{
		for(int j = 0;j<cc;j++)
		{
			s[i][j]-=m[i];
		}
		m[i] = 0.0;
	}
	
	MAT_D_2D t_eve(0.0,dim,rr);

	{
		//����Э�������
		MAT_D_2D comat(0.0,rr,rr);
		covariancemat(comat,s,m);

		
		//MAT_D_2D eve(0.0,rr,rr);
		//std::vector<double> eva(rr);
		//jacobi(comat,eva,eve,nrot);
		//eigsrt(eva,eve);
		////��ά
		//MAT_D_2D _eve(0.0,rr,dim);
		//for(int i = 0;i<rr;i++)
		//	for(int j = 0;j<dim;j++)
		//		_eve[i][j] = eve[i][j];

		std::vector<double> d(rr);
		std::vector<double> e(rr);
		tred2(comat,d,e);
		tqli(d,e,comat);
		eigsrt(d,comat);

		//��ά
		MAT_D_2D _eve(0.0,rr,dim);
		for(int i = 0;i<rr;i++)
			for(int j = 0;j<dim;j++)
				_eve[i][j] = comat[i][j];

		//�Խ�ά��õ��������ת��
		
		transpose(t_eve,_eve);
	}

	//����˷�
	MAT_D_2D * pr = new MAT_D_2D(0.0,dim,cc);
	multi(*pr,t_eve,s);
	return pr;
}