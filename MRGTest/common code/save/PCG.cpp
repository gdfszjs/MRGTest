//#include "stdafx.h"
#include "PCG.h"
#include <cmath>
PCG::PCG(double * &sa_,unsigned long *&ija_,double *&x_,double *&b_,const unsigned long n_,
											const double tol_,const int itmax_)
{
	sa=sa_;
	ija=ija_;
	x=x_;
	b=b_;
	n=n_;
	sb=0;
	tol=tol_;
	itmax=itmax_;
	flag=0;
	nz=0;
}
PCG::PCG(double * &sa_,unsigned long *&ija_,double *&x_,double *&b_,const unsigned long n_,
								const double tol_,const int itmax_,const int flag_,const int nz_)
{
	sa=sa_;
	ija=ija_;
	x=x_;
	b=b_;
	n=n_;
	sb=0;
	tol=tol_;
	itmax=itmax_;
	flag=flag_;
	nz=nz_;
}
PCG::~PCG()
{
	if(sb)
	{
		delete [] sb;
		sb=0;
	}
}
bool PCG::Choldc(double *&sa1,unsigned long *& ija1,const unsigned long n)
{
	unsigned long i,j,k,t,l,m;
	double sum;
//	//////////////////////////////���Ƚ���Ԥ����ʹ����������ȶ�
//	for(i=0;i<n;i++)
//	{
//		sum=sa1[i]*sa1[i];
//		for(k=ija1[i];k<ija1[i+1];k++)
//			sum+=sa1[k]*sa1[k];
//		sa1[i]=sqrt(sum);
//	}
//
	for(i=0;i<n;i++)//��ÿһ�зֱ���д���
	{
		sum=sa1[i];//�ӶԽ��ߵ�Ԫ�ؿ�ʼ
		for(k=ija1[i];ija1[k]<i && k<ija1[i+1];k++)
			sum-=sa1[k]*sa1[k];

		if(sum>0)
			sa1[i]=sqrt(sum);
	//	else 
	//	{
	//		::AfxMessageBox("����˹���ֽ���󣡾��󲻷���Ҫ��");
	//		return false;
			
	//	}
		/////////////////�ԶԽ����Ժ��Ԫ�ش���
	for(t=ija1[i];t<ija1[i+1];t++)
		if(ija1[t]>i)
			break;
		for(j=t;j<ija1[i+1];j++)//��i�У���j��Ԫ��
		{
			sum=sa1[j];
			for(k=ija1[i];ija1[k]<i && k<ija1[i+1];k++)
				for( l=ija1[ija1[j]];ija1[l]<i && l<ija1[ija1[j]+1];l++)//�к�Ϊija[j],��Ϊ�к�
					if(ija1[l]==ija1[k])
					{
						sum-=sa1[k]*sa1[l];//����ֵ���ʱ�����ۼ�
					}
			for( m=ija1[ija1[j]];m<ija1[ija1[j]+1];m++)//�к�Ϊija[j],��Ϊ�к�
				if(ija1[m]==i)
					break;
			sa1[m]=(sa1[i]!=0? sum/sa1[i]:sum);                                                                                                                                                                                                        
		}

	}	
	return true;
}

void PCG::matrix_SSOR(double * &sa,unsigned long * &ija,const double w,const unsigned long n,double * &sa1)
{
    unsigned long i,k;
    for(i=0;i<n;i++)
    {
        for(k=ija[i];k<ija[i+1] && ija[k]<i;k++)
                sa1[k]=+w*sa[k];
        sa1[i]=sa[i]+w*sa[i];
    }
    for(i=0;i<n;i++)
    {
        for(k=ija[i];k<ija[i+1] && ija[k]<i;k++)
               {
                    sa1[k]/=sqrt(sa[ija[k]]);
                    sa1[k]/=sqrt(w*(2-w));
                }    
        sa1[i]/=sqrt(sa[i]);
        sa1[i]/=sqrt(w*(2-w));
    }
}
////����dsprsax��sa ija ����������x,�����b��,n ��ʾԭ����Ϊn*n
void PCG::dsprsax(double *sa,unsigned long *ija,double *x,double *b ,unsigned long n)
{
	unsigned long i,k;
		for(i=0;i<n;i++){
			b[i]=sa[i]*x[i];
			for(k=ija[i];k<ija[i+1];k++)
				b[i] += sa[k]*x[ija[k]];
		}
}
void PCG::asolve(double *sa1,unsigned long *ija,double *x,double *b,unsigned long n)
{
	if(flag==0)
	{
		unsigned long i;
		for(i=0;i<n;i++)
			x[i]=(sa[i]!= 0.0?b[i]/sa[i]:b[i]);
	}
if(flag==1 || flag==2)//SSOR �� Cholsky��ͬ 
{	
	unsigned long i,k;
	double sum;
/*	for(i=0;i<n;i++)
	{
		sum=sa1[i]*sa1[i];
		for(k=ija[i];k<ija[i+1] && ija[k]<i;k++)
			sum+=sa1[k]*sa1[k];
		sa1[i]=sqrt(sum);
	}
*/
	for(i=0;i<n;i++)
	{
		sum=b[i];
		for(k=ija[i];k<ija[i+1] && ija[k]<i;k++)
				sum-=sa1[k]*x[ija[k]];
		x[i]=sum/sa1[i];
	}
	for(i=n-1;i>=0;i--)
	{
		x[i]=x[i]/sa1[i];
			for(k=ija[i];k<ija[i+1] && ija[k]<i;k++)
					x[ija[k]]-=sa1[k]*x[i];
	}
}

}

int PCG::linbcg(void)
{
	if(flag==1 )
	{
		if(nz==0)
			return 0;
			
		if(sb)
		  delete [] sb;
		  
		sb=new double[nz];
		for(unsigned long i=0;i<nz;i++)
			sb[i]=sa[i];
		
		if(!Choldc(sb,ija,n))
			return 0;
	}
	if(flag==2)
	{
		if(nz==0)
			return 0;
			
		if(sb)
		  delete [] sb;
		
		sb=new double[nz];
		const double w=1.3;//���ɳ�����ȡ1.3
		matrix_SSOR(sa,ija,w,n,sb);
	}
	
		unsigned long j;
		double ak,akden,bk,err;
		double *p = new double[n];
		double *r = new double[n];
		double *z = new double[n];
		double *w=new double[n];
		double rol;//��
		double roll;//��'
		double rol0;//��0
		//��r=Ax
		dsprsax(sa,ija,x,r,n);
		for(j=0;j<n;j++)     //��r=b-Ax
			r[j] = b[j]-r[j];
		//z0��ֵ Mz=r
		asolve(sb,ija,z,r,n);
		for(j=0;j<n;j++)//p1=z0
			p[j] = z[j];
		for(rol=0.0,j=0;j<n;j++)//��0=rt*z0
			rol+= z[j]*r[j];
		rol0=rol;
		roll=rol;
		iter=1;
		while(iter<=itmax)
		{
			//atimes(n,p,w,0);//��w=Ap
			dsprsax(sa,ija,p,w,n);
			for(akden=0.0,j=0;j<n;j++) //��=��/pt*w
				akden += p[j]*w[j];
			ak = roll/akden;

			for(j=0;j<n;j++)//x=x+��p  r=r-��w
			{
				x[j] += ak*p[j];
				r[j] -= ak*w[j];
			}

			//asolve(n,r,z,0);//���Mz=r
			asolve(sb,ija,z,r,n);
							//���rol=rt*z
			for(rol=0.0,j=0;j<n;j++)
				rol+= z[j]*r[j];
				bk = rol/roll;//��=��/��'
			for(j=0;j<n;j++)//p=��*z+p
				p[j] = bk*p[j]+z[j];
		
			roll =rol;
			err=rol/rol0;//��=��/��0
			if (err <= tol){
			//	fprintf(stderr,"���������ܵ�������Ϊ�� %d\n",iter);
				break;
			}
			else
				iter++;
		}
//		if(iter>itmax)
//			;
//				fprintf(stderr,"����������������������ϣ�\n");
//		fprintf(stderr,"�����ݶȷ�������\n");
				//MessageBox("����������������������ϣ�",NULL,MB_OK);
		delete [] p;
		delete [] r;
		delete [] z;
		delete [] w;
		
		return iter;
}

