
#pragma once
//���� 2008.12 Chu Yu.
class PCG //ǰ�����������ݶ� ��Ax=b��� ����AΪ��������Ϊ�Գơ������ľ��� 
{
private:
	double *sa,*sb,*x,*b;
	unsigned long *ija;								//What's this?
	double tol;//��������
	int itmax;//��������
	int iter;//ʵ�ʵ�������
	unsigned long n;//����߳� 
	int flag;									//What's this? Precondition factor
	unsigned long nz;//ϡ������з���Ԫ����Ŀ		
public:
	~PCG(void);
	//���캯�������б����£�
	//sa,ija ����ľ���A����ϡ����ʽ�洢        //How? sa is used to store values and ija is a auxiliary array.
	//x,�������Ľ��
	//b�����ұߵľ���b
	//n A������
	//sb ijb Ԥ��������
	// tol ��������
	// itmax ����������
	//iter ʵ�ʵ�������
	PCG(double * &sa_,unsigned long *&ija_,double *&x_,double *&b_,const unsigned long n_,
											const double tol_,const int itmax_);
	PCG(double * &sa_,unsigned long *&ija_,double *&x_,double *&b_,const unsigned long n_,
								const double tol_,const int itmax_,const int flag_,const int nz_);//��flagΪ0ʱ��sb��ijb һ��Ҫ��ֵ��
	
	//Ԥ��������
	//flag=0,ΪĬ�ϵ��ſɱȴ������ӣ��������Խ�Ԫ��ΪԤ��������
	//flag=2��Ϊ���ɳ�Ԥ��������,sa1 Ϊ�ֽ�Ľ�� 
	void matrix_SSOR( double * &sa, unsigned long * &ija,const double w,const unsigned long n,double * &sa1);	
	//flag=1 Ϊ����˹������ȫ�ֽ⣬�ֽ�Ľ��������sa1��
	bool Choldc(double *&sa1, unsigned long *& ija1,const unsigned long n);								

	//����PCG ����
	int linbcg();

	//linbcg ���õ��ĺ��� ����b=Ax
	//sa ija ΪA ��nΪ������
	void asolve(double *sa,unsigned long *ija,double *x,double *b, unsigned long n);

	///linbcg �õ��ĺ��� �����Mx=b ���Ϊx ����MΪԤ�������ӣ�Ϊ����ȫCholsky�ֽ�Ľ����Ϊ�����ǲ���
	void dsprsax(double *sa,unsigned long *ija,double *x,double *b ,unsigned long n);
};
