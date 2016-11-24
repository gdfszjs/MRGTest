/*****************************
	
	�������߳�
	2008/3/30 by Yongwei Nie
	
	������2008/3/31

	1. ��CWinThread::m_bAutoDelete��Ϊfalse
	2. ResumeWorkerThread �� SuspendWorkerThread ���Ϸ���ֵ (long)
	3. GetExitCode()��ȡ�߳��˳���ԭ�����̲߳�δ����������STILL_ACTIVE
	
*****************************/
#pragma once

class CWorkerThread
{
public:
	typedef UINT (*WorkerThreadProc)(LPVOID pParam);
private:
	CWinThread * pthread;
	WorkerThreadProc workerproc;
	int nPriority;
	UINT nStackSize;
	DWORD dwCreateFlags;
	LPSECURITY_ATTRIBUTES lpSecurityAttrs;
	
public:
	explicit CWorkerThread(WorkerThreadProc _workerproc,int pri = THREAD_PRIORITY_NORMAL,UINT siz = 0,
		DWORD fla = CREATE_SUSPENDED,LPSECURITY_ATTRIBUTES attrs = NULL);
	void BeginWorkerThread(LPVOID pParam);
    long ResumeWorkerThread();
	long SuspendWorkerThread();
	long GetExitCode();
	~CWorkerThread(void);
};
