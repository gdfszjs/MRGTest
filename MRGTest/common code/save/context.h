#pragma  once
#include <stdio.h>

// ��ȡ����ȫ�ֱ���
class Context
{
public:
	//enum MODE { RENDER_MESHDATA3D_POINT = 0,RENDER_MESHDATA3D_MESH,RENDER_VIDEO_2D};
	enum MODE { RENDER_VIDEO_STREAM_SOURCE = 0, RENDER_VIDEO_STREAM_TUBE, RENDER_TUBE, RENDER_TUBE_MANAGER};
	Context()
	{
		_anglex = 0;
		_angley = 0;
		_tx = _ty = _tz = 0;
		_sx = _sy = _sz = 1;
		_mesh_vis_pos = 0;
		_winx = 0;
		_winy = 0;
		_rightButtonInteractionType = 0;
		//memset(_cursorPos,0,sizeof(char)*100);
		setCursorPos(0,0);
		_frame = 0;
		//_mode = RENDER_MESHDATA3D_POINT;
		_mode = RENDER_VIDEO_STREAM_SOURCE;
	}
	/*Context(double anglex,double angley, double tz, double mesh_vis_pos,int winx,int winy,int rbit) :
		_anglex(anglex),_angley(angley),_tz(tz),_mesh_vis_pos(mesh_vis_pos),_winx(winx),_winy(winy),_rightButtonInteractionType(rbit)
	{
	
	}*/
	~Context()
	{

	}
	
	double getAnglex()	   { return _anglex; }
	double getAngley()	   { return _angley; };
	
	double getTx()		   { return _tx; };
	double getTy()		   { return _ty; };
	double getTz()		   { return _tz; };
		
	double getSx()		   { return _sx; };
	double getSy()		   { return _sy; };
	double getSz()		   { return _sz; };


	double getMeshVisPos() { return _mesh_vis_pos; }
	int    getWinX()       { return _winx; }
	int    getWinY()       { return _winy; }
	int    getRightButtonInteractionType() {return _rightButtonInteractionType;}
	char*  getCursorPos()  { return _cursorPos; }
	int    getFrame()      { return _frame; }
	MODE    getMode()       { return _mode; }
	
	void setAnglex(double d)     { _anglex = d; }
	void setAngley(double d)     { _angley = d; }
	
	void setTx(double d)         { _tx = d; }
	void setTy(double d)         { _ty = d; }
	void setTz(double d)         { _tz = d; }
	
	void setSx(double d)         { _sx = d; }
	void setSy(double d)         { _sy = d; }
	void setSz(double d)         { _sz = d; }

	void setMeshVisPos(double d) { _mesh_vis_pos = d; }
	void setWinX(int n)          { _winx = n; }
	void setWinY(int n)          { _winy = n; }
	/*void setRightButtonInteractionType(int t) { 
		if(_mode == RENDER_MESHDATA3D_POINT ||
			_mode == RENDER_MESHDATA3D_MESH)
			_rightButtonInteractionType = t; 
	}*/
	void setCursorPos(int x,int y) { sprintf_s(_cursorPos,100,"%d,%d",x,y); }
	void setFrame(int n)         { _frame = n; }
	//void setMode()          { _mode = MODE((_mode + 1)%(RENDER_VIDEO_2D+1)); }
	void setMode()          { _mode = MODE((_mode + 1)%(RENDER_TUBE_MANAGER+1)); }
	void setMode(MODE m)          { _mode = m; }


	void ToBe() { pos = be; }
	void ToEd() { pos = ed; }
	void SetBe(int be) { this->be = be; }
	void SetEd(int ed) { this->ed = ed; }
	int  GetPos()  { return pos; }
	void Next() { _check_(pos >= be && pos < ed); pos++; }
	void Back() { _check_(pos > be && pos <= ed); pos--; }


private:
	// ��ת
	double _anglex,_angley;
	
	// ƽ��
	double _tx,_ty,_tz;

	// ����
	double _sx,_sy,_sz;

	// ��ʾmesh��λ��
	double _mesh_vis_pos;

	// ���ڴ�С
	int _winx,_winy;

	// ��ǰ�������ͣ��ұ߽�/ָ���� ��
	int _rightButtonInteractionType;

	// ����ƶ�ʱ��λ��
	char _cursorPos[100];

	// ��ǰ֡ 
	int _frame;
	
	// ״̬
	MODE _mode;

	// ���������
	int pos;
	int be;
	int ed;
};