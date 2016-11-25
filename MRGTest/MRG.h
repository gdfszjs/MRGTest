#pragma once
#include <LibIV/libiv.h>
#include <string>
#include <vector>
using namespace std;

typedef LibIV::Math::v2d Range;
typedef LibIV::Math::v2d Rmessage;


typedef vector<vector<vector<int>>> LevelRangeNodeNumber;



//the number of the cutrange must be the of 2;
const static int rangenumber = 3;
//the weighting of the area and length parameters
const static int w = 0.5;

struct Location_Message
{
	int level_num;
	int range_num;
	int node_num;
};


struct MRGNode
{
	//the range that the node in
	Range range;
	//the node message
	Rmessage rmessage;

	int adjacentNodeNumber;
	Location_Message ** adjacentNodes;

	Location_Message * parentNode;

	int childNodeNumber;
	Location_Message ** childNodes;

	//record which pattern that node belongs to
	int treenum;
	//record which part that node belongs to
	int partnum;
	//record the branch
	MRGNode* Xnum = NULL;

	//record which level that node belongs to
	int nlevel;
	//record which range that node belongs to
	int nrange;
	//record which number that node belongs to
	int nnumber;
};

struct NodePair
{
	MRGNode * m;
	MRGNode * n;
};

typedef LibIV::Memory::Array::Array1D<MRGNode**> MRFGraph;
typedef LibIV::Memory::Array::Array1D<MRGNode***> MRFGraph2;
typedef LibIV::Memory::Array::Array1D<Location_Message*> NodeCandidate;
typedef vector<MRGNode*> MatchingFactory;

struct MRGPattern
{
	//the whole reeb graph graph2[a][b][c]
	//a-- the graph level
	//b-- the graph range
	//c-- the graph node
	MRFGraph2 graph2;
	//the node number in all level and range
	LevelRangeNodeNumber levelrangenodenumber;
	//the candidate node array for the corse graph creation
	NodeCandidate node_array;
	//the range that the graph has for every level
	vector<int> level_array;
	//the node number in one graph creation process
	int node_number = 0;

};



class MRG
{
public:
	MRG();
	~MRG();

	/*
		function fot the graph creation process
	*/
	//initialize the MRG
	void initialize(MRGPattern *e, string es);
	// allocate Memory For Graph2
	void allocateMemoryForGraph2(MRGPattern *e, string es);
	// free Memory Of Graph
	void freeMemoryOfGraph();
	// create the exact MRGNode in one range of the finest level
	void creatTestGraphinLevel(int range, MRGPattern *e, string es);
	// create the finest level through all range
	void creatTestAllGraph(MRGPattern *e, string es);
	// create other level through the lower level
	void creatTestOtherGraph(int level, MRGPattern *e, string es);
	// change the finest child and the corest parent
	void changeFinestChildandCorestParent(MRGPattern *e);
	// show the MRGNode that we create in one level
	void showMRGCaandidateArray(MRGPattern *e);
	// show MRGGaph&Node in exact Level and exact Range
	void showMRGGaphinLevelandRange(MRGPattern *e, int l, int r, int n);
	// init GraphNode with Parent
	void initGraphNodewithParent(MRGNode * e, v2d range, v2d message, int adjacentNodeNumber, Location_Message ** adjacentNodes, int treenum, int partnum, Location_Message * parent, int nlevel, int nrange, int nnumber);
	// init GraphNode with Child
	void initGraphNodewithChild(MRGNode * e, v2d range, v2d message, int adjacentNodeNumber, Location_Message ** adjacentNodes, int treenum, int partnum, int childnum,Location_Message ** child, int nlevel, int nrange, int nnumber);
	// init Finest GraphNode
	void initFinestGraphNode(MRGNode * e, v2d range, v2d message, int adjacentNodeNumber, Location_Message ** adjacentNodes, int treenum, int partnum, int nlevel, int nrange, int nnumber);
	// show the levelrangenodenumber in one MRGPattern
	void showLevelRangeNodeNumber(MRGPattern *e);
	//show all the message in a node
	void showAllMesasgeinNode(ofstream f1, MRGNode *e);
	//check MRGPattern
	void CheckMRGPattern(MRGPattern *e, string ss);

	/*
		function fot the graph comparation
	*/
	void matchingAlgorithm(MRGPattern * a, MRGPattern * b);
	double sim(Rmessage m, Rmessage n);
	double loss(Rmessage m, Rmessage n);
	bool com(MRGNode *m, MRGNode *n);
	Rmessage adj(Range c, MRGPattern *e, MRGNode *m);
	void initialization(MRGPattern *a, MRGPattern *b);
	void matching(MRGPattern * a, MRGPattern * b);
	NodePair* findmatchingpair(MRGPattern * a, MRGPattern * b);
	bool compareMLIST(MRGNode *m, MRGNode *n);
	bool findParentinMLIST(MRGNode *m, MRGNode *n);
	void spreadMLIST(MRGPattern *a, MRGPattern *b, MRGNode *m, MRGNode *n);
	void extendMLISTinOneDirection(MRGPattern * a,MRGNode *m,int direcion);
	int mat(MRGPattern *a, MRGPattern *b,MRGNode *m, MRGNode *n);
	int calculateResult();
	void printNLIST();

private:

	MatchingFactory NLIST;
	vector<NodePair*> MPAIR;
};

