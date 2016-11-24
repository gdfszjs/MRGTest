#include "MRG.h"

#include <sstream>
#include <fstream>
#include <iostream>
#include <algorithm>


using namespace std;

MRG::MRG()
{

}

MRG::~MRG()
{
	freeMemoryOfGraph();
}


/*
	function fot the graph creation process
*/

void MRG::initialize(MRGPattern *e, string es)
{
	allocateMemoryForGraph2(e, es);
	creatTestAllGraph(e, es);
	//showMRGCaandidateArray();

	for (int i = 1; i < e->level_array.size(); i++)
	{
		vector<vector<int>> temp_level_range;
		e->levelrangenodenumber.push_back(temp_level_range);
		creatTestOtherGraph(i, e, es);
	}

	changeFinestChildandCorestParent(e);
}

void MRG::allocateMemoryForGraph2(MRGPattern *e, string es)
{
	//init the whole finest MRGGraph
	for (int i = 0; i < rangenumber; i++)
	{
		e->level_array.push_back((2 << (rangenumber - i - 1)));
	}
	e->level_array.push_back(1);
	e->graph2.set(rangenumber + 1);
	e->graph2.fill(0);

	//the number og the finest graph node
	int finest_node_num = 0;

	e->graph2[0] = new MRGNode**[2 << (rangenumber - 1)];

	stringstream ss;
	string line;
	ss << es << "rangenum.txt";
	ifstream inFile(ss.str());

	while (getline(inFile, line))
	{
		char * strc = new char[strlen(line.c_str()) + 1];
		strcpy_s(strc, strlen(line.c_str()) + 1, line.c_str());
		char  *p = NULL, *pNext = NULL;

		p = strtok_s(strc, " ", &pNext);
		int index = 0;
		while (p != NULL)
		{
			e->graph2[0][index++] = new MRGNode*[atoi(p)];
			finest_node_num = finest_node_num + atoi(p);
			p = strtok_s(NULL, " ", &pNext);
		}

	}

	//init the candidate node array for the corse graph
	//cout << finest_node_num << endl;
	e->node_array.set(finest_node_num);
	for (int i = 0; i < finest_node_num; i++)
	{
		e->node_array[i] = new Location_Message;
	}
}

void MRG::freeMemoryOfGraph()
{

}

void MRG::creatTestGraphinLevel(int range, MRGPattern *e, string es)
{
	//record the node number int one range


	int graph_level = e->level_array.at(0) - (range + 1);

	//read the graph node message from the file
	stringstream ss;
	string line;
	ss << es << "0 " << range << " " << (range + 1) << ".txt";
	//cout << ss.str() << endl;

	//get the node number and init the rangearray
	ifstream inFile(ss.str());
	getline(inFile, line);
	istringstream iss(line);
	int nodenum;
	iss >> nodenum;

	//set the Node number in the selected range;
	int i = 0;

	//read the node message to create the reeb graph
	while (getline(inFile, line))
	{
		//MRGNode's adjanceNode index;
		vector<int> adjancenum;

		//MRGNode's area and length
		v2d as_l = _v2d_(0.0, 0.0);

		//MRGNode's treeindex and partindex
		int treeindex;
		int partindex;

		//read time->flag to the data message
		int index = 1;

		//read the data
		char * strc = new char[strlen(line.c_str()) + 1];
		strcpy_s(strc, strlen(line.c_str()) + 1, line.c_str());
		char  *p = NULL, *pNext = NULL;
		p = strtok_s(strc, " ", &pNext);
		while (p != NULL)
		{
			//the arae and length message
			if (index == 1 || index == 2)
			{
				as_l[index - 1] = atof(p);
			}
			//the treeindex message
			else if (index == 3)
			{
				treeindex = atoi(p);
			}
			//the partindex message
			else if (index == 4)
			{
				partindex = atoi(p);
			}
			//the adjance node number
			else
			{
				adjancenum.push_back(atoi(p));
			}
			p = strtok_s(NULL, " ", &pNext);
			index++;
		}
		//init the adjancenode message array
		Location_Message ** tempa = new Location_Message*[adjancenum.size()];
		for (int j = 0; j < adjancenum.size(); j++)
		{
			if (adjancenum.at(j) < 0)
			{
				//tempa[j] = graph2[0][graph_level - 1][0 - adjancenum.at(j)];
				tempa[j] = new Location_Message;
				tempa[j]->level_num = 0;
				tempa[j]->range_num = graph_level - 1;
				tempa[j]->node_num = 0 - adjancenum.at(j) - 1;
			}
			else
			{
				//tempa[j] = graph2[0][graph_level + 1][adjancenum.at(j)];
				tempa[j] = new Location_Message;
				tempa[j]->level_num = 0;
				tempa[j]->range_num = graph_level + 1;
				tempa[j]->node_num = adjancenum.at(j) - 1;
			}
		}
		//set the graph node
		e->graph2[0][graph_level][i] = new MRGNode;
		initFinestGraphNode(e->graph2[0][graph_level][i], _v2d_(1.0 / e->level_array.at(0) * range, 1.0 / e->level_array.at(0) * (range + 1)), as_l, adjancenum.size(), tempa, treeindex, partindex, 0, graph_level, i);

		//put the graph into the candidate array for the corser graph
		e->node_array[e->node_number]->range_num = graph_level;
		e->node_array[e->node_number]->node_num = i;

		e->node_number++;
		i++;
	}

	e->levelrangenodenumber.at(0).at(e->level_array.at(0) - (range + 1)).push_back(i);
	showMRGGaphinLevelandRange(e,0, graph_level, nodenum);
}

void MRG::creatTestAllGraph(MRGPattern *e, string es)
{
	vector<vector<int>> temp_level_range;
	e->levelrangenodenumber.push_back(temp_level_range);
	for (int i = 0; i < e->level_array.at(0); i++)
	{
		vector<int> range_node_number;
		e->levelrangenodenumber.at(0).push_back(range_node_number);
		creatTestGraphinLevel(e->level_array.at(0) - (i + 1), e, es);
	}
}

void MRG::creatTestOtherGraph(int level, MRGPattern *e, string es)
{

	//record the node that we create in one process
	int new_node_number = 0;

	vector<vector<MRGNode*>> temp_node_array;
	for (int i = 0; i < e->level_array.at(level); i++)
	{
		vector<MRGNode*> a;
		temp_node_array.push_back(a);
	}
	vector<v2d> selected_node;
	int quit_flag = 0;

	//showMRGCaandidateArray(e);

	for (int i = 0; i < e->node_array.size(); i++)
	{
		//get rid of the selected node and walk through it
		for (int j = 0; j < selected_node.size(); j++)
		{
			if (e->node_array[i]->range_num == selected_node.at(j)[0] && e->node_array[i]->node_num == selected_node.at(j)[1])
			{
				quit_flag = 1;
				break;
			}
			else;
		}
		if (quit_flag == 1)
		{
			quit_flag = 0;
			continue;
		}
		//store the selected node as the child node
		vector<MRGNode*> child_array;
		//the selected node
		//cout << "selected node:" << node_array[i]->range_num << node_array[i]->node_num << endl;
		MRGNode *p = e->graph2[0][e->node_array[i]->range_num][e->node_array[i]->node_num];
		selected_node.push_back(_v2d_(e->node_array[i]->range_num, e->node_array[i]->node_num));
		//new range
		int new_range = e->node_array[i]->range_num;
		//new message
		int new_message_x = p->rmessage[0];
		int new_message_y = p->rmessage[1];
		//new treenum and part num
		int treenum = p->treenum;
		int partnum = p->partnum;
		//create a new Node as a parent node
		Location_Message *parentLM = new Location_Message;
		MRGNode *parent = new MRGNode;
		p->parentNode = parentLM;
		//append the p Node into the child_array
		child_array.push_back(p);
		//find the adjance candidate
		//0 means that you need to find your partner in the upper range	
		if (new_range % 2 == 0)
		{
			for (int j = 0; j < p->adjacentNodeNumber; j++)
			{
				if (p->adjacentNodes[j]->range_num % 2 == 1 && p->adjacentNodes[j]->range_num > new_range)
				{
					MRGNode *p_can = e->graph2[0][p->adjacentNodes[j]->range_num][p->adjacentNodes[j]->node_num];
					selected_node.push_back(_v2d_(p->adjacentNodes[j]->range_num, p->adjacentNodes[j]->node_num));
					//calculate the new message for the parent node
					new_message_x = new_message_x + p_can->rmessage[0];
					new_message_y = new_message_y + p_can->rmessage[1];
					p_can->parentNode = parentLM;
					//append the node into the child_array
					child_array.push_back(p_can);
				}
			}
			new_range = (e->level_array.at(level - 1) - new_range - 1) / 2;
			Location_Message** child = new Location_Message*[child_array.size()];
			for (int i = 0; i < child_array.size(); i++)
			{
				child[i] = new Location_Message;
				child[i]->level_num = child_array.at(i)->nlevel;
				child[i]->range_num = child_array.at(i)->nrange;
				child[i]->node_num = child_array.at(i)->nnumber;
			}
			initGraphNodewithChild(parent, _v2d_(1.0 / e->level_array.at(level) * new_range, 1.0 / e->level_array.at(level) * (new_range + 1)), _v2d_(new_message_x, new_message_y), 0, NULL, p->treenum, p->partnum, child, level, e->level_array.at(level) - new_range - 1, temp_node_array.at(e->level_array.at(level) - new_range - 1).size());
			//get a new node
			new_node_number++;
			parentLM->level_num = level;
			parentLM->range_num = e->level_array.at(level) - new_range - 1;
			parentLM->node_num = temp_node_array.at(e->level_array.at(level) - new_range - 1).size();
			temp_node_array.at(e->level_array.at(level) - new_range - 1).push_back(parent);
		}
		//1 means that you need to find your partner in the lower range	
		else
		{
			for (int j = 0; j < p->adjacentNodeNumber; j++)
			{
				if (p->adjacentNodes[j]->range_num % 2 == 0 && p->adjacentNodes[j]->range_num < new_range)
				{
					MRGNode *p_can = e->graph2[0][p->adjacentNodes[j]->range_num][p->adjacentNodes[j]->node_num];
					selected_node.push_back(_v2d_(p->adjacentNodes[j]->range_num, p->adjacentNodes[j]->node_num));
					//calculate the new message for the parent node
					new_message_x = new_message_x + p_can->rmessage[0];
					new_message_y = new_message_y + p_can->rmessage[1];
					p_can->parentNode = parentLM;
					//append the node into the child_array
					child_array.push_back(p_can);
				}
			}
			new_range = new_range / 2;
			Location_Message** child = new Location_Message*[child_array.size()];
			for (int i = 0; i < child_array.size(); i++)
			{
				child[i] = new Location_Message;
				child[i]->level_num = child_array.at(i)->nlevel;
				child[i]->range_num = child_array.at(i)->nrange;
				child[i]->node_num = child_array.at(i)->nnumber;
			}
			initGraphNodewithChild(parent, _v2d_(1.0 / e->level_array.at(level) * new_range, 1.0 / e->level_array.at(level) * (new_range + 1)), _v2d_(new_message_x, new_message_y), 0, NULL, p->treenum, p->partnum, child, level, e->level_array.at(level) - new_range - 1, temp_node_array.at(e->level_array.at(level) - new_range - 1).size());
			//get a new node
			new_node_number++;
			parentLM->level_num = level;
			parentLM->range_num = e->level_array.at(level) - new_range - 1;
			parentLM->node_num = temp_node_array.at(e->level_array.at(level) - new_range - 1).size();
			temp_node_array.at(e->level_array.at(level) - new_range - 1).push_back(parent);
		}
	}

	//renew the node array for the next graph creation
	e->node_array.erase();
	//cout << "new_node_number:" << new_node_number << endl;
	e->node_array.set(new_node_number);
	for (int i = 0; i < new_node_number; i++)
	{
		e->node_array[i] = new Location_Message;
	}

	//create adjance-relationship between the created node
	int node_num = 0;
	e->graph2[level] = new MRGNode**[e->level_array.at(level)];
	for (int i = 0; i < temp_node_array.size(); i++)
	{

		e->graph2[level][i] = new MRGNode*[temp_node_array.at(i).size()];
		vector<int> range_node_number;
		range_node_number.push_back(temp_node_array.at(i).size());
		e->levelrangenodenumber.at(level).push_back(range_node_number);
		for (int j = 0; j < temp_node_array.at(i).size(); j++)
		{
			//cout << i << " " << j << endl;
			e->graph2[level][i][j] = temp_node_array.at(i).at(j);
			e->node_array[node_num + j]->range_num = i;
			e->node_array[node_num + j]->node_num = j;
		}
		node_num = node_num + temp_node_array.at(i).size();
		showMRGGaphinLevelandRange(e, level, i, temp_node_array.at(i).size());
	}

	//find the adjance node for every node
	vector<Location_Message*> adjance_node;
	for (int i = 0; i < temp_node_array.size(); i++)
	{
		for (int j = 0; j < temp_node_array.at(i).size(); j++)
		{
			MRGNode *p = e->graph2[level][i][j];

			//cout << p->partnum << endl;

			//find the adjance node in upper range
			if (i < (temp_node_array.size() - 1))
			{
				for (int k = 0; k < temp_node_array.at(i + 1).size(); k++)
				{
					//cout << "candidate:" << i + 1 << k << graph2[level][i + 1][k]->partnum << endl;

					if (e->graph2[level][i + 1][k]->partnum == p->partnum)
					{
						Location_Message *lm = new Location_Message;
						lm->range_num = i + 1;
						lm->node_num = k;
						adjance_node.push_back(lm);
					}
				}
			}
			//find the adjance node in lower range
			if (i > 0)
			{
				for (int k = 0; k < temp_node_array.at(i - 1).size(); k++)
				{
					//cout << "candidate:" << i - 1 << k << graph2[level][i - 1][k]->partnum << endl;

					if (e->graph2[level][i - 1][k]->partnum == p->partnum)
					{
						Location_Message *lm = new Location_Message;
						lm->range_num = i - 1;
						lm->node_num = k;
						adjance_node.push_back(lm);
					}
				}
			}
			//set the adjance node to every selected node
			Location_Message ** tempa = new Location_Message*[adjance_node.size()];
			for (int k = 0; k < adjance_node.size(); k++)
			{
				tempa[k] = adjance_node.at(k);
			}
			p->adjacentNodes = tempa;
			adjance_node.clear();

		}
	}


}

void MRG::changeFinestChildandCorestParent(MRGPattern * e)
{
	//the nodes in finest level have no childs
	for (int i = 0; i < e->levelrangenodenumber.at(0).size(); i++)
	{
		for (int j = 0; j < e->levelrangenodenumber.at(0).at(i).at(0); j++)
		{
			MRGNode *p = e->graph2[0][i][j];
			p->childNodeNumber = 0;
			p->childNodes = NULL;
		}
	}

	//the nodes in corest level have no parents
	for (int i = 0; i < e->levelrangenodenumber.at(e->levelrangenodenumber.size() - 1).size(); i++)
	{
		for (int j = 0; j < e->levelrangenodenumber.at(0).at(i).at(0); j++)
		{
			MRGNode *p = e->graph2[0][i][j];
			p->parentNode = NULL;
		}
	}
}

void MRG::showMRGCaandidateArray(MRGPattern *e)
{
	for (int i = 0; i < e->node_array.size(); i++)
	{
		cout << "Caandidate Node: " << e->node_array[i]->range_num << e->node_array[i]->node_num << endl;
	}
}

void MRG::showMRGGaphinLevelandRange(MRGPattern *e, int l, int r, int n)
{
	for (int i = 0; i < n; i++)
	{
		//cout << "range: " << e->graph2[l][r][i]->range[0] << " " << e->graph2[l][r][i]->range[1] << endl;
		cout << "message: " << e->graph2[l][r][i]->rmessage[0] << " " << e->graph2[l][r][i]->rmessage[1] << endl;
		cout << "treenum&partnum:" << e->graph2[l][r][i]->treenum << " " << e->graph2[l][r][i]->partnum << endl;
		cout << "nlevel&&nrange&&nnumber:" << e->graph2[l][r][i]->nlevel << " " << e->graph2[l][r][i]->nrange << " " << e->graph2[l][r][i]->nnumber << endl;
	}
}

void MRG::initGraphNodewithParent(MRGNode * e, v2d range, v2d message, int adjacentNodeNumber, Location_Message** adjacentNodes, int treenum, int partnum, Location_Message* parent, int nlevel, int nrange, int nnumber)
{
	//e = new MRGNode;
	e->range = range;
	e->rmessage = message;
	e->adjacentNodeNumber = adjacentNodeNumber;
	e->adjacentNodes = adjacentNodes;
	e->parentNode = parent;
	//e->childNodes = NULL;
	e->treenum = treenum;
	e->partnum = partnum;
	e->nlevel = nlevel;
	e->nrange = nrange;
	e->nnumber = nnumber;
}

void MRG::initGraphNodewithChild(MRGNode * e, v2d range, v2d message, int adjacentNodeNumber, Location_Message** adjacentNodes, int treenum, int partnum, Location_Message** child, int nlevel, int nrange, int nnumber)
{
	//e = new MRGNode;
	e->range = range;
	e->rmessage = message;
	e->adjacentNodeNumber = adjacentNodeNumber;
	e->adjacentNodes = adjacentNodes;
	//e->parentNode = NULL;
	e->childNodes = child;
	e->treenum = treenum;
	e->partnum = partnum;
	e->nlevel = nlevel;
	e->nrange = nrange;
	e->nnumber = nnumber;
}

void MRG::initFinestGraphNode(MRGNode * e, v2d range, v2d message, int adjacentNodeNumber, Location_Message** adjacentNodes, int treenum, int partnum, int nlevel, int nrange, int nnumber)
{
	//e = new MRGNode;
	e->range = range;
	e->rmessage = message;
	e->adjacentNodeNumber = adjacentNodeNumber;
	e->adjacentNodes = adjacentNodes;
	e->parentNode = NULL;
	e->childNodes = NULL;
	e->treenum = treenum;
	e->partnum = partnum;
	e->nlevel = nlevel;
	e->nrange = nrange;
	e->nnumber = nnumber;
}

void MRG::showLevelRangeNodeNumber(MRGPattern *e)
{
	for (int i = 0; i < e->levelrangenodenumber.size(); i++)
	{
		cout << "graph " << i << " :" << endl;
		for (int j = 0; j < e->levelrangenodenumber.at(i).size(); j++)
		{
			cout << "	range " << j << ": " << e->levelrangenodenumber.at(i).at(j).at(0) << " nodes" << endl;
		}
	}
	cout << endl;
}

void MRG::matchingAlgorithm(MRGPattern * a, MRGPattern * b)
{
	initialization(a, b);
	printNLIST();
	matching(a, b);
	//calculateResult();
}

/*
	function fot the graph comparation
*/
double MRG::sim(Rmessage m, Rmessage n)
{
	return 1.0 * w * max(m[0], n[0]) + 1.0 * (1 - w) * max(m[1], n[1]);
}

double MRG::loss(Rmessage m, Rmessage n)
{
	return  0.5 * (sim(m, m) + sim(n, n)) - sim(m, n);
}

bool MRG::com(MRGNode * m, MRGNode * n)
{
	if (m->nlevel == n->nlevel)
	{
		if (m->nrange == n->nrange)
		{
			if (m->nnumber == n->nnumber)
			{
				return true;
			}
			else return false;
		}
		else return false;
	}
	else return false;
}

Rmessage MRG::adj(Range c, MRGPattern *e,MRGNode *m)
{
	int rm1 = 0;
	int rm2 = 0;
	for (int i = 0; i < m->adjacentNodeNumber; i++)
	{
		MRGNode * p = e->graph2[m->adjacentNodes[i]->level_num][m->adjacentNodes[i]->range_num][m->adjacentNodes[i]->node_num];
		if (p->range[0] == c[0] && p->range[1] == c[1])
		{
			for (int j = 0; j < NLIST.size(); j++)
			{
				MRGNode * q = NLIST.at(j);
				if (com(p, q))
				{
					rm1 += q->rmessage[0];
					rm2 += q->rmessage[1];
				}
			}
		}
	}
	return _v2d_(rm1,rm2);
}

void MRG::initialization(MRGPattern * a, MRGPattern * b)
{
	cout << "put a nodes" << endl;
	for (int i = 0; i < a->levelrangenodenumber.at(a->levelrangenodenumber.size() - 1).at(0).at(0); i++)
	{
		NLIST.push_back(a->graph2[a->levelrangenodenumber.size() - 1][0][i]);
	}
	cout << "put b nodes" << endl;
	for (int i = 0; i < b->levelrangenodenumber.at(b->levelrangenodenumber.size() - 1).at(0).at(0); i++)
	{
		NLIST.push_back(b->graph2[b->levelrangenodenumber.size() - 1][0][i]);
	}
}

void MRG::matching(MRGPattern * a, MRGPattern * b)
{
	while (NLIST.size() != 0)
	{
		NodePair *s = findmatchingpair(a, b);
		MPAIR.push_back(s);
		if (s->m->childNodeNumber != 0)
		{
			for (int i = 0; i < s->m->childNodeNumber; i++)
			{
				MRGNode * p = a->graph2[s->m->childNodes[i]->level_num][s->m->childNodes[i]->range_num][s->m->childNodes[i]->node_num];
				NLIST.push_back(p);
			}
		}
		if (s->n->childNodeNumber != 0)
		{
			for (int i = 0; i < s->n->childNodeNumber; i++)
			{
				MRGNode * q = b->graph2[s->m->childNodes[i]->level_num][s->m->childNodes[i]->range_num][s->m->childNodes[i]->node_num];
				NLIST.push_back(q);
			}
		}
	}
}

NodePair* MRG::findmatchingpair(MRGPattern * a, MRGPattern * b)
{
	double result = 0;
	int select_index = 0;
	MRGNode * temp_m = NULL;
	MRGNode * selected_n = NULL;
	MRGNode * candidate_n = NULL;
	vector<int> temp_n;

	//find node m
	for (int i = 0; i < NLIST.size(); i++)
	{
		double temp_result = sim(NLIST.at(i)->rmessage, NLIST.at(i)->rmessage);
		if (temp_result > result)
		{
			result = temp_result;
			select_index = i;
		}
	}

	//remove node m
	temp_m = NLIST.at(select_index);
	std::vector<MRGNode*>::iterator it = NLIST.begin() + select_index;
	NLIST.erase(it);

	for (int i = 0; i < NLIST.size(); i++)
	{
		//比较条件有4个
		//1.m和n有same μn-range
		//2.m和n的MLIST相同 
		//3.m和n的parent也得匹配
		//4.m和n不属于同一个MRG
		candidate_n = NLIST.at(i);
		if (candidate_n->range[0] == temp_m->range[0] && candidate_n->range[1] == temp_m->range[1])
		{
			if (compareMLIST(temp_m,candidate_n))
			{
				if (candidate_n->treenum != temp_m->treenum)
				{
					if (findParentinMLIST(candidate_n, temp_m))
					{
						temp_n.push_back(i);
					}
				}
			}
		}
		
	}

	//match failed
	if (temp_n.size() == 0)
	{
		cout << "no matching Node!" << endl;
	}
	//only one matched
	else if (temp_n.size() == 1)
	{
		selected_n = NLIST.at(temp_n.at(0));

		it = NLIST.begin() + temp_n.at(0);
		NLIST.erase(it);
	}
	//lots of candidates
	else
	{
		int mat_result = 0;
		for (int i = 0; i < temp_n.size(); i++)
		{
			int temp_mat_result = mat(a,b,temp_m, NLIST.at(temp_n.at(i)));
			if (temp_mat_result > mat_result)
			{
				mat_result = temp_mat_result;
				select_index = i;
			}
		}
		selected_n = NLIST.at(temp_n.at(select_index));

		it = NLIST.begin() + select_index;
		NLIST.erase(it);
	}

	//remove node n
	

	//return NM pair
	NodePair *result_pair = new NodePair;
	result_pair->m = temp_m;
	result_pair->n = selected_n;
	extendMLIST(a,temp_m); 
	extendMLIST(b, selected_n);
	return result_pair;
}

bool MRG::compareMLIST(MRGNode * m, MRGNode * n)
{
	if (m->Xnum == NULL && n->Xnum == NULL)
	{
		return true;
	}
	if (findParentinMLIST(m->Xnum,n->Xnum))
	{
		return true;
	}
	else return false;
}

bool MRG::findParentinMLIST(MRGNode * m, MRGNode * n)
{
	if (m->parentNode == NULL && n->parentNode == NULL)
	{
		return true;
	}

	for (int i = 0; i < MPAIR.size(); i++)
	{
		MRGNode *c_m = MPAIR.at(i)->m;
		MRGNode *c_n = MPAIR.at(i)->n;
		if (m->treenum == c_m->treenum)
		{
			if (m->nlevel == c_m->nlevel && m->nrange == c_m->nrange && m->nnumber == c_m->nnumber)
			{
				if (n->treenum == c_n->treenum)
				{
					if (n->nlevel == c_n->nlevel && n->nrange == c_n->nrange && n->nnumber == c_n->nnumber)
					{
						return true;
					}
				}
			}
		}
	}
	return false;
}

void MRG::extendMLIST(MRGPattern * a,MRGNode * m)
{
	int nodenum = 0;
	int index = 0;
	//upper branch
	for (int i = 0; i < m->adjacentNodeNumber; i++)
	{
		if (m->adjacentNodes[i]->range_num == m->nrange + 1)
		{
			nodenum++;
			index = i;
			a->graph2[m->adjacentNodes[i]->level_num][m->adjacentNodes[i]->range_num][m->adjacentNodes[i]->node_num]->Xnum = m;
		}
	}
	if (nodenum == 1)
	{
		extendMLIST(a, a->graph2[m->adjacentNodes[index]->level_num][m->adjacentNodes[index]->range_num][m->adjacentNodes[index]->node_num]);
	}
	if (nodenum == 0 && nodenum > 1)
	{
		return;
	}
	//lower branch
	for (int i = 0; i < m->adjacentNodeNumber; i++)
	{
		if (m->adjacentNodes[i]->range_num == m->nrange - 1)
		{
			nodenum++;
			index = i;
			a->graph2[m->adjacentNodes[i]->level_num][m->adjacentNodes[i]->range_num][m->adjacentNodes[i]->node_num]->Xnum = m;
		}
	}
	if (nodenum == 1)
	{
		extendMLIST(a, a->graph2[m->adjacentNodes[index]->level_num][m->adjacentNodes[index]->range_num][m->adjacentNodes[index]->node_num]);
	}
	if (nodenum == 0 && nodenum > 1)
	{
		return;
	}
}

int MRG::mat(MRGPattern *a, MRGPattern *b,MRGNode * m, MRGNode * n)
{
	double LossMN = loss(m->rmessage, n->rmessage);
	double LossADJ = 0;
	for (int i = 0; i < 2; i++)
	{
		//upper range
		double r1 = m->range[0] + 1.0 / a->levelrangenodenumber.at(m->nlevel).size();
		double r2 = m->range[1] + 1.0 / a->levelrangenodenumber.at(m->nlevel).size();
		LossADJ += loss(adj(_v2d_(r1,r2), a, m), adj(_v2d_(r1, r2), b, n));
		//lower range
		r1 = m->range[0] - 1.0 / a->levelrangenodenumber.at(m->nlevel).size();
		r2 = m->range[1] - 1.0 / a->levelrangenodenumber.at(m->nlevel).size();
		LossADJ += loss(adj(_v2d_(r1, r2), a, m), adj(_v2d_(r1, r2), b, n));
	}
	return -LossMN - LossADJ;
}

int MRG::calculateResult()
{
	int result = 0;
	for (int i = 0; i < MPAIR.size(); i++)
	{
		result += sim(MPAIR.at(i)->m->rmessage, MPAIR.at(i)->n->rmessage);
	}
	return result;
}

void MRG::printNLIST()
{
	for (int i = 0; i < NLIST.size(); i++)
	{
		//cout << "range: " << e->graph2[l][r][i]->range[0] << " " << e->graph2[l][r][i]->range[1] << endl;
		cout << "message: " << NLIST.at(i)->rmessage[0] << " " << NLIST.at(i)->rmessage[1] << endl;
		cout << "treenum&partnum:" << NLIST.at(i)->treenum << " " << NLIST.at(i)->partnum << endl;
		cout << "nlevel&&nrange&&nnumber:" << NLIST.at(i)->nlevel << " " << NLIST.at(i)->nrange << " " << NLIST.at(i)->nnumber << endl;
	}
}
