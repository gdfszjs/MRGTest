#include <iostream>
#include <vector>
#include <LibIV/libiv.h>
#include <iostream>
#include "MRG.h"

using namespace std;

int main()
{
	MRG mrg;
	//one parrtern
	cout << "creating MRG1..."<< endl;
	MRGPattern *a = new MRGPattern;
	string e = "D:\\MRGFIle\\1\\1\\";
	mrg.initialize(a, e);
	mrg.CheckMRGPattern(a,e);
	//mrg.showLevelRangeNodeNumber(a);

	//the other pattern
	cout << "creating MRG2..." << endl;
	MRGPattern *b = new MRGPattern;
	e = "D:\\MRGFIle\\2\\1\\";
	mrg.initialize(b, e);
	mrg.CheckMRGPattern(b, e);
	//mrg.showLevelRangeNodeNumber(b);

	//match two pattern
	cout << "matching algorithm..." << endl;
	mrg.matchingAlgorithm(a, b);
	system("pause");
}