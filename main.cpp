/*
 * main.cpp
 *
 *  Created on: Mar 6, 2016
 *      Author: RLi
 */

#include "Sat.h"

using namespace std;

int main() {

    // construct sat problem
    Sat sat;

    double prob = 1;

    // begin clock
	double start, end, seconds;
	start = clock();

	vector<int> numActive = sat.numActive;
    vector<int> numPos = sat.numPos;
    vector<int> numNeg = sat.numNeg;

    prob = sat.solveSSAT();

	// end clock
	end = clock();
	seconds = (end - start) / CLOCKS_PER_SEC;

	//sat.printInfo();

	cout << "Total UCP: " << sat.returnUCP() << endl;
	cout << "Total PVE: " << sat.returnPVE() << endl;
	cout << "Total VOS: " << sat.returnVOS() << endl;
	cout << "Success probability is " << prob << endl;
	cout << "Solution Time (CPU secs): " << seconds;

}
