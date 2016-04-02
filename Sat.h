/*
 * Sat.h
 *
 *  Created on: Mar 6, 2016
 *      Author: RLi
 */

#ifndef SAT_H_
#define SAT_H_

#include <fstream>
#include <iostream>
#include <cstdlib>
#include <string>
#include <vector>
#include <sstream>
#include <algorithm>
#include <stack>
#include <iterator>

using namespace std;

class Sat {


    public:
	// function prototypes
    Sat();
    void printInfo();
    void printClauses();
    void updateNums(int varAssigned, int sign);
    vector<int> * getClause(int num);
    int reduceUnitClause();
    double pureVarSolve();
    double solveSSAT();
    int RAND();
    int MOMS();
    int MAX_SAT(vector<int> numPos, vector<int> numNeg);

    int returnUCP() {return UCP;}
    int returnPVE() {return PVE;}
    int returnVOS() {return VOS;}


    // variables
    vector<int> numActive;
    vector<int> assignment;
    vector<int> numPos;
    vector<int> numNeg;

    private:
    vector<double> varProb;
    vector<vector<int> > posClauses;
    vector<vector<int> > negClauses;
    vector<vector<int> > clauses;
    int UCP;
    int PVE;
    int VOS;

};

#endif /* SAT_H_ */
