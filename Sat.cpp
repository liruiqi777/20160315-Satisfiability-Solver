/*
 * Sat.cpp
 *
 *  Created on: Mar 11, 2016
 *      Author: RLi
 */


#include "Sat.h"

using namespace std;

// constructor
// reads all relevant information from SSAT filename into data structures
Sat::Sat() {
	UCP = 0;
	PVE = 0;
	VOS = 0;

    string fileName;
    ifstream SSATStream;

    // open file
    do {
        cout << "Enter a SSAT problem file: ";
        cin >> fileName;
        SSATStream.open(fileName.c_str(), ios::in);
    } while (!SSATStream.good());

    string line;

    //check that file is not empty
    if (!getline(SSATStream, line)) {
        cerr << "See error: SSAT file is empty" << endl;
        return;
    }

    // skip block of comments
    getline(SSATStream, line);
    while (!line.empty()) {
        getline(SSATStream, line);
    }

    // skip clause and variable numbers
    getline(SSATStream, line);
    while (!line.empty()) {
        getline(SSATStream, line);
    }

    // skip next two lines to access variables
    getline(SSATStream, line);
    getline(SSATStream, line);

    int varNum;
    double variableProb;

    // push all probabilities into varProb vector
    while (!line.empty()) {
        istringstream iss(line);
        iss >> varNum;
        iss >> variableProb;
        if (variableProb == -1) {
            variableProb = 1;
        }
        varProb.push_back(variableProb);
        getline(SSATStream, line);
    }

    // initialize clause info with proper size
    int n = varProb.size();
    posClauses.resize(n);
    negClauses.resize(n);
    numPos.resize(n);
    numNeg.resize(n);

    // skip next two lines to access clauses
    getline(SSATStream, line);
    getline(SSATStream, line);

    // gather information from clauses
    int item;
    int index;
    int numVars = 0;
    int clauseCount = 0;
    vector<int> newClause;

    while (!line.empty()) {

        ++clauseCount;

        istringstream iss(line);
        iss >> item;

        // add the rest of the line to newClause
        while (item != 0) {

            // add to clause
            newClause.push_back(item);

            // add variable to Var structure
            index = abs(item) - 1;
            if (item > 0) {

                (posClauses.at(index)).push_back(clauseCount - 1);
                numPos.at(index) ++;

            } else {

                (negClauses.at(index)).push_back(clauseCount - 1);
                numNeg.at(index) ++;
            }

            iss >> item;
        }

        // determine number of variables in the clause
        numVars = newClause.size();
        numActive.push_back(numVars);   // all are unassigned

        // push newClause onto clauses vector
        clauses.push_back(newClause);

        // get next line of file
        getline(SSATStream, line);
        newClause.clear();

    }

    // initialize assignment vector
    for (int i = 0; i < varProb.size(); i++) {

        assignment.push_back(0);
    }
}

// prints info
void Sat::printInfo() {

    // print # active variables in each clause
    cout << "numActive:" << endl;
    for (int i = 0; i < numActive.size(); i++) {
        cout << numActive.at(i) << endl;
    }

    // print current assignments
    cout << endl;
    cout << "Current assignments:" << endl;
    for (int i = 0; i < assignment.size(); i++) {
        cout << assignment.at(i) << endl;
    }

    // print numPos
    cout << endl;
    cout << "numPos:" << endl;
    for (int i = 0; i < numPos.size(); i++) {
        cout << numPos.at(i) << endl;
    }

}

// returns pointer to single clause vector
vector<int> *Sat::getClause(int num) {

    return &(clauses.at(num));
}

// print clauses
void Sat::printClauses() {

    vector<int> * curr;
    for (int i = 0; i < clauses.size(); i++) {

        curr = getClause(i);
        for (int j = 0; j < curr->size(); j++) {

            cout << curr->at(j) << " ";
        }
        cout << endl;
    }
}

int Sat::reduceUnitClause(){
	// return the index of the unit variable
	for (int i = 0; i < numActive.size(); i ++){

		// if the clause has only one active variable
		if (numActive[i] == 1){
			vector<int> clause = clauses[i];

			// once a unit var is found
			// return the index of the var
			for (int j = 0; j < clause.size(); j++){
				int varIndex = abs(clause[j]) - 1;
				if (assignment[varIndex] == 0){
					return clause[j];
				}
			}
		}
	}
	// if no unit var exists, return 0
	return 0;
}

double Sat::pureVarSolve(){
	int varNum = 0;
	double pos = 0, neg = 0;

	// go through the assignment vector to find an unassigned variable
	for(vector<int>::iterator iter = assignment.begin();
			iter != assignment.end(); iter++, varNum++){

		// unassigned variable found
		if (*iter == 0 && varProb.at(varNum) == 1){

			// extract its pos and neg clauses
			vector<int> posClauseOfVar = posClauses.at(varNum);
			vector<int> negClauseOfVar = negClauses.at(varNum);
			int posCount = 0, negCount = 0;

			// go through the positive clauses and see if there are
			// any more active clauses the var is in
			for(vector<int>::iterator itPos = posClauseOfVar.begin();
					itPos != posClauseOfVar.end(); itPos++){

				// if there is, increase posCount and break the for loop
				if (numActive.at(*itPos) != 0){
					posCount ++;
					break;
				}
			}

			// go through the negative clauses and see if there are
			// any more active clauses the var is in
			for(vector<int>::iterator itNeg = negClauseOfVar.begin();
					itNeg != negClauseOfVar.end(); itNeg++){

				// if there is, increase negCount and break the for loop
				if (numActive.at(*itNeg) != 0){
					negCount ++;

					// if posCount and negCount are both bigger than 1
					// then it is not a pure var. go to the next
					// unassigned variable
					if (posCount != 0){
						goto label;
					}
					break;
				}
			}

			// backup the numActive
			vector<int> numActiveCopy = numActive;

			// set the variable positive
			if (posCount != 0){
				assignment.at(varNum) = 1;

				// adjust the numActive according to the positive assignment
				for (vector<int>::iterator iterPos = posClauseOfVar.begin();
						iterPos != posClauseOfVar.end(); iterPos++){
					numActive.at(*iterPos) = 0;
				}

				// call solveSSAT() recursively to find the satisfiability
				// when the variable is assigned positive
				pos = solveSSAT();
				// reassign numActive to its original values
				numActive = numActiveCopy;

				// if posCount is not 0 and negCount is 0, then
				// reset the assignment and return the probability
				if (negCount == 0){
					assignment.at(varNum) = 0;
					return pos;
				}
			}


			// Just set the var to negative
			assignment.at(varNum) = -1;

			// adjust the numActive according to the negative assignment
			for (vector<int>::iterator itNeg = negClauseOfVar.begin();
					itNeg != negClauseOfVar.end(); itNeg++){
				numActive.at(*itNeg) = 0;
			}

			// call solveSSAT() recursively to find the satisfiability
			// when the variable is assigned negative
			neg = solveSSAT();
			// reassign numActive to its original values
			numActive = numActiveCopy;

			// remove the assignment to the variable
			assignment.at(varNum) = 0;
			// if the variable only appears negatively
			// return the neg probability
			if (negCount != 0 && posCount == 0){
				return neg;
			}

			// if it can be both way
			return neg + pos;
		}

		// the label corresponding to the "goto"
		label:;

	}
	// return -1 if there is no pure variable
	return -1;
}

int Sat::MOMS(){

	// create an array to keep track of how many times
	// each variable appears
	int *A = (int*) malloc (assignment.size() * (sizeof(int)));
	for (int i = 0; i < assignment.size(); i++){
		A[i] = 0;
	}

	int varIndex;
	int clauseIndex = 0;
	// go through all the clauses that have less than 4
	// variables active and is not satisfied already
	for (vector<int>::iterator iter = numActive.begin();
			iter != numActive.end(); iter++, clauseIndex++){
		if (*iter < 4 && *iter != 0){

			// extract the clause
			vector<int> clause = clauses[clauseIndex];

			// go through the clauses and increase the count
			// for the corresponding variable
			for (vector<int>::iterator
				it = clause.begin(); it != clause.end(); it++){
				varIndex = abs(int(*it)) - 1;
				if (assignment[varIndex] == 0){
					A[varIndex] ++;
				}
			}
		}
	}
	int maxVar = 0;
	varIndex = -1;

	// find the variable that has appeared the most times
	for (int i = 0; i < assignment.size(); i++){
		if (A[i] > maxVar){
			maxVar = A[i];
			varIndex = i;
		}
	}

	return varIndex;
}

double Sat::solveSSAT(){

	//BRUTE FORCE METHOD
	//iterate through all variables
	for(vector<int>::iterator iter = numActive.begin();
			iter != numActive.end(); iter++){

		// if an unsatisfied clause is identified
		// break the loop and start solving
		// else return probability of 1
		if(*iter != 0){
			break;
		} else if (iter + 1 == numActive.end()){
			return 1;
		}
	}

	// back up the original numActive
	vector<int> numActiveCopy = numActive;

	double prob = 1.0;

//	// solve using pureVarSolve
//	prob = pureVarSolve();
//	// if the prob return is -1, then there is no pure
//	// variable. Otherwise return the probability
//	if (prob != -1){
//		PVE ++;
//		return prob;
//	}

	// find the unit variable (the index returned as it appears in the unit clause)
//	int unitVar = reduceUnitClause();
	// uncomment the following when choose not to use reduceUnitClause
	int unitVar = 0;

	if (unitVar != 0){
		UCP ++;
		int sign = 0;

		// if the unitVar appeared positively
		// then assign positive to that var
		if (unitVar > 0){
			sign = 1;
			unitVar = unitVar -1;
			assignment.at(unitVar) = 1;
			//remove all clauses where variable i appears positively
			for (int y = 0; y < posClauses.at(unitVar).size(); y++) {
				int clauseToReduce = posClauses.at(unitVar).at(y);
				numActive.at(clauseToReduce) = 0;
			}
			// if i appears negatively in a clause, reduce the number of
			// numActiveCopy by 1. If it is already 1 before reducing,
			// then the clause cannot be satisfied, so return 0
			for (int x = 0; x < negClauses.at(unitVar).size(); x++) {
				int clauseToReduce = negClauses.at(unitVar).at(x);
				if (numActive.at(clauseToReduce) == 1){
					numActive = numActiveCopy;
					assignment.at(unitVar) = 0;
					return 0;
				} else if (numActive.at(clauseToReduce) != 0){
					numActive.at(clauseToReduce) --;
				}
			}

			// call solveSSAT() recursively
			prob = solveSSAT();

			// reset numActive and assignment
			numActive = numActiveCopy;
			assignment.at(unitVar) = 0;

			return prob * varProb[unitVar];

		} else {
			// assign negative to the variable
			sign = -1;

			unitVar = -unitVar -1;
			assignment.at(unitVar) = -1;
			//remove all clauses where variable i appears negatively
			for (int y = 0; y < negClauses.at(unitVar).size(); y++) {
				int clauseToReduce = negClauses.at(unitVar).at(y);
				numActive.at(clauseToReduce) = 0;
			}
			// if i appears positively in a clause, reduce the number of
			// numActiveassginedCopy by 1. If it is already 1 before reducing,
			// then the clause cannot be satisfied, so return 0
			for (int x = 0; x < posClauses.at(unitVar).size(); x++) {
				int clauseToReduce = posClauses.at(unitVar).at(x);
				if (numActive.at(clauseToReduce) == 1){
					numActive = numActiveCopy;
					assignment.at(unitVar) = 0;
					return 0;
				} else if (numActive.at(clauseToReduce) != 0){
					numActive.at(clauseToReduce) --;
				}
			}

			// call solveSSAT recursively with the updated assignments
			prob = solveSSAT();

			//reset assignment and numActive
			numActive = numActiveCopy;
			assignment.at(unitVar) = 0;

			if (varProb.at(unitVar) == 1){
				return prob;
			} else{
				return prob * (1 - varProb[unitVar]);
			}
		}
	}

	// Using the Heuristics
//	int varNum = MAX_SAT(numPos, numNeg);
//	int varNum = RAND();
//  int varNum = MOMS();
	int varNum = -1;

	if (varNum == -1){
		varNum = 0;
		for (vector<int>::iterator iter = assignment.begin(); iter != assignment.end(); iter++, varNum++) {
			if (*iter == 0){
				break;
			}
		}
	}



	//make a copy of numActiveassigned

	VOS += 2;

	double neg = -1, pos = -1;

	assignment.at(varNum) = -1;
	//remove all clauses where variable i appears negatively
	for (int y = 0; y < negClauses.at(varNum).size(); y++) {
		int clauseToReduce = negClauses.at(varNum).at(y);
		numActive.at(clauseToReduce) = 0;
	}
	//if i appears positively in a clause, reduce the number of
	//numActiveassginedCopy by 1
	for (int x = 0; x < posClauses.at(varNum).size(); x++) {
		int clauseToReduce = posClauses.at(varNum).at(x);
		if (numActive.at(clauseToReduce) == 1){
			neg = 0;
			break;
		} else if (numActive.at(clauseToReduce) != 0){
			numActive.at(clauseToReduce) --;
		}
	}

	//make recursive call
	if (neg == -1){
		neg = solveSSAT();
	}

	//reset numActiveassigned
	numActive = numActiveCopy;

	assignment.at(varNum) = 1;
	//remove all clauses where variable i appears positively
	for (int y = 0; y < posClauses.at(varNum).size(); y++) {
		int clauseToReduce = posClauses.at(varNum).at(y);
		numActive.at(clauseToReduce) = 0;
	}
	//if i appears negatively in a clause, reduce the number of
	//numActiveassginedCopy by 1
	for (int x = 0; x < negClauses.at(varNum).size(); x++) {
		int clauseToReduce = negClauses.at(varNum).at(x);
		if (numActive.at(clauseToReduce) == 1){
			pos = 0;
			break;
		} else if (numActive.at(clauseToReduce) != 0){
			numActive.at(clauseToReduce) --;
		}
	}

	if (pos == -1){
		pos = solveSSAT();
	}
	numActive = numActiveCopy;

	//reset assignment
	assignment.at(varNum) = 0;

	if (varProb.at(varNum) == 1){
		if (pos > neg){
			return pos;
		}
		return neg;
	}

	return pos * varProb.at(varNum) + neg * (1-varProb.at(varNum));
}

//RANDOM variable choosing heuristic
int Sat::RAND() {

    int firstIndexOfBlock = 0;
    int lastIndexOfBlock = 0;

    //use blockType to check if you're still in the same block of variables
    //-1 corresponds to a chance variable block
    //1 corresponds to a choice variable block
    int blockType = 0;
    bool blockStarted = false;

    int varNum = 0;
    for (vector<int>::iterator iter = assignment.begin();
            iter != assignment.end(); iter++, varNum++) {

        //check if variable has yet to be assigned
        //and if it's the first variable in the block
        if (*iter == 0 && !blockStarted) {
            //first time it's unassigned:
            firstIndexOfBlock = varNum;
            blockStarted = true;
            if (varProb.at(varNum) == 1) {
                blockType = 1;
            } else {
                blockType = -1;
            }
        }

        //check if the block has switched from choice to chance or chance to choice
        //if not, continue on down the block
        if (varProb.at(varNum) != 1 && blockType == 1) {
            lastIndexOfBlock = varNum - 1;
            break;
        } else if (varProb.at(varNum) == 1 && blockType == -1) {
            lastIndexOfBlock = varNum - 1;
            break;
        }
    }
    if (blockType == 0)
        return 0;

    //calculate the last size of the block
    int blockSize = (lastIndexOfBlock - firstIndexOfBlock);

    //pick a random index within the block and return it
    varNum = rand() % ((blockSize + 1) + firstIndexOfBlock);

    //if randNum happens to be an already assigned index,
    //move it until it isn't.
    while (assignment.at(varNum) != 0) {
        if (varNum == lastIndexOfBlock + 1)
            varNum = firstIndexOfBlock;
        else
            varNum++;
    }

//  cout << varNum << endl;
    return varNum;
}

// takes variable assignment and looks in every newly satisfied clause
// decrements numPos and numNeg accordingly for every variable
// IMPORTANT: must call BEFORE we set numActive to zero at every satisfied clauseIndex
void Sat::updateNums(int varAssigned, int sign) {

    if (sign != 1 && sign != 1) {
        return;
    }

    switch (sign) {

    case 1:
        // iterate through posClauses.at(varIndex)
        for (int i = 0; i < posClauses.at(varAssigned).size(); i++) {

            int clauseIndex = posClauses.at(varAssigned).at(i);
            int var;
            int varInactive;

            if (numActive.at(clauseIndex) != 0) {

                // iterate through each clause
                for (int j = 0; j < getClause(clauseIndex)->size(); j++) {

                    var = getClause(clauseIndex)->at(j);
                    varInactive = abs(var) - 1;
                    if (assignment.at(varInactive) != 0) {

                        if (var > 0) {

                            if (numPos.at(varInactive) > 0) {
                                numPos.at(varInactive) --;      // # of positive appearances is reduced by 1
                            }
                        }
                        else if (var < 0) {

                            if (numNeg.at(varInactive) > 0) {
                                numNeg.at(varInactive) --;      // # of negative appearances is reduced by 1
                            }
                        }
                    }
                }
            }
        }
        break;

        case -1:
        // iterate through negClauses.at(varIndex)
        for (int i = 0; i < negClauses.at(varAssigned).size(); i++) {

            int clauseIndex = posClauses.at(varAssigned).at(i);
            int var;
            int varInactive;

            if (numActive.at(clauseIndex) != 0) {

                // iterate through each clause
                for (int j = 0; j < getClause(clauseIndex)->size(); j++) {

                    var = getClause(clauseIndex)->at(j);
                    varInactive = abs(var) - 1;
                    if (assignment.at(varInactive) != 0) {

                        if (var > 0) {

                            if (numPos.at(varInactive) > 0) {
                                numPos.at(varInactive) --;      // # of positive appearances is reduced by 1
                            }
                        }
                        else if (var < 0) {

                            if (numNeg.at(varInactive) > 0) {
                                numNeg.at(varInactive) --;      // # of negative appearances is reduced by 1
                            }
                        }
                    }
                }
            }
        }
    }
}

// HEURISTIC: choose the variable that appears in the most clauses with a single sign
int Sat::MAX_SAT(vector<int> numPos, vector<int> numNeg) {

    // find first unassigned variable
    int blockStart = 0;
    for (int i = 0; i < assignment.size(); i++) {

        if (assignment.at(i
        		) == 0) {
            blockStart = i;
            break;
        }
    }


    bool chance;
    // determine chance/choice
    if (varProb.at(blockStart) != 1) {
        chance = true;
    }
    else{
        chance = false;
    }

    int blockEnd = blockStart;

    // determine the limits of the block
    for (int i = blockStart; i < varProb.size(); i++) {

        if (assignment.at(i) == 0) {

            if (chance == true && varProb.at(i) == 1) {

                chance = false;
                blockEnd = i-1;
                break;
            }
            else if (chance == false && varProb.at(i) != 1) {
                chance = true;
                blockEnd = i - 1;
                break;
            }
        }
    }

    // select var that appears the same way in max number of clauses
    int max = 0;
    int maxVar = 0;
    for (int i = blockStart; i <= blockEnd; i++) {

        if (numPos.at(i) > max) {
            if (assignment.at(i) == 0) {
                max = numPos.at(i);
                maxVar = i;
            }
        }
        if (numNeg.at(i) > max) {
            if (assignment.at(i) == 0) {
                max = numPos.at(i);
                maxVar = i;
            }
        }
    }

    return maxVar;

}
