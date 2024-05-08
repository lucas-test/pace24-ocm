#ifndef SOLVER_BRUTEFORCE_H
#define SOLVER_BRUTEFORCE_H

#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <iostream>
#include <algorithm> 

using namespace std;

int solver_bruteforce(const vector<vector<int>>& adj, const bool& verbose);
int auxBF(
    vector<int>& to_insert,
    vector<int>& order,
    const vector<vector<int>>& pair_crossings,
    const bool verbose
    );


#endif 