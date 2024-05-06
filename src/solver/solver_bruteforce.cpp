#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <iostream>
#include <algorithm> 
#include <list>

#include "common.h"

using namespace std;

int auxBF(
    vector<int>& to_insert,
    vector<int>& order,
    const vector<vector<int>>& pair_crossings,
    const bool verbose
    ) {
        if (to_insert.size() == 0){
            int nb_cr = nb_crossings_from_order2(order, pair_crossings);
            if (nb_cr <= 1223){
                print(order);
            }
            if (verbose){
                print(order);
                cout << nb_cr << endl;
            }
            return nb_cr;
        } else {

            int min_bad_cr = 100000;

            for (int i = 0; i < to_insert.size(); ++i){
                int v = to_insert[i];
                order.push_back(v);
                to_insert.erase(to_insert.begin() + i);
                int bad_cr = auxBF( to_insert, order, pair_crossings, verbose);
                to_insert.insert(to_insert.begin() + i, v);
                order.pop_back();

                if (bad_cr < min_bad_cr){
                    min_bad_cr = bad_cr;
                }
            }
            return min_bad_cr;
        }
}



int solver_bruteforce(const vector<vector<int>>& adj, const bool& verbose){
    if (verbose){
        cout << "# solver bruteforce" << endl;
        cout << "nb vertices: " << adj.size() << endl;
    }
   


    vector<vector<int>> pair_crossings(adj.size());
    for (int i = 0; i < adj.size(); ++i){
        vector<int> crossings_i(adj.size());
        pair_crossings[i] = crossings_i;
        for (int j = 0; j < adj.size(); ++j){
            pair_crossings[i][j] = crossings_between_pair(adj[i], adj[j]).first;
        }
    }

    vector<int> to_insert;
    for (int i = 0; i < adj.size(); ++i){
        to_insert.push_back(i);
    }

    vector<int> order;

    int min_cr = auxBF(to_insert, order, pair_crossings, false);

    if (verbose){
        print(order);
        cout << "min crossings: " << min_cr << endl;
    }
    return min_cr;
}