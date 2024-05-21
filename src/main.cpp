#include <iostream>
#include <vector>
#include <string>
#include <iomanip> 
#include <algorithm>
#include <random>
#include <chrono>
#include <unordered_set>

#include "solver/common.h"
#include "solver/solver1.h"
#include "solver/solver1b.h"
#include "solver/solver2.h"
#include "solver/solver3.h"
#include "solver/solver4.h"
#include "solver/solver_bruteforce.h"

using namespace std;


void search_random(void);
void compare_greedy_insertion_orders(const vector<vector<int>>& adj);
pair<vector<vector<int>>, int> load_stdin() ;


int main(int argc, char* argv[]) {

    pair<vector<vector<int>>, int> r = load_stdin();
    vector<vector<int>> adj = r.first;
    int num1 = r.second;

    vector<int> final_order =  solver4(adj, false);

    for (int i = 0; i < final_order.size(); ++i){
        cout << final_order[i] + num1 + 1 << endl;
    }

    return 0;
}


void compare_greedy_insertion_orders(const vector<vector<int>>& adj){

    // // Natural
    // int greedy_nat = nb_crossings(adj, greedy_sequential(adj));
    // // cout << "greedy: " << greedy_nat << "\n";

    // // Sort by increasing degree
    // sort(adj.begin(), adj.end(), [](const vector<int>& a, const vector<int>& b) {
    //     return a.size() < b.size();
    // });
    // int greedy_inc_deg = nb_crossings(adj, greedy_sequential(adj));
    // // cout << "greedy increasing degree: " << greedy_inc_deg << "\n";

    // // Sort by decreasing degree
    // sort(adj.begin(), adj.end(), [](const vector<int>& a, const vector<int>& b) {
    //     return a.size() > b.size();
    // });
    // int greedy_dec_deg = nb_crossings(adj, greedy_sequential(adj));
    // // cout << "greedy decreasing degree: " << greedy_dec_deg << "\n";

    // // Sort by increasing max neighbor
    // sort(adj.begin(), adj.end(), [](const vector<int>& a, const vector<int>& b) {
    //    if (a.empty() || b.empty()) {
    //         return a.size() < b.size();
    //     } else {
    //         return *max_element(a.begin(), a.end()) < *max_element(b.begin(), b.end());
    //     }
    // });
    // vector<int> posGIMX = greedy_sequential(adj);
    // int greedy_inc_max_neighbor = nb_crossings(adj, greedy_sequential(adj));
    // // cout << "greedy increasing max neighbor: " << greedy_inc_max_neighbor  << "\n";
    
    // // vector<int> orderGIMX(posGIMX.size());
    // // for (int i = 0; i < posGIMX.size(); i++) {
    // //     orderGIMX[posGIMX[i]] = i ;
    // // }
    // // cout << nb_crossings_from_order(adj, orderGIMX) << "\n";
    // // print(orderGIMX);

    //  // Sort by decreasing max neighbor
    // sort(adj.begin(), adj.end(), [](const vector<int>& a, const vector<int>& b) {
    //    if (a.empty() || b.empty()) {
    //         return a.size() < b.size();
    //     } else {
    //         return *max_element(a.begin(), a.end()) > *max_element(b.begin(), b.end());
    //     }
    // });
    // int greedy_dec_max_neighbor = nb_crossings(adj, greedy_sequential(adj));
    // // cout <<  "greedy decreasing max neighbor: " << greedy_dec_max_neighbor << "\n";
   
    // cout << lower_bound(adj) << " " << greedy_nat << " " << greedy_inc_deg << " " << greedy_dec_deg << " " << greedy_inc_max_neighbor << " " << greedy_dec_max_neighbor << "\n";

}


/**
 * @brief for finding a complicated graphs with not that much vertices: it should be done in some seconds by solver1
 * 
 */
void search_random(void){
    srand(static_cast<unsigned>(time(0)));

    for (int i = 0; i < 100000 ; i ++){

        double randomValue = static_cast<double>(rand()) / static_cast<double>(RAND_MAX);

        vector<vector<int>> adj = generate_random_adj(11, 13 + (i % 30), 0.2 + ((double) (i % 10))/10.);
        reduce_degree_0(adj);
        vector<int> pos = greedy_sequential(adj);
        if (lower_bound(adj) != nb_crossings(adj, pos)){
            // cout << "try " << i << endl;

            auto start = chrono::high_resolution_clock::now();
            int min_cr = solver1b(adj, false);
            auto end = chrono::high_resolution_clock::now();
            chrono::duration<double> elapsed = end - start;

            
            // cout << elapsed.count() << endl;
            if (elapsed.count() > 0.00002){
                int min_cr_brute = solver_bruteforce(adj, false);
                if (min_cr_brute != min_cr){
                    cout << "------" << endl;
                    cout << "bug: brute=" << min_cr_brute << " solver=" << min_cr << endl;
                    cout << "p=" << randomValue << endl;
                    print_adj(adj);
                    print_gr_format(adj);
                    cout << "Lower bound: " << lower_bound(adj) << "\n";
                    cout << "Greedy crossings: " << nb_crossings(adj, pos) << endl;
                    cout << "Execution time of solver1: " << elapsed.count() << " seconds\n";
                } else {
                    // cout << "ok " << lower_bound(adj) << " " << min_cr << " " << nb_crossings(adj, pos) << endl;
                }
               
            }

           

            
        }
    }
}


pair<vector<vector<int>>, int> load_stdin() {
    string line;
    vector<vector<int>> adj;
    int num1, num2, num3;
    string type, ocr;

    getline(cin, line);
    while (!line.empty() && line[0] == 'c') {
        getline(cin, line);
    }

    // First line contains num1, num2, num3
    istringstream lineStream(line);
    lineStream >> type >> ocr >> num1 >> num2 >> num3;
    adj.resize(num2, vector<int>());

    // Read the rest of the input
    while (getline(cin, line)) {
        if (!line.empty() && line[0] == 'c') continue;

        istringstream lineStream(line);
        int firstNumber, secondNumber;
        if (lineStream >> firstNumber >> secondNumber) {
            auto it = upper_bound(adj[secondNumber - num1 - 1].begin(), adj[secondNumber - num1 - 1].end(), firstNumber - 1);
            adj[secondNumber - num1 - 1].insert(it, firstNumber - 1);
        }
    }

    return make_pair(adj, num1);
}