#include <iostream>
#include <vector>
#include <string>
#include <iomanip> 
#include <algorithm>
#include <random>
#include <chrono>
#include <unordered_set>

#include "common.h"
#include "solver1.h"
#include "solver1a.h"
#include "solver1b.h"
#include "solver_bruteforce.h"

using namespace std;

void search_random(void);
void compare_greedy_insertion_orders(const vector<vector<int>>& adj);


int main(int argc, char* argv[]) {
    
    // search_random();


    if (argc != 2) {
        cerr << "Usage: " << argv[0] << " <filename>\n";
        return 1;
    }

    cout << "\n" << argv[1] << endl;


    vector<vector<int>> adj = load_file(argv[1]);

    reduce_degree_0(adj);
    // cout << adj.size() << "\n";


    // Solver1 test

    // cout << "lol: " << nb_crossings(adj, greedy_sequential(adj)) << "\n";

    // sort(adj.begin(), adj.end(), [](const vector<int>& a, const vector<int>& b) {
    //     return a.size() < b.size();
    // });

    //     cout << "lol: " << nb_crossings(adj, greedy_sequential(adj)) << "\n";


    // sort(adj.begin(), adj.end(), [](const vector<int>& a, const vector<int>& b) {
    //     return a.size() > b.size();
    // });

    //     cout << "lol: " << nb_crossings(adj, greedy_sequential(adj)) << "\n";


    // sort(adj.begin(), adj.end(), [](const vector<int>& a, const vector<int>& b) {
    //    if (a.empty() || b.empty()) {
    //         return a.size() < b.size();
    //     } else {
    //         return *max_element(a.begin(), a.end()) > *max_element(b.begin(), b.end());
    //     }
    // });
    // print_adj(adj);

    


    // Try to find a better greedy insertion order by starting with the triangles
    // Not so much better ... or not better
    // int lb = lower_bound(adj);
    // cout<< "lb: " << lb << "\n";
    // cout << "greedy lol: " << nb_crossings(adj, greedy_sequential(adj))-lb << "\n";

    // list<vector<int>> triangles = find_disjoint_3cycles(adj);
    // cout << "nb triangles: " << triangles.size() << "\n";

    // vector<vector<int>> adj2;
    // vector<int> seen(adj.size(), false);
    // for (const auto& triangle : triangles){
    //     adj2.push_back(adj[triangle[0]]);
    //     adj2.push_back(adj[triangle[1]]);
    //     adj2.push_back(adj[triangle[2]]);
    //     seen[triangle[0]] = true;
    //     seen[triangle[1]] = true;
    //     seen[triangle[2]] = true;
    // }
    // for (int i = 0; i < seen.size(); ++i){
    //     if (seen[i] == false){
    //         adj2.push_back(adj[i]);
    //     }
    // }

    // int greedy_inc_deg = nb_crossings(adj2, greedy_sequential(adj2));
    // cout << "greedy lol: " << greedy_inc_deg-lb << "\n";


    // auto digraph = compute_directed_graph(adj);
    // auto compo = scc(digraph.second, digraph.first);

    // print_adj(adj);

    // for (int i = 0; i < digraph.first.size(); ++i){
    //     if (digraph.first[i].size() == 0){
    //         cout << i << " is source" << endl;
    //     }
    //     // cout << i << " <-:";
    //     // for (int j = 0; j < digraph.first[i].size(); ++j){
    //     //     cout << digraph.first[i][j] << " ";
    //     // }
    //     // cout << endl;
    // }

    // for (int i = 0; i < digraph.second.size(); ++i){
    //     cout << i << " ->:";
    //     for (int j = 0; j < digraph.second[i].size(); ++j){
    //         cout << digraph.second[i][j] << " ";
    //     }
    //     cout << endl;
    // }

    // for (int i = 0; i < compo.size(); ++i){
    //     cout << i << ":";
    //     for (int j = 0; j < compo[i].size(); ++j){
    //         cout << compo[i][j] << " ";
    //     }
    //     cout << endl;
    // }

    // print_adj(adj);
    // vector<int> pos = greedy_sequential(adj);
    // cout << lower_bound(adj)  << " " << nb_crossings(adj, pos) << endl;

    // cout << solver_bruteforce(adj, false) << endl;
    // cout << endl;
    cout << solver1b(adj, true) << endl;

    // search_random();

    // vector<vector<int>> pair_crossings(adj.size());
    // for (int i = 0; i < adj.size(); ++i){
    //     vector<int> crossings_i(adj.size());
    //     pair_crossings[i] = crossings_i;
    //     for (int j = 0; j < adj.size(); ++j){
    //         pair_crossings[i][j] = crossings_between_pair(adj[i], adj[j]).first;
    //     }
    // }

    // auto digraph = compute_directed_graph(adj);
    // vector<int> vertices;
    // for (int i = 0; i < adj.size(); ++i){
    //     vertices.push_back(i);
    // }

    // unordered_set<int> subvertices = {174,82,246,289,236,302,129,55,204,384,299,83,316,34,374,17,109,312,355,228,64,381,106,376,11,372,232};

    // for (int i = 0; i < adj.size(); ++i){
    //     if (subvertices.find(i) != subvertices.end()){
    //         cout << i << " <- ";
    //         for (int x: digraph.first[i]){
    //             if (subvertices.find(x) != subvertices.end()){
    //                 cout << x << " ";
    //             }
    //         }
    //         cout << endl;
    //         cout << i << " -> ";
    //         for (int x: digraph.second[i]){
    //             if (subvertices.find(x) != subvertices.end()){
    //                 cout << x << " ";
    //             }
    //         }
    //         cout << endl;
    //     }
    // }

    // find_triangles_tree(pair_crossings, digraph.first, digraph.second, vertices);

    // find_edge_disjoint_triangles(pair_crossings, digraph.first, digraph.second, vertices);

    // list<vector<int>> triangles = find_disjoint_3cycles(adj);
    // int triangles_total = 0;
    // for (const auto& triangle: triangles){
    //     triangles_total += triangle[3];
    // }
    // cout << "triangles: " << triangles.size() << " total: " << triangles_total << "\n";

    // // Sort by increasing degree
    // sort(adj.begin(), adj.end(), [](const vector<int>& a, const vector<int>& b) {
    //     return a.size() < b.size();
    // });
    // int greedy_inc_deg = nb_crossings(adj, greedy_sequential(adj));
    // cout << "greedy increasing degree: " << greedy_inc_deg << "\n";

    // // Sort by decreasing degree
    // sort(adj.begin(), adj.end(), [](const vector<int>& a, const vector<int>& b) {
    //     return a.size() > b.size();
    // });
    // vector<int> greedy_dec_deg_pos = greedy_sequential(adj);
    // int greedy_dec_deg = nb_crossings(adj, greedy_sequential(adj));
    // cout << "greedy decreasing degree: " << greedy_dec_deg << "\n";
    
    // vector<int> greedy_dec_deg_order(greedy_dec_deg_pos.size());
    // for (int i = 0; i < greedy_dec_deg_pos.size(); i++) {
    //     greedy_dec_deg_order[greedy_dec_deg_pos[i]] = i ;
    // }
    // int greedy_dec_deg_order_nc = nb_crossings_from_order(adj, greedy_dec_deg_order);
    // cout << "greedy decreasing degree order: " << greedy_dec_deg_order_nc << "\n";


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