
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <iostream>
#include <algorithm> 
#include <list>

#include "common.h"

using namespace std;

// TODO use pair_crossings
vector<int> order_greedy_sequential_mask(const vector<vector<int>>& adj, const vector<bool>& mask) {
    int n = adj.size();
    vector<int> order;
    for (int i = 0; i < n; i++) {

        if (mask[i] == false) continue;
        int cr = 0;
        int minCr = 0;
        int minK = 0;
        for (int k = 0; k < order.size(); k++) {
            int j = order[k];
            auto [iLeft, iRight] = crossings_between_pair(adj[j], adj[i]);
            cr += iRight - iLeft;
            if (cr < minCr) {
                minCr = cr;
                minK = k + 1;
            }
        }
        order.insert(order.begin() + minK, i);
    }
    return order;
}




vector<int> order_greedy_sequential_mask2(const vector<vector<int>>& adj, const vector<bool>& mask, const vector<vector<int>>& pair_crossings, int& nb_bad_cr) {
    int n = adj.size();
    nb_bad_cr = 0;
    vector<int> order;
    for (int i = 0; i < n; i++) {
        if (mask[i] == false) continue;
        int cr = 0;
        int minCr = 0;
        int minK = 0;
        int if_leftmost_cr = 0;
        for (int k = 0; k < order.size(); k++) {
            int j = order[k];
            int iLeft = pair_crossings[j][i];
            int iRight = pair_crossings[i][j];
            if_leftmost_cr += iLeft - min(iLeft, iRight);
            // auto [iLeft, iRight] = crossings_between_pair(adj[j], adj[i]);
            cr += iRight - iLeft;
            if (cr < minCr) {
                minCr = cr;
                minK = k + 1;
            }
        }
        nb_bad_cr += if_leftmost_cr + minCr;
        order.insert(order.begin() + minK, i);
    }

    return order;
}















/**
 * @brief 
 * 
 * @param adj i -> adj[i][0] is increasing, j -> adj[i][j] is increasing
 * @param to_do sorted by adj[todo[i]][0]
 * @param order 
 * @param best_order 
 * @param best_order_nc 
 */
void aux(const vector<vector<int>>& adj, vector<int>& to_do, vector<int>& order, vector<int>& best_order, int& best_bad_cr, int& current_bad_cr, const vector<vector<int>>& pair_crossings, const pair<vector<vector<int>>, const vector<vector<int>>> digraph, vector<bool>& mask, int deep ){


    if (to_do.size() == 0){
        // End of branch
        if (current_bad_cr < best_bad_cr){
            cout << "better " << current_bad_cr << "\n";
            best_order = order;
            best_bad_cr = current_bad_cr;
        }

    } else {
        vector<vector<int>> in_neighbors = digraph.first;

        int min_indegree = 100000;
        for (int j = 0; j < to_do.size(); ++j){
            int indegree = 0;
            for (int k = 0; k < in_neighbors[to_do[j]].size(); ++k){
                if (mask[in_neighbors[to_do[j]][k]] ) indegree ++;
            }
            if (indegree < min_indegree){
                min_indegree = indegree;
            }
        }

        if (current_bad_cr + min_indegree >= best_bad_cr){
            return;
        }


        vector<vector<int>> compo = scc_mask(digraph.second, digraph.first, mask); // restreint à mask

        if (compo.size() >= 2){

            int c3 = 0;
            // cout << compo.size() << endl;
            for (int i = 0; i < compo.size(); ++i){
                if (compo[i].size() >= 3){
                    
                    // vector<bool> rmask(adj.size(), false);
                    // for (int j =0; j < compo[i].size(); ++j){
                    //     rmask[compo[i][j]] = true;
                    // }

                    // min_indegree = 100000;
                    // for (int j = 0; j < compo[i].size(); ++j){
                    //     int indegree = 0;
                    //     for (int k = 0; k < in_neighbors[compo[i][j]].size(); ++k){
                    //         if (rmask[in_neighbors[compo[i][j]][k]] ) indegree ++;
                    //     }
                    //     if (indegree < min_indegree){
                    //         min_indegree = indegree;
                    //     }
                    // }
                    // cout << compo[i].size() << " " << min_indegree << endl;
                    // c3 += min_indegree;
                    c3 ++;
                }
            }
            // if (c3 >= 1) cout << c3 << endl;
            if (current_bad_cr + c3 >= best_bad_cr) return;

            vector<int> sub_order;
            for (int i = 0; i < compo.size(); ++i){
                if (deep == 0 && compo[i].size() >= 2 ) cout << i << "/" << compo.size() << " size " << compo[i].size() <<  endl;
                if (compo[i].size() == 1) {
                    sub_order.push_back(compo[i][0]);
                    continue;
                }

                // tout recalculer

                sort(compo[i].begin(), compo[i].end(), [&adj](int a, int b) {
                    return adj[a].size() > adj[b].size();
                });


                vector<bool> rmask(adj.size(), false);
                for (int j =0; j < compo[i].size(); ++j){
                    rmask[compo[i][j]] = true;
                }

                vector<int> rorder;


                // V0
                vector<int> rbest_order;
                for (int i = 0; i < best_order.size(); ++i){
                    if (rmask[best_order[i]]) rbest_order.push_back(best_order[i]);
                }
                int rbest_order_nc = nb_crossings_from_order2(adj, rbest_order, pair_crossings); // compute in greedy algo
                int lb = lower_bound_mask(adj, pair_crossings, compo[i]);
                int rbest_bad_cr = rbest_order_nc - lb;

                // V1
                // vector<int> rbest_order = order_greedy_sequential_mask3(adj, compo[i], pair_crossings);
                // int rbest_order_nc = nb_crossings_from_order2(adj, rbest_order, pair_crossings); // compute in greedy algo
                // int lb = lower_bound_mask(adj, pair_crossings, compo[i]);
                // int rbest_bad_cr = rbest_order_nc - lb;

                // V2
                // int rbest_bad_cr = 0;
                // auto lol = order_greedy_sequential_mask2(adj, rmask, pair_crossings, rbest_bad_cr);
                // vector<int> rbest_order = lol;
                // int rbest_bad_cr = lol.second ;
                // if (rbest_order_nc != lol.second) {
                //     cout << "pb" << rbest_order_nc << " " << lol.second << endl;
                //     print(rbest_order);
                //     print(lol.first);
                
                // }
                // if (rbest_bad_cr2 != rbest_bad_cr) cout << "pb" << rbest_bad_cr2 << " " << rbest_bad_cr << endl;

                sort(compo[i].begin(), compo[i].end(), [&adj](int a, int b) {
                    return adj[a][0] < adj[b][0];
                });

                int rcurrent_bad_cr = 0;

                if (rbest_bad_cr > 0){
                    aux(adj, compo[i], rorder, rbest_order, rbest_bad_cr, rcurrent_bad_cr, pair_crossings, digraph, rmask , deep+1);
                }


                // Concatenate rbest_order
                current_bad_cr += rbest_bad_cr;
                // best_order.push_back(rbest_order);
                sub_order.insert(sub_order.end(), rbest_order.begin(), rbest_order.end());

                if (current_bad_cr >= best_bad_cr){
                    return;
                }


            }

            if (best_bad_cr > current_bad_cr){
                best_order = order;
                best_order.insert(best_order.end(), sub_order.begin(), sub_order.end());
                best_bad_cr = current_bad_cr;
            }

            return;
        } else {
            // cout << compo[0].size() << endl;
        }

        
        int triangles_w = find_disjoint_triangles(adj, to_do, pair_crossings);

        if (current_bad_cr + triangles_w >= best_bad_cr){
            return;
        }


        // Branch on to_do
        // cout << "branch " << to_do.size() << "\n";
        for (int i = 0; i < to_do.size(); ++i){
            int x = to_do[i];


            // Cut
            int c = 0;
            for (int j = 0; j < i ; ++j){
                if (adj[to_do[j]].back() < adj[x][0]){
                    c += adj[to_do[j]].size();
                } else if (adj[to_do[j]].back() == adj[x][0] ) {
                    c += adj[to_do[j]].size()-1;
                } 
            }
            int mind = adj[x].size();
            for (int j = i+1; j < to_do.size(); ++j){
                if (adj[to_do[j]].size() < mind){
                    mind = adj[to_do[j]].size();
                }
            }
            if (current_bad_cr + triangles_w + c*mind >= best_bad_cr){ // Because triangles crossings are disjoint from the c crossings
                break;
            }


           
            int new_current_bad_cr = current_bad_cr;
            for (int j=0; j < to_do.size(); j ++){
                if (adj[x].back() <= adj[to_do[j]][0]){
                    break;
                }
                if (j != i){
                    // pair<int, int> r = crossings_between_pair(adj[x], adj[to_do[j]]);
                    int rfirst = pair_crossings[x][to_do[j]];
                    int rsecond = pair_crossings[to_do[j]][x];
                    int diff = -(rfirst - rsecond);
                    if (diff > 0){
                        new_current_bad_cr += diff;
                    }
                }
            } 

            if (new_current_bad_cr  >= best_bad_cr){
                continue;
            }


            to_do.erase(to_do.begin() + i);
            order.push_back(x);
            mask[x] = false;
            aux(adj, to_do, order, best_order, best_bad_cr, new_current_bad_cr, pair_crossings, digraph, mask, deep+1);
            mask[x] = true;
            to_do.insert(to_do.begin() + i, x);
            order.pop_back();
        }
    }
}

int solver1a(const vector<vector<int>>& adj) {
    cout << "nb vertices: " << adj.size() << endl;

    vector<int> pos = greedy_sequential(adj);
    vector<int> best_order(pos.size());
    for (int i = 0; i < pos.size(); ++i) {
        best_order[pos[i]] = i ;
    }

    vector<int> to_do;
    for (int i = 0; i < adj.size(); ++i){
        to_do.push_back(i);
    }

    vector<vector<int>> pair_crossings(adj.size());
    for (int i = 0; i < adj.size(); ++i){
        vector<int> crossings_i(adj.size());
        pair_crossings[i] = crossings_i;
        for (int j = 0; j < adj.size(); ++j){
            pair_crossings[i][j] = crossings_between_pair(adj[i], adj[j]).first;
        }
    }

    

    auto digraph = compute_directed_graph(adj);

    vector<bool> mask(adj.size(), true);

    // print(best_order);
    // auto lol = order_greedy_sequential_mask2(adj, mask, pair_crossings);
    // print(lol.first);
    // cout << "lol ?" << lol.second << endl;

    // Sort to_do by increasing leftmost neighbor: adj[i][0]
    // If adj[i][0] == adj[j][0] sort by c(i,j)-c(j,i) (so that it follows the digraph order)
    // If c(i,j) == c(j,i) sort by rightmost neighbor: adj[i].back()
    sort(to_do.begin(), to_do.end(), [&adj](int a, int b) {
        if (adj[a][0] != adj[b][0]) {
            return adj[a][0] < adj[b][0];
        } else {
            pair<int, int> r = crossings_between_pair(adj[a], adj[b]);
            if (r.first == r.second){
                return adj[a].back() < adj[b].back();
            } else {
                return r.first > r.second;
            }
        }
    });

    int lb = lower_bound(adj);
    cout << "lower bound: " << lb << "\n";

    int best_order_nc = nb_crossings_from_order(adj, best_order);
    // cout << "greedy nb crossings: " << best_order_nc << "\n";
    int best_bad_cr = best_order_nc - lb;
     cout << "greedy bad crossings: " << best_bad_cr << "\n";

    if (best_bad_cr == 0){
        cout << "min bad crossings: " << best_bad_cr << "\n";
        return best_bad_cr;
    }

    int current_bad_cr = 0;

    vector<int> order;
    aux(adj, to_do, order, best_order, best_bad_cr, current_bad_cr, pair_crossings, digraph, mask, 0);

    // Print best order found
    // print(best_order);


    cout << "min bad crossings: " << best_bad_cr << "\n";
    cout << "min crossings: " << lb + best_bad_cr << "\n";

    // Check solution
    // print(best_order);
    // cout << best_order.size() << endl;
    // cout << adj.size() << endl;
    // cout << has_duplicates(best_order) << endl;
    // cout << nb_crossings_from_order(adj, best_order) << endl;

    return best_order_nc;
}