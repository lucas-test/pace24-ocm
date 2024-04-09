#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <iostream>
#include <algorithm> 
#include <list>

#include "common.h"

using namespace std;

/*
    order: [0,3,1,2]
    pos:   [0,2,3,1]
*/


#include <unordered_set>

/**
 * @brief check for duplicates in a vector.
 * Just used to check if the vectors in the algorithm are ok.
 * 
 * @param numbers 
 * @return true 
 * @return false 
 */
bool has_duplicates(const std::vector<int>& numbers) {
    std::unordered_set<int> seen;
    for (int number : numbers) {
        if (seen.count(number) > 0) {
            return true; // Duplicate found
        }
        seen.insert(number);
    }
    return false; // No duplicates found
}




int lower_bound_pre_order(const vector<vector<int>>& adj, const vector<int>& pre_order, const vector<int>& to_right) {
    int nb_crossings = nb_crossings_from_order(adj, pre_order);
    int n = to_right.size();
    for (int i = 0; i < n; ++i) {
        // crossings with pre_order
        for (int j = 0; j < pre_order.size(); ++j) {
            for (const auto& x : adj[to_right[i]]) {
                for (const auto& y : adj[pre_order[j]]) {
                    if ( x < y ){
                        nb_crossings ++;
                    }
                }
            }
        }
        int max_neighbor = adj[to_right[i]].back(); // adj[to_right[i]][adj[to_right[i]].size()-1];

        // Crossings between to_right[i] and to_right[j]
        for (int j = i + 1; j < n; ++j) {
            if (adj[to_right[j]][0] >= max_neighbor){
                break;
            }
            std::pair<int,int> cr = crossings_between_pair( adj[to_right[i]], adj[to_right[j]]);
            nb_crossings += std::min(cr.first, cr.second);
        }
    }
    return nb_crossings;
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
void aux(const vector<vector<int>>& adj, vector<int>& to_do, vector<int>& order, vector<int>& best_order, int& best_bad_cr, int& current_bad_cr, vector<vector<int>>& triangles, int& triangles_total, const vector<vector<int>>& pair_crossings ){


    if (to_do.size() == 0){
        // End of branch
        if (current_bad_cr < best_bad_cr){
            std::cout << "better " << current_bad_cr << "\n";
            best_order = order;
            best_bad_cr = current_bad_cr;
        }

    } else {

        if (current_bad_cr + triangles_total >= best_bad_cr){
            return;
        }

        // Search a source
        vector<int> is_source(to_do.size(), true);
        for (int i = 0; i < to_do.size(); ++i){
            if (is_source[i] == false) continue;
            int max_neighbor = adj[to_do[i]].back();
            if (triangles[to_do[i]][2] > 0) continue;
            bool no_more_source = false;
            for (int j = 0; j < to_do.size(); ++j){
                if (j < i && adj[to_do[j]].back() < adj[to_do[i]][0]){
                    no_more_source = true;
                    break;
                }
                if ( j == i) continue;
                if (j > i && adj[to_do[j]][0] >= max_neighbor){
                    break;
                }
                // std::pair<int, int> r = crossings_between_pair(adj[to_do[i]], adj[to_do[j]]);
                int rfirst = pair_crossings[to_do[i]][to_do[j]];
                int rsecond = pair_crossings[to_do[j]][to_do[i]];
                if (rfirst < rsecond){
                    is_source[i] = false;
                    break;
                } else if (rfirst > rsecond){
                    is_source[j] = false;
                }
            }

            if (no_more_source) break;

            if (is_source[i]){
                // std::cout << "todo[" << i << "] = " << to_do[i] << " is a source (degree=" << adj[to_do[i]].size() << "\n";
                int x = to_do[i];
                to_do.erase(to_do.begin() + i);
                order.push_back(x);
                aux(adj, to_do, order, best_order, best_bad_cr, current_bad_cr, triangles, triangles_total, pair_crossings);
                to_do.insert(to_do.begin() + i, x);
                order.pop_back();
                return;
            }

        }
        

        


        // Branch on to_do
        // std::cout << "branch " << to_do.size() << "\n";
        for (int i = 0; i < to_do.size(); ++i){
            int x = to_do[i];
           
            // Compute a lower bound on the number of edges that will be crossed if to_do[i] is chosen
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
            if (current_bad_cr + triangles_total + c*mind >= best_bad_cr){ // Because triangles crossings are disjoint from the c crossings
                break;
            }

            int new_current_bad_cr = current_bad_cr;
            for (int j=0; j < to_do.size(); j ++){
                if (adj[x].back() <= adj[to_do[j]][0]){
                    break;
                }
                if (j != i){
                    // std::pair<int, int> r = crossings_between_pair(adj[x], adj[to_do[j]]);
                    int rfirst = pair_crossings[x][to_do[j]];
                    int rsecond = pair_crossings[to_do[j]][x];
                    int diff = -(rfirst - rsecond);
                    if (diff > 0){
                        new_current_bad_cr += diff;
                    }
                }
            } 

            if (new_current_bad_cr + triangles_total >= best_bad_cr){
                continue;
            }


            // Search for a triangle containing x
            int x_triangle_weight = triangles[x][2];
            int new_weight = 0;
            int x2 = -1;


            to_do.erase(to_do.begin() + i);
            order.push_back(x);
            if (x_triangle_weight > 0){
                triangles_total -= x_triangle_weight;
                triangles[x][2] = 0;
                int x0 = triangles[x][0];
                int x1 = triangles[x][1];
                triangles[x0][2] = 0;
                triangles[x1][2] = 0;

                vector<int> new_triangle = find_disjoint_triangle(pair_crossings, to_do, x0 , x1, triangles );
                if (new_triangle[1] > 0){
                    new_weight = new_triangle[1];
                    x2 = new_triangle[0];
                    triangles[x0] = {x1, x2, new_weight};
                    triangles[x1] = {x0, x2, new_weight};
                    triangles[x2] = {x0, x1, new_weight};
                    triangles_total += new_weight;
                    // cout << new_triangle[1] << endl;
                }
            }
            
            aux(adj, to_do, order, best_order, best_bad_cr, new_current_bad_cr, triangles, triangles_total, pair_crossings);

            if (x_triangle_weight > 0){
                triangles_total += x_triangle_weight;
                triangles[x][2] = x_triangle_weight;
                triangles[triangles[x][0]][2] = x_triangle_weight;
                triangles[triangles[x][1]][2] = x_triangle_weight;

                if (new_weight > 0){
                    triangles_total -= new_weight;
                    int x0 = triangles[x][0];
                    int x1 = triangles[x][1];
                    triangles[x0] = {x, x1,  x_triangle_weight};
                    triangles[x1] = {x0, x,  x_triangle_weight};
                    triangles[x2] = {0,0,0};
                }
            }
            to_do.insert(to_do.begin() + i, x);
            order.pop_back();
        }
    }
}

int solver1(const vector<vector<int>>& adj) {
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

    // Sort to_do by increasing leftmost neighbor: adj[i][0]
    // If adj[i][0] == adj[j][0] sort by c(i,j)-c(j,i) (so that it follows the digraph order)
    // If c(i,j) == c(j,i) sort by rightmost neighbor: adj[i].back()
    std::sort(to_do.begin(), to_do.end(), [&adj](int a, int b) {
        if (adj[a][0] != adj[b][0]) {
            return adj[a][0] < adj[b][0];
        } else {
            std::pair<int, int> r = crossings_between_pair(adj[a], adj[b]);
            if (r.first == r.second){
                return adj[a].back() < adj[b].back();
            } else {
                return r.first > r.second;
            }
        }
    });

    int lb = lower_bound(adj);
    std::cout << "lower bound: " << lb << "\n";

    int best_order_nc = nb_crossings_from_order(adj, best_order);
    // std::cout << "greedy nb crossings: " << best_order_nc << "\n";
    int best_bad_cr = best_order_nc - lb;
     std::cout << "greedy bad crossings: " << best_bad_cr << "\n";

    if (best_bad_cr == 0){
        std::cout << "min bad crossings: " << best_bad_cr << "\n";
        return best_bad_cr;
    }

    int current_bad_cr = 0;


    // Triangles
    list<vector<int>> triangles = find_disjoint_3cycles(adj);
    int triangles_total = 0;
    for (const auto& triangle: triangles){
        triangles_total += triangle[3];
    }
    cout << "triangles: " << triangles.size() << " total: " << triangles_total << "\n";

    vector<vector<int>> triangles_adj(adj.size());
    for (int i = 0; i < adj.size(); i ++){
        triangles_adj[i] = {0,0,0};
        for (const auto& triangle: triangles){
            if (triangle[0] == i){
                triangles_adj[i] = {triangle[1], triangle[2], triangle[3]};
                break;
            }
            if (triangle[1] == i){
                triangles_adj[i] = {triangle[0], triangle[2], triangle[3]};
                break;
            }
            if (triangle[2] == i){
                triangles_adj[i] = {triangle[0], triangle[1], triangle[3]};
                break;
            }
        }
    }

    vector<int> order;
    aux(adj, to_do, order, best_order, best_bad_cr, current_bad_cr, triangles_adj, triangles_total, pair_crossings);

    // Print best order found
    // print(best_order);


    std::cout << "min bad crossings: " << best_bad_cr << "\n";
    std::cout << "min crossings: " << lb + best_bad_cr << "\n";


    return best_order_nc;
}