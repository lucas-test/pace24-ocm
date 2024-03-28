#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <iostream>
#include <algorithm> 

#include "common.h"

using std::vector;

/*
    order: [0,3,1,2]
    pos:   [0,2,3,1]
*/


int lower_bound_pre_order(const vector<vector<int>>& adj, const vector<int>& pre_order, const vector<int>& to_right) {
    int nb_crossings = nb_crossings_from_order(adj, pre_order);
    // std::cout << "lb pre order " << pre_order.size() << " " << to_right.size() << "\n";
    // print(pre_order);
    // std::cout << nb_crossings << "\n";
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

        // Crossings between to_right[i] and to_right[j]
        for (int j = i + 1; j < n; ++j) {
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
void aux(const vector<vector<int>>& adj, vector<int>& to_do, vector<int>& order, vector<int>& best_order, int& best_order_nc ){


    if (to_do.size() == 0){
        // std::cout << "end " <<"\n";
        // std::cout << "final order ";
        // print(order);
        // std::cout << "best order ";
        // print(best_order);

        int nc = nb_crossings_from_order(adj, order);
        // std::cout << "?? " << order.size() << " " << nc << "\n";

        if (nc < best_order_nc){
            std::cout << "better " << order.size() << " " << nc << "\n";
            // print(order);
            best_order = order; // deep copy
            best_order_nc = nc;
        }
    } else {

        // std::cout << to_do.size() << "\n";

        // std::cout << "lb:" << lower_bound_pre_order(adj, order, to_do) << "\n";

        // Check if there is a left degree 1
        int min_id = adj[to_do[0]][0]; 
        vector<int> left_ids;
        for (int i = 0; i < to_do.size(); i++){
            int x = to_do[i];
            if (adj[x][0] > min_id) break;

            if (adj[x].size() == 1){
                if (adj[x][0] < min_id){
                    min_id = adj[x][0];
                    left_ids.clear(); 
                    left_ids.push_back(i); 
                } else if (adj[x][0] == min_id){
                    left_ids.push_back(i); 
                }
            }
            if (adj[x].size() >= 2){
                if (adj[x][0] < min_id){ // adj[x][0] is the minimal neighbor
                    min_id = adj[x][0];
                    left_ids.clear();
                }
            }
        }

        if (left_ids.size() >= 1){
            std::cout << "left-1:  " << left_ids.size() << "\n";
            vector<int> left_vertices;
            for(int j = left_ids.size() - 1; j >= 0; j--){
                int left_id = left_ids[j];
                int x = to_do[left_id];
                left_vertices.push_back(x);
                to_do.erase(to_do.begin() + left_id);
                order.push_back(x);
            }
            aux(adj, to_do, order, best_order, best_order_nc);

            for( int j = 0 ; j < left_ids.size(); j++){
                int x = left_vertices[j];
                to_do.insert(to_do.begin() + left_ids[j], x);
                order.pop_back();
            }
            return;
        }

        // Search for the leftmost min degree-2
        min_id = adj[to_do[0]][0]; 
        int left2_id = -1;
        int second_neighbor = -1;
        for (int i = 0; i < to_do.size(); i++){
            int x = to_do[i];
            if (adj[x][0] > min_id) break;

            if (adj[x].size() == 2){
                int a = adj[x][0];
                int b = adj[x][1]; // by hypothesis adj[x][1] > adj[x][0]
                if (a < min_id){
                    min_id = a;
                    left2_id = i;
                    second_neighbor = b;
                } else if (a == min_id){
                    if (left2_id == -1 || b < second_neighbor){
                        left2_id = i;
                        second_neighbor = b;
                    }
                }
            }
            else if (adj[x].size() >= 1) {
                int a = adj[x][0]; // min of adj[x]
                int b = adj[x][adj[x].size()-1];

                if (a < min_id){
                    min_id = a;
                    left2_id = -1;
                    second_neighbor = -1;
                } else if (a == min_id && b < second_neighbor){
                    left2_id = -1;
                    second_neighbor = -1;
                }
            }
        }

        if (left2_id >= 0){
            int x = to_do[left2_id];
            // std::cout << "left-2 " << x << " degree: " << adj[to_do[left2_id]].size() << "\n";
            
            to_do.erase(to_do.begin() + left2_id);
            order.push_back(x);
            aux(adj, to_do, order, best_order, best_order_nc);
            to_do.insert(to_do.begin() + left2_id, x);
            order.pop_back();
            return;
        }


        // Check if to_do[0] is a source
        for (int i = 0; i < 5 && i < to_do.size(); i ++){
            int max_neighbor = adj[to_do[i]][adj[to_do[i]].size()-1];
            bool is_source = true;
            for (int j = 0; j < to_do.size(); j++){
                if ( j == i) continue;
                if (adj[to_do[j]][0] >= max_neighbor){
                    break;
                }
                std::pair<int, int> r = crossings_between_pair(adj[i], adj[j]);
                if (r.first > r.second){
                    is_source = false;
                    break;
                }
            }

            if (is_source){
                // std::cout << "todo[" << i << "] is a source\n";
                int x = to_do[i];
                to_do.erase(to_do.begin() + i);
                order.push_back(x);
                aux(adj, to_do, order, best_order, best_order_nc);
                to_do.insert(to_do.begin() + i, x);
                order.pop_back();
                return;
            }
        }
        

        int lb = lower_bound_pre_order(adj, order, to_do);
        // std::cout << "lower bound pre order " << lb << " " << order.size() << "\n";
        if (lb >= best_order_nc){
            // std::cout << "cut " << to_do.size() << " " << lb << " " << best_order_nc << "\n";
            return;
        }

        // Branch on to_do
        // std::cout << "branch " << to_do.size() << "\n";
        for (int i = 0; i < to_do.size(); i ++){
            int x = to_do[i];
            to_do.erase(to_do.begin() + i);
            order.push_back(x);
            aux(adj, to_do, order, best_order, best_order_nc);
            to_do.insert(to_do.begin() + i, x);
            order.pop_back();
        }
    }
}

int solver1(const vector<vector<int>>& adj) {
    std::cout << "loal\n";
    vector<int> pos = greedy_sequential(adj);
    vector<int> best_order(pos.size());
    for (int i = 0; i < pos.size(); i++) {
        best_order[pos[i]] = i ;
    }

    vector<int> to_do;
    for (int i = 0; i < adj.size(); i ++){
        to_do.push_back(i);
    }

    // Sort to_do by increasing adj[i][0]
    std::sort(to_do.begin(), to_do.end(), [&adj](int a, int b) {
        return adj[a][0] < adj[b][0];
    });


    // std::cout << nb_crossings(adj, pos) << "\n";
    // std::cout << nb_crossings_from_order(adj, best_order) << "\n";

    int best_order_nc = nb_crossings_from_order(adj, best_order);
    vector<int> order;

    // std::cout << "best order ";
    // for (const auto& elem : best_order) {
        // std::cout << elem << ' ';
    // }
    // std::cout << std::endl;
    std::cout << best_order_nc << "\n";


    aux(adj, to_do, order, best_order, best_order_nc  );

    // Print best order found
    // print_adj(adj);
    // std::cout << "best order ";
    // for (const auto& elem : best_order) {
    //     std::cout << elem << ' ';
    // }
    // std::cout << std::endl;

    std::cout << best_order_nc << "\n";

    return best_order_nc;
}