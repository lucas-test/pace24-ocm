
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <iostream>
#include <algorithm> 
#include <list>
#include <chrono>

#include "common.h"

using namespace std;






/**
 * @brief 
 * 
 * @param adj 
 * @param to_do 
 * @param order 
 * @param best_order 
 * @param best_order_nc 
 */
void aux3(const vector<vector<int>>& adj, 
    const vector<vector<int>>& twins,
    vector<int>& to_do, 
    vector<int>& order, 
    vector<int>& best_order, 
    int& best_bad_cr, 
    int& current_bad_cr, 
    const vector<vector<int>>& pair_crossings, 
    vector<vector<int>>& in_neighbors,
    vector<vector<int>>& out_neighbors, 
    vector<bool>& mask, 
    int depth, 
    vector<vector<int>>& triangles_adj, 
    int triangles_total,
    vector<bool>& excluded,
    const vector<vector<vector<int>>>& triangles ){


    if (to_do.size() == 0){
        if (current_bad_cr < best_bad_cr){
            cout << "better " << current_bad_cr << "\n";
            best_order = order;
            best_bad_cr = current_bad_cr;
        }
    } else {

        if (best_bad_cr - current_bad_cr <= 10){
            vector<int> indegrees;
            for (const int& v: to_do){
                indegrees.push_back(in_neighbors[v].size());
            }
            sort(indegrees.begin(), indegrees.end());
            int rind = 0;
            for (int i = 0; i < to_do.size(); ++i){
                if (indegrees[i] >= i){
                    rind += indegrees[i] - i;
                }
            }
            if (current_bad_cr + rind >= best_bad_cr) return;
        }
        

        vector<int> in_weights(adj.size());
        for( const int& x: to_do){
            for (const int& v: in_neighbors[x]){
                in_weights[x] += pair_crossings[v][x] - pair_crossings[x][v];
            }
        }

        sort(to_do.begin(), to_do.end(), [&in_weights](int a, int b) {
            if (in_weights[a] == in_weights[b]){
                return a > b;
            }
            return in_weights[a] < in_weights[b];
        });

        // cout << string(depth, '-');
        // print(to_do);

        if (in_weights[to_do[0]] == 0){

            vector<int> sources;

            for (int i = 0; i < to_do.size(); ++i){
                int x = to_do[i];
                if (in_weights[x] > 0){
                    break;
                }
                sources.push_back(x);
            }


            

            vector<vector<int>> old_in_neighbors(adj.size());
            for (const int& y: to_do){
                if (in_weights[y] > 0)
                    old_in_neighbors[y] = in_neighbors[y];
            }
            
            vector<vector<int>> old_out_neighbors(adj.size());
            for (const int& y: sources){
                old_out_neighbors[y] = out_neighbors[y];
            }

            for (const int& v: to_do){
                if (in_weights[v] > 0){
                    in_neighbors[v] = vector<int>();
                    for (const int& w: old_in_neighbors[v]){
                        if (in_weights[w] > 0){
                            in_neighbors[v].push_back(w);
                        }
                    }
                }
            }

            vector<int> new_to_do;
            for (size_t j = 0; j < to_do.size(); ++j){
                if (in_weights[to_do[j]] > 0) 
                    new_to_do.push_back(to_do[j]);
            }

            for (const int& y: sources){
                order.push_back(y);
                mask[y] = false;
                out_neighbors[y] = {};
            }

            vector<bool> next_excluded(adj.size(), false);

            aux3(adj, twins, new_to_do, order, best_order, best_bad_cr, current_bad_cr, pair_crossings, in_neighbors, out_neighbors, mask, depth+1, triangles_adj, 0, next_excluded, triangles);


            for (const int& y: sources){
                order.pop_back();
                mask[y] = true;
                out_neighbors[y] = old_out_neighbors[y];
            }


            for (const int& v: to_do){
                if (in_weights[v] > 0){
                    in_neighbors[v] = old_in_neighbors[v];
                }
            }
            


            return;
        }
       

        // Branch on to_do
        // cout << "branch " << to_do.size() << "\n";
        for (int i = 0; i < to_do.size(); ++i){
            int x = to_do[i];
            if (excluded[x]){
                // cout << string(depth, '-') << "branch exclude " << x << endl;
                continue;
            } 

            

            if (order.size() >= 2){
                int y = order[order.size()-1];
                int z = order[order.size()-2];
                int xz = pair_crossings[x][z] - pair_crossings[z][x];
                int w1 = pair_crossings[z][y] - pair_crossings[y][z];
                int w2 = pair_crossings[y][x] - pair_crossings[x][y];
                int mw = min(w1, w2);
                if (xz > mw){
                    // cout << string(depth, '-') << "branch cut non optimal last third " << x << endl;
                    continue;
                }

                // if ( xz == mw && xz == w1 && x < z){
                //     continue;
                // }
                // if (xz == mw && xz == w2 && x < y){
                //     continue;
                // }
            }

            if (order.size() >= 3){
                int a = order[order.size()-3];
                int b = order[order.size()-2];
                int c = order[order.size()-1];
                int bc = pair_crossings[b][c] - pair_crossings[c][b];
                int ca = pair_crossings[c][a] - pair_crossings[a][c];
                // int xa = pair_crossings[x][a] - pair_crossings[a][x];
                int xb = pair_crossings[x][b] - pair_crossings[b][x];

                int xa = pair_crossings[x][a] - pair_crossings[a][x];
                if (xa < 0) xa = 0;
                
                if ( ca >= 0 && xb >= 0  && bc < ca + (xb + xa)*twins[x].size() ){
                    // cout << string(depth, '-') << "branch cut non optimal last third " << x << endl;
                    // c x a b or c a x b is better than a b c x and has only one BA: bc
                    // we are sure that there are arcs ab bc and cd
                    continue;
                }

                int cx = pair_crossings[c][x] - pair_crossings[x][c];
                if ( ca <= 0 && xb >= 0  && cx <  xb + xa ){
                    // cout << string(depth, '-') << "branch cut non optimal last third " << x << endl;
                    // c x a b or c a x b is better than a b c x and has only one BA: bc
                    // we are sure that there are arcs ab bc and cd
                    continue;
                }
            }

            // Twins should be not be reconsidered later in the current iteration
            for (const int& y: twins[x]){
                excluded[y] = true;
            }


            int s = 0;
            bool optimal = true;
            for (int j = order.size()-1;  j >= 0; --j){
                int y = order[j];
                int xy = pair_crossings[x][y] - pair_crossings[y][x];
                if (xy > 0){
                    s -= xy;
                } else if (xy < 0){
                    s += -xy;
                }
                if (s < 0) {
                    optimal = false;
                    break;
                }
            }
            if (optimal == false) continue;

            

            int x_bad_cr = in_weights[x];


            // Only if to_do is sorted by increasing in_weights
            if (current_bad_cr + x_bad_cr >= best_bad_cr){
                // cout << string(depth, '-') << "branch mega cut " << x << endl;
                break; 
            } 
            

            int new_current_bad_cr = current_bad_cr + x_bad_cr*(twins[x].size());




            if (new_current_bad_cr  >= best_bad_cr){
                // cout << string(depth, '-') << "branch cut " << x << endl;
                continue;
            }

            
            // cout << string(depth, '-') << "branch " << x << " in neighbors: ";
            // print(in_neighbors[x]);

            vector<bool> next_excluded(adj.size(), false);
            for (const int& in_neighbor: in_neighbors[x]){
                next_excluded[in_neighbor] = true;
            }
            for (const int& y: to_do){
                if (pair_crossings[x][y] == pair_crossings[y][x] && y > x){
                    next_excluded[y] = true;
                }
            }
            

            // Remove x and its twins and update all the structures
            

            vector<int> in_neighbors_reinsert;
            for (const int& j: to_do) {
                if (pair_crossings[x][j] - pair_crossings[j][x] > 0){
                    in_neighbors_reinsert.push_back(j);
                    for (const int& y: twins[x]){
                        auto it2 = find(in_neighbors[j].begin(), in_neighbors[j].end(), y);
                        in_neighbors[j].erase(it2);
                    }
                }
                
            }
            vector<int> out_neighbors_reinsert;
            for (const int& j: to_do) {
                if (pair_crossings[j][x] - pair_crossings[x][j] > 0){
                    out_neighbors_reinsert.push_back(j);
                    for (const int& y: twins[x]){
                        auto it2 = find(out_neighbors[j].begin(), out_neighbors[j].end(), y);
                        out_neighbors[j].erase(it2);
                    }
                }
            }

            


            vector<int> x_in_neighbors = in_neighbors[x];
            vector<int> x_out_neighbors = out_neighbors[x];

            vector<int> new_to_do;
            for (size_t j = 0; j < to_do.size(); ++j){
                if (are_equal( adj[x], adj[to_do[j]]) == false) 
                    new_to_do.push_back(to_do[j]);
            }

            for (const int& y: twins[x]){
                order.push_back(y);
                mask[y] = false;
                in_neighbors[y] = {};
                out_neighbors[y] = {};
            }

            aux3(adj, twins, new_to_do, order, best_order, best_bad_cr, new_current_bad_cr, pair_crossings, in_neighbors, out_neighbors, mask, depth+1, triangles_adj, 0, next_excluded, triangles);


            for (const int& y: twins[x]){
                order.pop_back();
                mask[y] = true;
                in_neighbors[y] = x_in_neighbors;
                out_neighbors[y] = x_out_neighbors;
            }


            for (const int& v: in_neighbors_reinsert){
                for (const int& y: twins[x]){
                    in_neighbors[v].push_back(y);
                }
            }
            for (const int& v: out_neighbors_reinsert){
                for (const int& y: twins[x]){
                    out_neighbors[v].push_back(y);
                }
            }


            // if (in_weights[x] == 0){
            //     // cout << "break on source" << endl;
            //     break;
            // }
        }
    }
}

int solver3(const vector<vector<int>>& adj, bool verbose) {
    if (verbose){
        cout << "############" << endl;
        cout << "solver3 d" << endl;
        cout << "nb vertices: " << adj.size() << endl;
    }
   
    vector<vector<int>> pair_crossings = compute_pair_crossings(adj);
    auto digraph = compute_directed_graph(adj);
    auto in_neighbors = digraph.first;
    auto out_neighbors = digraph.second;
    vector<vector<int>> twins = compute_twins(adj);
    for (int i = 0; i < twins.size(); ++i){
        print(twins[i]);
    }

    int bad_cr_total = 0;

    auto components = scc(digraph.second, digraph.first );

    for (int i = 0; i < components.size(); ++i){
        auto component = components[i];
        if (component.size() == 1) continue;
        cout << "compo " << i << "/" << components.size() << " size: " << component.size() << endl;
        vector<int> best_order = order_greedy_sequential_mask3(component, pair_crossings);
        vector<bool> mask(adj.size(), false);

        for (const int&x: component){
            mask[x] = true;
        }

        int lb = lower_bound_mask(pair_crossings, component);
        if (verbose){
            cout << "nb unavoidable crossings: " << lb << "\n";
        }

        int best_order_nc = nb_crossings_from_order(adj, best_order);
        int best_bad_cr = best_order_nc - lb;

        if (best_bad_cr == 0){
            if (verbose){
                cout << "min bad crossings: " << best_bad_cr << "\n";
                cout << "min crossings: " << lb + best_bad_cr << "\n";
            }   
            continue;
        }

        // Triangles
        int triangles_total = 0;
        vector<vector<int>> triangles_adj(adj.size());
        vector<vector<vector<int>>> triangles(adj.size());
        vector<bool> excluded(adj.size(), false);

        int current_bad_cr = 0;

        if (verbose){
            cout << "initial best bad cr " << best_bad_cr << endl;
        }

        vector<vector<int>> sub_in_neighbors(in_neighbors.size());
        for (const int& v: component){
            for (const int& w: in_neighbors[v]){
                if (mask[w]){
                    sub_in_neighbors[v].push_back(w);
                }
            }
        }
        vector<vector<int>> sub_out_neighbors(out_neighbors.size());
        for (const int& v: component){
            for (const int& w: out_neighbors[v]){
                if (mask[w]){
                    sub_out_neighbors[v].push_back(w);
                }
            }
        }


        vector<int> order;
        aux3(adj, twins, component, order, best_order, best_bad_cr, current_bad_cr, pair_crossings, sub_in_neighbors, sub_out_neighbors, mask, 0, triangles_adj, triangles_total, excluded, triangles);

        bad_cr_total += best_bad_cr;

        if (verbose){
            cout << "min bad crossings: " << best_bad_cr << "\n";
            cout << "min crossings: " << lb + best_bad_cr << "\n";
            print(best_order);
        }
    }

    cout << bad_cr_total << endl;

    // vector<int> to_do;
    // for (int i = 0; i < adj.size(); ++i){
    //     to_do.push_back(i);
    // }

    // vector<int> best_order = order_greedy_sequential_mask3(to_do, pair_crossings);


    // vector<bool> mask(adj.size(), true);

    // int lb = lower_bound(adj);
    // if (verbose){
    //     cout << "nb unavoidable crossings: " << lb << "\n";
    // }

    // int best_order_nc = nb_crossings_from_order(adj, best_order);
    // int best_bad_cr = best_order_nc - lb;

    // if (best_bad_cr == 0){
    //     if (verbose){
    //         cout << "min bad crossings: " << best_bad_cr << "\n";
    //         cout << "min crossings: " << lb + best_bad_cr << "\n";
    //     }   
    //     return best_bad_cr;
    // }

    // // Triangles
    // int triangles_total = 0;

    // vector<vector<int>> triangles_adj(adj.size());
    
    // vector<vector<vector<int>>> triangles(adj.size());
    // for (int x = 0; x < adj.size(); ++x){
    //     for (const int& y: digraph.first[x]){
    //         int yx = pair_crossings[y][x] - pair_crossings[x][y];
    //         for (const int& z: digraph.second[x]){
    //             if (pair_crossings[z][y] - pair_crossings[y][z] > 0){
    //                 int zy = pair_crossings[z][y] - pair_crossings[y][z];
    //                 int xz = pair_crossings[x][z] - pair_crossings[z][x];
    //                 int w = min(yx, min(zy, xz));
    //                 triangles[x].push_back({y,z,w});
    //             }
    //         }
    //     }
    //     sort( triangles[x].begin(), triangles[x].end(), [](vector<int> a, vector<int> b){
    //         return a[2] > b[2];
    //     });
    // }


    



    // vector<bool> excluded(adj.size(), false);

    // int current_bad_cr = 0;

    // if (verbose){
    //     cout << "initial best bad cr " << best_bad_cr << endl;
    // }


    // vector<int> order;
    // aux3(adj, twins, to_do, order, best_order, best_bad_cr, current_bad_cr, pair_crossings, digraph.first, digraph.second, mask, 0, triangles_adj, triangles_total, excluded, triangles);


    // if (verbose){
    //     cout << "min bad crossings: " << best_bad_cr << "\n";
    //     cout << "min crossings: " << lb + best_bad_cr << "\n";
    //     print(best_order);
    // }


    // Check solution
    // print(best_order);
    // cout << "best order size: " << best_order.size() << endl;
    // cout << "nb vertices: " << adj.size() << endl;
    // cout << "has duplicates? " << has_duplicates(best_order) << endl;
    // cout << nb_crossings_from_order(adj, best_order) << endl;
    // cout << nb_crossings_from_order(adj, best_order)-lb << endl;

    // return lb + best_bad_cr;
    return bad_cr_total;
}