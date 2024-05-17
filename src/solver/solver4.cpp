
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
void aux4(const vector<vector<int>>& adj, 
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
    const vector<vector<vector<int>>>& triangles,
    vector<int>& in_weights,
    int last_source ){


    if (to_do.size() == 0){
        if (current_bad_cr < best_bad_cr){
            cout << "better " << current_bad_cr << "\n";
            best_order = order;
            best_bad_cr = current_bad_cr;
        }
    } else {

        // if (best_bad_cr - current_bad_cr <= 30){
        //     vector<int> t(to_do.size(), numeric_limits<int>::max());
        //     for (const int& v: to_do){
        //         vector<int> v_in_weights;
        //         for (const int& w: in_neighbors[v]){
        //             v_in_weights.push_back(pair_crossings[w][v] - pair_crossings[v][w]);
        //         }
        //         sort(v_in_weights.begin(), v_in_weights.end(), [](int a, int b) {
        //             return a > b;
        //         });
        //         int accu = 0;
        //         for (int i = v_in_weights.size()-1; i >= 0; --i){
        //             accu += v_in_weights[i];
        //             t[i] = min(t[i], accu);
        //         }
        //         for (int i = v_in_weights.size(); i < to_do.size(); ++i){
        //             t[i] = 0;
        //         }
        //     }

        //     int lower_bound0 = 0;
        //     for (int i = 0; i < t.size() && t[i] > 0; ++i){
        //         lower_bound0 += t[i];
        //     }
        //     if (current_bad_cr + lower_bound0 >= best_bad_cr) return;
        // }

        // if (best_bad_cr - current_bad_cr <= 20){
        //     vector<int> indegrees;
        //     for (const int& v: to_do){
        //         int vind = 0;
        //         for (const int&w: in_neighbors[v]){
        //             if (mask[w]) ++vind;
        //         }
        //         // indegrees.push_back(in_neighbors[v].size());
        //         indegrees.push_back(vind);
        //     }
        //     sort(indegrees.begin(), indegrees.end());
        //     int rind = 0;
        //     for (int i = 0; i < to_do.size()/2; ++i){
        //         if (indegrees[i] >= i){
        //             rind += indegrees[i] - i;
        //         }
        //     }
        //     if (current_bad_cr + rind >= best_bad_cr) return;
        // }
        

        // vector<int> in_weights(adj.size());
        // for( const int& x: to_do){
        //     for (const int& v: in_neighbors[x]){
        //         if (mask[v])
        //         in_weights[x] += pair_crossings[v][x] - pair_crossings[x][v];
        //     }
        // }

        sort(to_do.begin(), to_do.end(), [&in_weights](int a, int b) {
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

            // cout << string(depth, '-') << "sources ";
            // print(sources);
            

            // vector<vector<int>> old_in_neighbors(adj.size());
            // for (const int& y: to_do){
            //     if (in_weights[y] > 0)
            //         old_in_neighbors[y] = in_neighbors[y];
            // }
            
            // vector<vector<int>> old_out_neighbors(adj.size());
            // for (const int& y: sources){
            //     old_out_neighbors[y] = out_neighbors[y];
            // }

            // for (const int& v: to_do){
            //     if (in_weights[v] > 0){
            //         in_neighbors[v] = vector<int>();
            //         for (const int& w: old_in_neighbors[v]){
            //             if (in_weights[w] > 0){
            //                 in_neighbors[v].push_back(w);
            //             }
            //         }
            //     }
            // }

            vector<int> new_to_do;
            for (const int& y: to_do){
                if (in_weights[y] > 0) 
                    new_to_do.push_back(y);
            }

            for (const int& y: sources){
                for (const int& z: out_neighbors[y]){
                    if (mask[z])
                    in_weights[z] -= pair_crossings[y][z] - pair_crossings[z][y];
                }
            }

            for (const int& y: sources){
                order.push_back(y);
                mask[y] = false;
                // out_neighbors[y] = {};
            }

            // vector<bool> next_excluded(adj.size(), false);

            aux4(adj, twins, new_to_do, order, best_order, best_bad_cr, current_bad_cr, pair_crossings, in_neighbors, out_neighbors, mask, depth+1, triangles_adj, 0, excluded, triangles, in_weights, 0);


            for (const int& y: sources){
                order.pop_back();
                mask[y] = true;
                // out_neighbors[y] = old_out_neighbors[y];
            }

            for (const int& y: sources){
                for (const int& z: out_neighbors[y]){
                    if (mask[z])
                    in_weights[z] += pair_crossings[y][z] - pair_crossings[z][y];
                }
            }


            // for (const int& v: to_do){
            //     if (in_weights[v] > 0){
            //         in_neighbors[v] = old_in_neighbors[v];
            //     }
            // }
            


            return;
        }
       


        // Branch on to_do
        // cout << "branch " << to_do.size() << "\n";
        for (int i = 0; i < to_do.size(); ++i){
            int x = to_do[i];
            // if (excluded[x]){
            //     // cout << string(depth, '-') << "branch exclude " << x << endl;
            //     continue;
            // } 

            if (order.size() >= 1){
                int y = order[order.size()-1];
                if (pair_crossings[x][y] - pair_crossings[y][x] > 0){
                    // cout << string(depth, '-') << "branch exclude out-neighbor " << x << endl;
                    continue;
                } 
                else if ( last_source >= 1 && pair_crossings[x][y] == pair_crossings[y][x] && y > x) {
                    // cout << string(depth, '-') << "branch exclude == " << x << " after " << y << endl;
                    continue;
                }
            }

            int x_bad_cr = in_weights[x];

            // Only if to_do is sorted by increasing in_weights
            if (current_bad_cr + x_bad_cr >= best_bad_cr){
                // cout << string(depth, '-') << "branch break " << x << endl;
                break; 
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

                if (last_source >= 2){
                    if (xz == w1 && xz == w2){
                        if (!(x > y && x > z)){
                            continue;
                        }
                    }
                    if (xz == mw && xz == w1){
                        if (!(x > z)){
                            continue;
                        }
                    } else if ( xz == mw && xz == w2 && !(x > y)){
                        continue;
                    }
                }
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
                
                if ( ca >= 0 && xb >= 0  && bc < ca + (xb + xa) ){
                    // cout << string(depth, '-') << "branch cut non optimal last fourth " << x << endl;
                    // c x a b or c a x b is better than a b c x and has only one BA: bc
                    // we are sure that there are arcs ab bc and cd
                    continue;
                }

                if ( last_source >= 3 && ca >= 0 && xb >= 0 && bc == ca + (xb + xa)){
                    continue;
                }

                int cx = pair_crossings[c][x] - pair_crossings[x][c];
                if ( ca <= 0 && xb >= 0  && cx <  xb + xa ){
                    // cout << string(depth, '-') << "branch cut non optimal last fourth " << x << endl;

                    // cout << string(depth, '-') << "branch cut non optimal last third " << x << endl;
                    // c x a b or c a x b is better than a b c x and has only one BA: bc
                    // we are sure that there are arcs ab bc and cd
                    continue;
                }

                // if (last_source >= 3 && ca <= 0 && xb > 0 && xa > 0 && cx == xb + xa){
                //     continue;
                // }
            }

            // if (order.size() >= 4){
            //     int a = order[order.size()-4];
            //     int b = order[order.size()-3];
            //     int c = order[order.size()-2];
            //     int d = order[order.size()-1];

            //     int db = pair_crossings[d][b] - pair_crossings[b][d];
            //     int ca = pair_crossings[c][a] - pair_crossings[a][c];
            //     // int xa = pair_crossings[x][a] - pair_crossings[a][x];
            //     int xc = pair_crossings[x][c] - pair_crossings[c][x];
            //     int ax = pair_crossings[a][x] - pair_crossings[x][a];

            //     int cd = pair_crossings[c][d] - pair_crossings[d][c];

            //     int xb = pair_crossings[x][b] - pair_crossings[b][x];
            //     if (xb < 0) xb = 0;
            //     int da = pair_crossings[d][a] - pair_crossings[a][d];
            //     if (da < 0) da = 0;

            //     if (ax >= 0 && db >= 0 && ca >= 0 && xc >= 0 && ca + cd < db + ca + xc + (xb + da)){
            //         // a b c d x is not optimal
            //         // (a d) (b x) c is better (paranthesis indicates that we can swap a and d depending on the arc ad)
            //         // same with bx
            //         continue;
            //     }

            // }


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
            if (optimal == false){
                // cout << string(depth, '-') << "continue non optimal " << x << endl;
                continue;
            } 

            

            
            

            int new_current_bad_cr = current_bad_cr + x_bad_cr;


            
            // cout << string(depth, '-') << "branch " << x << " in neighbors: ";
            // print(in_neighbors[x]);

            

            // Remove x  and update all the structures
            

            // for (const int& j: out_neighbors[x]){
            //     // in_weights[j] -= pair_crossings[x][j] - pair_crossings[j][x];
            //     auto it2 = find(in_neighbors[j].begin(), in_neighbors[j].end(), x);
            //     in_neighbors[j].erase(it2);
            // }
            // for (const int& j: in_neighbors[x]){
            //     auto it2 = find(out_neighbors[j].begin(), out_neighbors[j].end(), x);
            //     out_neighbors[j].erase(it2);
            // }
            

            


            // vector<int> x_in_neighbors = in_neighbors[x];
            // vector<int> x_out_neighbors = out_neighbors[x];

            vector<int> new_to_do;
            for (const int& y: to_do){
                if (y != x)
                    new_to_do.push_back(y);
            }

            for (const int& v: out_neighbors[x]){
                if (mask[v])
                    in_weights[v] -= pair_crossings[x][v] - pair_crossings[v][x];
            }

            order.push_back(x);
            mask[x] = false;
            // in_neighbors[x] = {};
            // out_neighbors[x] = {};

            aux4(adj, twins, new_to_do, order, best_order, best_bad_cr, new_current_bad_cr, pair_crossings, in_neighbors, out_neighbors, mask, depth+1, triangles_adj, 0, excluded, triangles, in_weights, last_source+1);

            order.pop_back();
            mask[x] = true;
            // in_neighbors[x] = x_in_neighbors;
            // out_neighbors[x] = x_out_neighbors;

            for (const int& v: out_neighbors[x]){
                if (mask[v])
                    in_weights[v] += pair_crossings[x][v] - pair_crossings[v][x];
            }

            // for (const int& v: out_neighbors[x]){
            //     // in_weights[v] += pair_crossings[x][v] - pair_crossings[v][x];
            //     in_neighbors[v].push_back(x);
            // }
            // for (const int& v: in_neighbors[x]){
            //     out_neighbors[v].push_back(x);
            // }


        }
    }
}

int solver4(const vector<vector<int>>& adj, bool verbose) {
    if (verbose){
        cout << "############" << endl;
        cout << "solver4 A" << endl;
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

    print_adj(adj);


    int bad_cr_total = 0;

    auto components = scc(digraph.second, digraph.first );

    for (int i = 0; i < components.size(); ++i){
        auto raw_component = components[i];
        if (raw_component.size() == 1) continue;
        cout << endl;
        cout << "compo " << i << "/" << components.size() << " size: " << raw_component.size() << endl;
        
        // Twin reduction
        vector<int> component;
        for (const int& x: raw_component){
            if (twins[x][0] == x){ // twins[x] is sorted increasing
                component.push_back(x);
            }
        }

        for (int j = 0; j < component.size(); ++j){
            int x = component[j];
            for (int k = 0; k < component.size(); ++k){
                int y = component[k];
                pair_crossings[x][y] *= twins[x].size() * twins[y].size();
            }
        }

        cout << "reduced compo size: " << component.size() << endl;
        

        vector<int> best_order = order_greedy_sequential_mask3(component, pair_crossings);
        vector<bool> mask(adj.size(), false);

        for (const int&x: component){
            mask[x] = true;
        }

        int lb = lower_bound_mask(pair_crossings, component);
        if (verbose){
            cout << "nb unavoidable crossings: " << lb << "\n";
        }

        int best_order_nc = nb_crossings_from_order2(best_order, pair_crossings);
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

        vector<int> in_weights(adj.size());
        for( const int& x: component){
            for (const int& v: sub_in_neighbors[x]){
                in_weights[x] += pair_crossings[v][x] - pair_crossings[x][v];
            }
        }


        vector<int> order;
        aux4(adj, twins, component, order, best_order, best_bad_cr, current_bad_cr, pair_crossings, sub_in_neighbors, sub_out_neighbors, mask, 0, triangles_adj, triangles_total, excluded, triangles, in_weights, adj.size());

        bad_cr_total += best_bad_cr;

        if (verbose){
            cout << "min bad crossings: " << best_bad_cr << "\n";
            cout << "min crossings: " << lb + best_bad_cr << "\n";
            print(best_order);
        }
    }

    cout << bad_cr_total << endl;

    // return lb + best_bad_cr;
    return bad_cr_total;
}