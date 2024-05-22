
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
    int s,
    vector<int>& order, 
    vector<int>& best_order, 
    int& best_bad_cr, 
    int& current_bad_cr, 
    const vector<vector<int>>& pair_crossings, 
    const vector<vector<int>>& in_neighbors,
    const vector<vector<int>>& out_neighbors, 
    vector<bool>& mask, 
    int depth, 
    const vector<vector<int>>& triangles_adj, 
    const vector<int>& triangles_weight,
    vector<bool>& triangles_available,
    int triangles_total,
    vector<bool>& excluded,
    const vector<vector<vector<int>>>& triangles,
    vector<int>& in_weights,
    int last_source,
    vector<int>& sources,
    vector<int>& pos ){

        // cout << endl;
        // cout << "aux s=" << s << endl;
        // print(order);
        // print(to_do);
        // for (int i = 0; i < s; ++i){
        //     cout << i << " " << to_do[i] << " " << pos[to_do[i]] << endl;
        // }

    if (s == 0){
        if (current_bad_cr < best_bad_cr){
            // cout << "better " << current_bad_cr << "\n";
            best_order = order;
            best_bad_cr = current_bad_cr;
            // if ( has_duplicates(best_order)){
            //     cout << "order has duplicates" << endl;
            //     print(order);
            // }
        }
    } else {

        if (current_bad_cr + triangles_total >= best_bad_cr){
            return;
        }

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

        // sort(to_do.begin(), to_do.end(), [&in_weights](int a, int b) {
        //     return in_weights[a] < in_weights[b];
        // });

        // cout << string(depth, '-') << current_bad_cr << " todo " ;
        // print(to_do);

        // print(in_weights);
        // for (int j=0; j < adj.size(); ++j){
        //     cout << "mask " << j << " " << mask[j] << endl;
        // }

        // vector<int> sources;

        // for (int i = 0; i < s; ++i){
        //     int x = to_do[i];
        //     if (in_weights[x] == 0)
        //         sources.push_back(i);
        // }

        if (sources.size() > 0){

            sort(sources.begin(), sources.end());

            // cout << endl;
            // cout << "sources ";
            // print(sources);


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

            // vector<int> new_to_do;
            // for (int j = 0; j < s; ++j){
            //     int y = to_do[j];
            //     if (in_weights[y] > 0) 
            //         new_to_do.push_back(y);
            // }


            int sources_size = sources.size();
            vector<int> sources2;

            vector<int> sources_values;
            for (const int& id: sources){
                sources_values.push_back(to_do[id]);
            }



            // cout << "---" << endl;
            // print(order);
            // print(sources);


            for (int j = sources.size()-1; j >= 0; --j){
                int k = sources[j];
                int y = to_do[k];
                order.push_back(y);
                mask[y] = false;
                // out_neighbors[y] = {};
            }

            for (int j = sources.size()-1; j >= 0; --j){
                int k = sources[j];
                pos[to_do[s-1-(sources.size()-1-j)]] = k;
                swap( to_do[k], to_do[s-1-(sources.size()-1-j)]);

            }


          
            for (const int& y: sources_values){
                for (const int& z: out_neighbors[y]){
                    if (mask[z]){
                        in_weights[z] -= pair_crossings[y][z] - pair_crossings[z][y];
                        if (in_weights[z] == 0){
                            sources2.push_back(pos[z]);
                        }
                    }
                }
            }


            aux4(adj, twins, to_do, s- sources.size(), order, best_order, best_bad_cr, current_bad_cr, pair_crossings, in_neighbors, out_neighbors, mask, depth+1, triangles_adj, triangles_weight, triangles_available, triangles_total, excluded, triangles, in_weights, 0, sources2, pos);

            for (int j = 0; j < sources.size(); ++j){
                int k = sources[j];
                swap( to_do[k], to_do[s-1-(sources.size()-1-j)]);
            }

            for (int j = 0; j < sources.size(); ++j){
                pos[to_do[s-1-(sources.size()-1-j)]] = s-1-(sources.size()-1-j);
            }

            for (int j = sources.size()-1; j >= 0; --j){
                int k = sources[j];
                int y = to_do[k];
                order.pop_back();
                mask[y] = true;
                // out_neighbors[y] = old_out_neighbors[y];
            }


            for (const int& j: sources){
                int y = to_do[j];
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
       

        vector<int> triangles_to_reinsert;

        // Branch on to_do
        // cout << "branch " << to_do.size() << "\n";
        for (int i = 0; i < s; ++i){
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
                // cout << string(depth, '-') << "branch continue " << x << endl;
                // break; 
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

                
                if (last_source >= 2){
                    if (xz == w1 && xz == w2){
                        if (!(x > y && x > z)){
                            // cout << string(depth, '-') << "brafdfdnch cut non optimal last third " << x << endl;
                            continue;
                        }
                    }
                    if (xz == mw && xz == w1){
                        if (!(x > z)){
                            // cout << string(depth, '-') << "braaaaaaanch cut non optimal last third " << x << endl;
                            continue;
                        }
                    } else if ( xz == mw && xz == w2 && !(x > y)){
                        // cout << string(depth, '-') << "branch cut non opfdfdfdtimal last third " << x << endl;
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

                int ab = pair_crossings[a][b] - pair_crossings[b][a]; 
                int cx = pair_crossings[c][x] - pair_crossings[x][c];
                

                if ( last_source >= 3 && ab > 0 && bc > 0 && cx > 0 &&
                ca > 0 && xb > 0 && xa > 0 &&
                  bc == ca + (xb + xa)){
                     continue;
                }
            }

            if (order.size() >= 4){
                int a = order[order.size()-4];
                int b = order[order.size()-3];
                int c = order[order.size()-2];
                int d = order[order.size()-1];

                int db = pair_crossings[d][b] - pair_crossings[b][d];
                int ca = pair_crossings[c][a] - pair_crossings[a][c];
                // int xa = pair_crossings[x][a] - pair_crossings[a][x];
                int xc = pair_crossings[x][c] - pair_crossings[c][x];
                int ax = pair_crossings[a][x] - pair_crossings[x][a];

                int cd = pair_crossings[c][d] - pair_crossings[d][c];

                int xb = pair_crossings[x][b] - pair_crossings[b][x];
                if (xb < 0) xb = 0;
                int da = pair_crossings[d][a] - pair_crossings[a][d];
                if (da < 0) da = 0;

                if (ax >= 0 && db >= 0 && ca >= 0 && xc >= 0 && ca + cd < db + ca + xc + (xb + da)){
                    // a b c d x is not optimal
                    // (a d) (b x) c is better (paranthesis indicates that we can swap a and d depending on the arc ad)
                    // same with bx
                    continue;
                }

            }


            int sum = 0;
            bool optimal = true;
            int jmin = ((int)(order.size()) - 50 >= 0) ? order.size() - 50 : 0;
            for (int j = order.size()-1;  j >= jmin; --j){
                int y = order[j];
                int xy = pair_crossings[x][y] - pair_crossings[y][x];
                if (xy > 0){
                    sum -= xy;
                } else if (xy < 0){
                    sum += -xy;
                }
                if (sum < 0) {
                    optimal = false;
                    break;
                }
            }
            if (optimal == false){
                // cout << string(depth, '-') << "continue non optimal " << x << endl;
                continue;
            } 

            

            
            

            int new_current_bad_cr = current_bad_cr + x_bad_cr;
            

            
            // cout << string(depth, '-') << "branch " << x << " " << x_bad_cr << " "<< new_current_bad_cr <<  " in neighbors: ";
            // print(in_neighbors[x]);
            // print(in_weights);

            

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

            // vector<int> new_to_do;
            // for (const int& y: to_do){
            //     if (y != x)
            //         new_to_do.push_back(y);
            // }

            
            sources.clear();
            // vector<int> sources2;

            for (const int& v: out_neighbors[x]){
                if (mask[v]){
                    in_weights[v] -= pair_crossings[x][v] - pair_crossings[v][x];
                    if (in_weights[v] == 0){
                        if (pos[v] == s-1){
                            sources.push_back(i); // because of the swap
                        } else {
                            sources.push_back(pos[v]);
                        }
                    }
                }
            }

            // cout << "--- branch " << "x="<<  x << " i=" << i <<  endl;
            // print(sources);

            triangles_to_reinsert.clear();
            int new_triangles_total = triangles_total;
            for (const int& triangle_id: triangles_adj[x]){
                if (triangles_available[triangle_id]){
                    new_triangles_total -= triangles_weight[triangle_id];
                    triangles_available[triangle_id] = false;
                    triangles_to_reinsert.push_back(triangle_id);
                }
            }

            order.push_back(x);
            mask[x] = false;
            // in_neighbors[x] = {};
            // out_neighbors[x] = {};
            
            pos[to_do[s-1]] = i;
            pos[to_do[i]] = s-1;

            swap(to_do[i], to_do[s-1]);

        


            aux4(adj, twins, to_do, s-1, order, best_order, best_bad_cr, new_current_bad_cr, pair_crossings, in_neighbors, out_neighbors, mask, depth+1, triangles_adj, triangles_weight, triangles_available, new_triangles_total, excluded, triangles, in_weights, last_source+1, sources, pos);

            swap(to_do[i], to_do[s-1]);
            pos[to_do[s-1]] = s-1;
            pos[to_do[i]] = i;

            order.pop_back();
            mask[x] = true;
            // in_neighbors[x] = x_in_neighbors;
            // out_neighbors[x] = x_out_neighbors;

            for (const int& triangle_id: triangles_to_reinsert){
                triangles_available[triangle_id] = true;
            }

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

vector<int> solver4(const vector<vector<int>>& adj, bool verbose) {
    if (verbose){
        cout << "############" << endl;
        cout << "solver4 E" << endl;
        cout << "nb vertices: " << adj.size() << endl;
    }

    vector<int> final_order;

   
    vector<vector<int>> pair_crossings = compute_pair_crossings(adj);
    auto digraph = compute_directed_graph(adj);
    auto in_neighbors = digraph.first;
    auto out_neighbors = digraph.second;
    vector<vector<int>> twins = compute_twins(adj);
    // for (int i = 0; i < twins.size(); ++i){
    //     print(twins[i]);
    // }

    // print_adj(adj);

    // to_dot(out_neighbors, pair_crossings);


    int bad_cr_total = 0;

    auto components = scc(digraph.second, digraph.first );

    for (int i = 0; i < components.size(); ++i){


        auto raw_component = components[i];
        if (raw_component.size() == 1){
            final_order.push_back(raw_component[0]);
            continue;
        }
        if (verbose){
            cout << endl;
            cout << "compo " << i << "/" << components.size() << " size: " << raw_component.size() << endl;
        }
        
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

        if (verbose){
            cout << "reduced compo size: " << component.size() << endl;
        }
        

        // Structure
        vector<bool> mask(adj.size(), false);
        for (const int&x: component){
            mask[x] = true;
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

        // sort(component.begin(), component.end(), [&in_weights](auto a, auto b){
        //     return in_weights[a] < in_weights[b];
        // });


        // Upper bound
        vector<int> best_order = order_greedy_sequential_mask3(component, pair_crossings);

        

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

            for (int j = 0; j < best_order.size(); ++j){
                for (const int& twin: twins[best_order[j]]){
                    final_order.push_back(twin);
                }
            }
            continue;
        }

        // Obsolete and useless
        vector<bool> excluded(adj.size(), false);

        // Triangles
        vector<vector<vector<int>>> triangles(adj.size());
        auto r = find_edge_disjoint_triangles_greedy3(pair_crossings, sub_in_neighbors, sub_out_neighbors, component);
        vector<vector<int>> triangles_adj = r.first;
        vector<int> triangles_weight = r.second;
        vector<bool> triangles_available(triangles_weight.size(), true);

        int triangles_total = 0;
        for (const int& w: triangles_weight){
            triangles_total += w;
        }

        if (verbose){
            cout << "triangles: " << triangles_weight.size() << " total_weight: " << triangles_total << endl;
        }
        // for (const int& x: component){
        //     cout << "triangles_adj[" << x << "] ";
        //     print(triangles_adj[x]);
        // }
        // for (int j = 0; j < triangles_weight.size(); ++j){
        //     printf("triangles_weight[%d] %d\n", j, triangles_weight[j]);
        // }



        int current_bad_cr = 0;

        if (verbose){
            cout << "initial best bad cr " << best_bad_cr << endl;
        }

        

        int s = component.size();

        vector<int> pos(adj.size());
        for(int j=0; j < component.size(); ++j){
            pos[component[j]] = j;
        }
        
        vector<int> sources;
        vector<int> order;

        aux4(adj, twins, component, s, order, best_order, best_bad_cr, current_bad_cr, pair_crossings, sub_in_neighbors, sub_out_neighbors, mask, 0, triangles_adj, triangles_weight, triangles_available, triangles_total, excluded, triangles, in_weights, adj.size(), sources, pos);

        bad_cr_total += best_bad_cr;

        if (verbose){
            cout << "min bad crossings: " << best_bad_cr << "\n";
            cout << "min crossings: " << lb + best_bad_cr << "\n";
            print(best_order);
        }

        for (int j = 0; j < best_order.size(); ++j){
            for (const int& twin: twins[best_order[j]]){
                final_order.push_back(twin);
            }
        }
    }

    if (verbose){
        cout << "min bad crossings: " << bad_cr_total << endl;
    }

    // return lb + best_bad_cr;
    return final_order;
}