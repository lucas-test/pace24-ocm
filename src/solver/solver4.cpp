
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


pair<vector<int>,int> iterative_greedy_order(const vector<int>& component, const vector<vector<int>>& pair_crossings, int lb){
    vector<int> order = order_greedy_sequential_mask3(component, pair_crossings);
    int bad_cr = nb_crossings_from_order2(order, pair_crossings) - lb;
    // cout << "init bad cr: " << bad_cr << "\n";
    while (true){
        vector<int> new_order = order_greedy_sequential_mask3(order, pair_crossings);
        int new_bad_cr = nb_crossings_from_order2(new_order, pair_crossings) - lb;
        if (new_bad_cr < bad_cr){
            order = new_order;
            bad_cr = new_bad_cr;
            // cout << "better bad cr: " << bad_cr << "\n";
        } else {
            return make_pair(order, bad_cr);
        }
    }
}



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
    const vector<int>& branch_order,
    vector<int>& to_do, 
    int s,
    vector<int>& order, // left order
    vector<int>& right_order,
    vector<int>& best_order, 
    int& best_bad_cr, 
    int& current_bad_cr, 
    const vector<vector<int>>& pair_diff, 
    const vector<vector<int>>& pair_crossings,
    const vector<vector<int>>& in_neighbors,
    const vector<vector<int>>& out_neighbors, 
    vector<bool>& mask, 
    int depth, 
    const vector<vector<int>>& triangles_adj, 
    const vector<int>& triangles_weight,
    vector<bool>& triangles_available,
    int triangles_total,
    vector<int>& in_weights,
    vector<int>& out_weights,
    int last_source,
    vector<int>& sources,
    vector<int>& pos,
    vector<int>& min_left_sum, // for every vertex to insert
    vector<int>& right_sum, // for every vertex already inserted (in order)
    vector<int>& indegrees,
    vector<int>& indegrees1_indices ){

        // cout << endl;
        // cout << "aux s=" << s << endl;
        // print(order);
        // print(to_do);
        // for (int i = 0; i < s; ++i){
        //     cout << i << " " << to_do[i] << " " << pos[to_do[i]] << endl;
        // }

    if (s == 0){
        if (current_bad_cr < best_bad_cr){
            best_order = order;
            best_bad_cr = current_bad_cr;

            int lb = lower_bound_mask(pair_crossings, best_order);

            auto r = iterative_greedy_order(best_order, pair_crossings, lb);
            best_order = r.first;
            best_bad_cr = r.second;
            // cout << "better " << best_bad_cr << "\n";


            // if ( has_duplicates(best_order)){
            //     cout << "order has duplicates" << endl;
            //     print(order);
            // }
        }
    } else {

        if (current_bad_cr + triangles_total >= best_bad_cr){
            return;
        }


        


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
                if (min_left_sum[to_do[id]] < 0) return;
            }

            // Optimize this ? only 1 iteration over order ?
            for (const int& v: order){
                int sum = 0;
                for (const int& source: sources_values){
                    sum += pair_diff[v][source];
                }
                if (right_sum[v] + sum < 0){
                    return;
                }
            }

            for (const int& v: order){
                for (const int& source: sources_values){
                    right_sum[v] += pair_diff[v][source];
                }
            }


            vector<int> new_min_left_sum(adj.size());
            for (int j = 0; j < s; ++j){
                int x = to_do[j];
                if (in_weights[x] == 0) continue;
                int sum = 0;
                for (const int& source: sources_values){
                    sum -= pair_diff[x][source];
                }
                if (sum < min_left_sum[x] + sum){
                    new_min_left_sum[x] = sum;
                } else {
                    new_min_left_sum[x] = min_left_sum[x] + sum;
                }
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
                        indegrees[z] --;
                        in_weights[z] -= pair_diff[y][z] ;
                        if (in_weights[z] == 0){
                            sources2.push_back(pos[z]);
                        }
                    }
                }
            }

            vector<int> new_indegrees1_indices;
            for (int k = 0; k < s-sources.size(); ++k){
                int v = to_do[k];
                if (indegrees[v] == 1){
                    new_indegrees1_indices.push_back(k);
                }
            }


            aux4(adj, branch_order, to_do, s- sources.size(), order, right_order, best_order, best_bad_cr, current_bad_cr, pair_diff, pair_crossings, in_neighbors, out_neighbors, mask, depth+1, triangles_adj, triangles_weight, triangles_available, triangles_total, in_weights, out_weights, 0, sources2, pos, new_min_left_sum, right_sum, indegrees, new_indegrees1_indices);

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

             for (const int& v: order){
                for (const int& source: sources_values){
                    right_sum[v] -= pair_diff[v][source];
                }
            }

            for (const int& j: sources){
                int y = to_do[j];
                for (const int& z: out_neighbors[y]){
                    if (mask[z]){
                        indegrees[z] ++;
                        in_weights[z] += pair_diff[y][z] ;
                    }
                }
            }


            // for (const int& v: to_do){
            //     if (in_weights[v] > 0){
            //         in_neighbors[v] = old_in_neighbors[v];
            //     }
            // }
            


            return;
        }

        // if (indegrees1_indices.size() == 3){
        //     int x = to_do[indegrees1_indices[0]];
        //     int y = to_do[indegrees1_indices[1]];
        //     int z = to_do[indegrees1_indices[2]];


        //     if (pair_diff[x][y] > 0){
        //         if (pair_diff[y][z] > 0 && pair_diff[z][x] > 0){
        //     cout << x << " " << y << " " << z << endl;

        //             cout << "triangle dominating" << endl;
        //         }
        //     } else if (pair_diff[x][y] < 0){
        //         if (pair_diff[z][y] > 0 && pair_diff[x][z] > 0){
        //     cout << x << " " << y << " " << z << endl;

        //             cout << "triangle dominating" << endl;
        //         }
        //     }
        // }




        // Search for sinks
        // for (int i = 0; i < s; ++i){
        //     int x = to_do[i];

        //     if (out_weights[x] == 0){
        //         cout << "sink " << x << " depth: " << depth<< endl;
        //     }
        // }









        

        // vector<int> vertices;
        // for (int i = 0; i < s ; ++i){
        //     vertices.push_back(to_do[i]);
        // }

        // vector<vector<int>> current_in_neighbors(in_neighbors.size());
        // for (const int& v: vertices){
        //     for (const int& w: in_neighbors[v]){
        //         if (mask[w]){
        //             current_in_neighbors[v].push_back(w);
        //         }
        //     }
        // }
        // vector<vector<int>> current_out_neighbors(out_neighbors.size());
        // for (const int& v: vertices){
        //     for (const int& w: out_neighbors[v]){
        //         if (mask[w]){
        //             current_out_neighbors[v].push_back(w);
        //         }
        //     }
        // }

        // vector<vector<int>> components = scc_sub_digraph(current_out_neighbors, current_in_neighbors, vertices);

        // if (components.size() >= 2){
        //     // cout << "s=" << s << endl;
        //     // print(to_do);
        //     // print(vertices);

        //     if (components[0].size() == 3){
        //         cout << "scc: "<< components.size() << " depth: " << depth << endl;
        //         vector<int> components_sizes;
        //         for (const vector<int>& component: components){
        //             components_sizes.push_back(component.size());
        //         }
        //         // sort(components_sizes.begin(), components_sizes.end());
        //         print(components_sizes);
        //     }
        // }



       

        vector<int> triangles_to_reinsert;
        vector<int> new_min_left_sum(adj.size());
        vector<int> new_indegrees1_indices;

        // Branch on to_do
        // cout << "branch " << to_do.size() << "\n";
        // for (int ii = 0; ii < s; ++ii){
            // int ix = ii;

        for (int i = 0; i < s; ++i){
            int x = to_do[i];
            if (-min_left_sum[x] > in_weights[x]){
                return;
            }
        }
        
        for(int ii = 0; ii < branch_order.size(); ++ii){
            if (mask[branch_order[ii]] == false) continue;
            int ix = pos[branch_order[ii]];

            int x = to_do[ix];

            if (min_left_sum[x] < 0){
                continue;
            }

            // inserting x is not optimal, because putting it at the end is better
            if (in_weights[x] > out_weights[x]){
                continue;
            }

            // bool optimal = true;
            // for (const int& v: order){
            //     if (right_sum[v] + pair_diff[v][x] < 0){ // v can be moved after x
            //         optimal = false;
            //         break;
            //     }
            // }
            // if (optimal == false){
            //     continue;
            // }

            if (order.size() >= 1){
                int y = order[order.size()-1];
                if (pair_diff[x][y] > 0){
                    // cout << string(depth, '-') << "branch exclude out-neighbor " << x << endl;
                    continue;
                } 
                else if ( last_source >= 1 && pair_diff[x][y] == 0 && y > x) {
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
                int xz = pair_diff[x][z];
                int w1 = pair_diff[z][y];
                int w2 = pair_diff[y][x];
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
                int bc = pair_diff[b][c];
                int ca = pair_diff[c][a] ;
                // int xa = pair_diff[x][a];
                int xb = pair_diff[x][b] ;

                int xa = pair_diff[x][a];
                if (xa < 0) xa = 0;
                
                if ( ca >= 0 && xb >= 0  && bc < ca + (xb + xa) ){
                    // cout << string(depth, '-') << "branch cut non optimal last fourth " << x << endl;
                    // c x a b or c a x b is better than a b c x and has only one BA: bc
                    // we are sure that there are arcs ab bc and cd
                    continue;
                }


                int ab = pair_diff[a][b];
                xa = pair_diff[x][a];
                if ( xa >= 0 && ca >= 0 && ab < xa + ca){
                    // a b c x not optimal
                    // b c x a is better
                    continue;
                }

                int cx = pair_diff[c][x] ;
                
                if (last_source >= 3 && bc == ca + (xb + xa) && (xb >0 || (xb == 0 && x < b)  ) ){
                    continue;
                }

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

                int db = pair_diff[d][b] ;
                int ca = pair_diff[c][a] ;
                int xc = pair_diff[x][c] ;
                int ax = pair_diff[a][x] ;
                int cd = pair_diff[c][d];

                int xb = pair_diff[x][b];
                if (xb < 0) xb = 0;
                int da = pair_diff[d][a];
                if (da < 0) da = 0;

                if (ax >= 0 && db >= 0 && ca >= 0 && xc >= 0 && ca + cd < db + ca + xc + (xb + da)){
                    // a b c d x is not optimal
                    // (a d) (b x) c is better (paranthesis indicates that we can swap a and d depending on the arc ad)
                    // same with bx
                    continue;
                }

                xb = pair_diff[x][b];
                int bc = pair_diff[b][c];
                int cx = pair_diff[c][x];
                if (xb >= 0 && db >= 0 && ca >= 0 && ax >= 0 && cx >= 0 && 
                 bc < ca + db + xb){
                    // a b c d x is not optimal
                    // c (a d) x b is better
                    // paranthesis indicates that we can swap a and d
                    continue;
                }

                da = pair_diff[d][a];
                int ab = pair_diff[a][b];
                int ac = pair_diff[a][c];
                int xa = pair_diff[x][a];
                if ( ac >= 0 && xa >= 0 && da >= 0 && ab + ac < xa + da ){
                    // a b c d x not optimal
                    // b c d x a is strictly better
                    continue;
                }


            }

            if (order.size() >= 2){
                int d = order[order.size()-1];
                int c = order[order.size()-2];
                int sumd = 0;
                int sumc = 0;
                bool optimal = true;
                int jmin = ((int)(order.size()) - 50 >= 0) ? order.size() - 50 : 0;
                for (int j = order.size()-2;  j >= jmin; --j){
                    int y = order[j];
                    sumd -= pair_diff[x][y];
                    sumd += -pair_diff[d][y];
                    if (sumd < 0){
                        optimal = false;
                        break;
                    }
                    if (j <= order.size()-3){
                        sumc -= pair_diff[x][y] + pair_diff[d][y] + pair_diff[c][y];
                        if (sumc < 0){
                            optimal = false;
                            break;
                        }
                    }
                }
                if (optimal == false){
                    // cout << string(depth, '-') << "continue non optimal " << x << endl;
                    continue;
                } 
            }

            

            

            
            

            int new_current_bad_cr = current_bad_cr + x_bad_cr;
            

            
            // cout << string(depth, '-') << "branch " << x << " " << x_bad_cr << " "<< new_current_bad_cr <<  " in neighbors: ";
            // print(in_neighbors[x]);
            // print(in_weights);

            

            // Remove x  and update all the structures

            bool optimal = true;
            for (const int& v: order){
                right_sum[v] += pair_diff[v][x];
                if (right_sum[v] < 0){
                    optimal = false;
                }
            }
            if (optimal == false){
                // Reverse changes on right_sum and continue
                for (const int& v: order){
                    right_sum[v] -= pair_diff[v][x];
                }
                continue;
            }

            for (int j = 0; j < s; ++j){
                int y = to_do[j];
                if (x != y){
                    if (-pair_diff[y][x] < min_left_sum[y] - pair_diff[y][x]){
                        new_min_left_sum[y] = -pair_diff[y][x];
                    } else {
                        new_min_left_sum[y] = min_left_sum[y] - pair_diff[y][x];
                    }
                }
            }

            
            

            // for (const int& j: out_neighbors[x]){
            //     // in_weights[j] -= pair_diff[x][j] ;
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

            

            // Recompute sources
            sources.clear();
            for (const int& v: out_neighbors[x]){
                if (mask[v]){
                    in_weights[v] -= pair_diff[x][v];
                    if (in_weights[v] == 0){
                        if (pos[v] == s-1){
                            sources.push_back(ix); // because of the swap
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
            
            pos[to_do[s-1]] = ix;
            pos[to_do[ix]] = s-1;

            swap(to_do[ix], to_do[s-1]);

             // Recompute indegrees1 vertices
            new_indegrees1_indices.clear();
            for (const int& v: out_neighbors[x]){
                if (mask[v]){
                    -- indegrees[v];
                    // if (indegrees[v] == 1){
                    //     new_indegrees1_indices.push_back(pos[v]);
                    // }
                }
            }
            for (int k = 0; k < s-1; ++k){
                int v = to_do[k];
                if (indegrees[v] == 1){
                    new_indegrees1_indices.push_back(k);
                }
            }

            for (const int& v: in_neighbors[x]){
                if (mask[v]){
                    out_weights[v] -= pair_diff[v][x];
                }
            }


            aux4(adj, branch_order, to_do, s-1, order, right_order, best_order, best_bad_cr, new_current_bad_cr, pair_diff, pair_crossings, in_neighbors, out_neighbors, mask, depth+1, triangles_adj, triangles_weight, triangles_available, new_triangles_total,  in_weights, out_weights, last_source+1, sources, pos, new_min_left_sum, right_sum, indegrees, new_indegrees1_indices);

            swap(to_do[ix], to_do[s-1]);
            pos[to_do[s-1]] = s-1;
            pos[to_do[ix]] = ix;

            order.pop_back();
            mask[x] = true;
            // in_neighbors[x] = x_in_neighbors;
            // out_neighbors[x] = x_out_neighbors;

            for (const int& triangle_id: triangles_to_reinsert){
                triangles_available[triangle_id] = true;
            }

            for (const int& v: out_neighbors[x]){
                if (mask[v]) {
                    indegrees[v] ++;
                    in_weights[v] += pair_diff[x][v];
                }
            }

            for (const int& v: in_neighbors[x]){
                if (mask[v]){
                    out_weights[v] += pair_diff[v][x];
                }
            }

            for (const int& v: order){
                right_sum[v] -= pair_diff[v][x];
            }

            // for (const int& v: out_neighbors[x]){
            //     // in_weights[v] += pair_diff[x][v] ;
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

        vector<int> out_weights(adj.size());
        for( const int& x: component){
            for (const int& v: sub_out_neighbors[x]){
                out_weights[x] += pair_crossings[x][v] - pair_crossings[v][x];
            }
        }

        // vector<int> aha;
        // for (const int& x: component){
        //     aha.push_back(in_weights[x]);
        // }
        // sort(aha.begin(), aha.end());
        // print(aha);

        

        int lb = lower_bound_mask(pair_crossings, component);
        if (verbose){
            cout << "nb unavoidable crossings: " << lb << "\n";
        }

        if (verbose){
            cout << "bad cr lower bound: " << bad_cr_lower_bound(pair_crossings, sub_in_neighbors, sub_out_neighbors, component) << endl;
        }


        // Upper bound
        // vector<int> greedy_order_natural = order_greedy_sequential_mask3(component, pair_crossings);
        // int natural_bad_cr = nb_crossings_from_order2(greedy_order_natural, pair_crossings) - lb;
        // cout << "natural bad cr: " << natural_bad_cr << "\n";
        
        auto r_natural = iterative_greedy_order(component, pair_crossings, lb);
        vector<int> greedy_order_natural = r_natural.first;
        int natural_bad_cr = r_natural.second;
        if (verbose)
        cout << "natural bad cr: " << natural_bad_cr << "\n";


        // Decreasing degree
        sort(component.begin(), component.end(), [&adj](auto a, auto b){
            return adj[a].size() > adj[b].size();
        });
        // vector<int> greedy_order_DD = order_greedy_sequential_mask3(component, pair_crossings);
        // int DD_bad_cr = nb_crossings_from_order2(greedy_order_DD, pair_crossings) - lb;
        // cout << "decreasing degree bad cr: " << DD_bad_cr << "\n";

        auto r_DD = iterative_greedy_order(component, pair_crossings, lb);
        vector<int> greedy_order_DD = r_DD.first;
        int DD_bad_cr = r_DD.second;
        if (verbose)
        cout << "decr degree bad cr: " << DD_bad_cr << "\n";

        // Increasing indegree
        sort(component.begin(), component.end(), [&sub_in_neighbors](auto a, auto b){
            return sub_in_neighbors[a].size() < sub_in_neighbors[b].size();
        });
        auto r_II = iterative_greedy_order(component, pair_crossings, lb);
        vector<int> greedy_order_incr_indegree = r_II.first;
        int incr_indegree_bad_cr = r_II.second;
        if (verbose)
        cout << "incr indegree bad cr: " << incr_indegree_bad_cr << "\n";

        // Increasing rightmost adjacency
        sort(component.begin(), component.end(), [&adj](auto a, auto b){
            return adj[a].back() < adj[b].back();
        });
        auto r_IRA = iterative_greedy_order(component, pair_crossings, lb);
        vector<int> greedy_order_IRA = r_IRA.first;
        int IRA_bad_cr = r_IRA.second;
        if (verbose)
        cout << "incr rightmost adjacency bad cr: " << IRA_bad_cr << "\n";

        // Increasing total in-weights
        sort(component.begin(), component.end(), [&in_weights, &out_weights](auto a, auto b){
            return in_weights[a]-out_weights[a] < in_weights[b]-out_weights[b];
        });
        auto r_IW = iterative_greedy_order(component, pair_crossings, lb);
        vector<int> best_order = r_IW.first;
        int best_bad_cr = r_IW.second;
        // vector<int> best_order = order_greedy_sequential_mask3(component, pair_crossings);
        // int best_bad_cr = nb_crossings_from_order2(best_order, pair_crossings) - lb;

        if (incr_indegree_bad_cr < best_bad_cr){
            best_order = greedy_order_incr_indegree;
            best_bad_cr = incr_indegree_bad_cr;
        }

        if (DD_bad_cr < best_bad_cr){
            best_order = greedy_order_DD;
            best_bad_cr = DD_bad_cr;
        }

        if (natural_bad_cr < best_bad_cr){
            best_order = greedy_order_natural;
            best_bad_cr = natural_bad_cr;
        }

        if (IRA_bad_cr < best_bad_cr){
            best_order = greedy_order_IRA;
            best_bad_cr = IRA_bad_cr;
        }

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


        // Triangles
        auto r = find_edge_disjoint_triangles_greedy3(pair_crossings, sub_in_neighbors, sub_out_neighbors, best_order);
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


        vector<int> sub_indegrees(adj.size());
        vector<int> sub_indegrees1_indices;
        for (const int& x: component){
            sub_indegrees[x] = sub_in_neighbors[x].size();
            if (sub_indegrees[x] == 1){
                sub_indegrees1_indices.push_back(pos[x]);
            }
        }
        // sort(sub_indegrees.begin(), sub_indegrees.end());
        // cout << "sub indegrees"<< endl;
        // print(sub_indegrees);
        
        vector<int> sources;
        vector<int> order;

        vector<vector<int>> pair_diff(pair_crossings.size());
        for (const int& x: component){
            pair_diff[x] = vector<int>(pair_crossings.size());
            for (const int& y: component){
                pair_diff[x][y] = pair_crossings[x][y] - pair_crossings[y][x];
            }
        }

        vector<int> branch_order = best_order;
        vector<int> min_left_sum(adj.size());
        vector<int> right_sum(adj.size());
        vector<int> right_order;

        aux4(adj, branch_order, component, s, order, right_order, best_order, best_bad_cr, current_bad_cr, pair_diff, pair_crossings, sub_in_neighbors, sub_out_neighbors, mask, 0, triangles_adj, triangles_weight, triangles_available, triangles_total,  in_weights, out_weights, adj.size(), sources, pos, min_left_sum, right_sum, sub_indegrees, sub_indegrees1_indices);

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