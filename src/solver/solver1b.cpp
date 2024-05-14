
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
 * @param adj i -> adj[i][0] is increasing, j -> adj[i][j] is increasing
 * @param to_do sorted by adj[todo[i]][0]
 * @param order 
 * @param best_order 
 * @param best_order_nc 
 */
void aux1b(const vector<vector<int>>& adj, 
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
    vector<vector<int>>& triangles_adj, int triangles_total,
    vector<bool>& excluded,
    const vector<vector<vector<int>>>& triangles ){


    if (to_do.size() == 0){
        // End of branch
        if (current_bad_cr < best_bad_cr){
            cout << "better " << current_bad_cr << "\n";
            best_order = order;
            best_bad_cr = current_bad_cr;
        }

    } else {


        

        // auto r = find_edge_disjoint_triangles(pair_crossings, in_neighbors, out_neighbors, to_do);
        // triangles_total = r.first;
        // vector<vector<vector<int>>> triangles_adj2 = r.second;

        // triangles_total =  bad_cr_lower_bound(pair_crossings, in_neighbors, out_neighbors, to_do);

        // if (r.first < triangles_total){
        //     triangles_total = r.first;
        // }

        auto r2 = find_edge_disjoint_triangles_greedy(pair_crossings, in_neighbors, out_neighbors, to_do );
        triangles_total = r2.first;
       
        
        // vector<vector<vector<int>>> triangles_adj2 = r2.second;

        // auto r4 = find_triangles_from_order(pair_crossings, in_neighbors, best_order, to_do);
        // triangles_total = r4;
        // cout << r2.first - r4 << endl;

        // auto r3 = find_edge_disjoint_triangles_greedy2(mask, to_do, triangles);
        // triangles_total = r3.first;
        // cout << r.first << " " << r2.first << endl;
        
        // cout << to_do.size() << " cur=" << current_bad_cr << " " << triangles_total << " " << r.first << " best=" << best_bad_cr << endl;

        // cout << to_do.size() << " cur=" << current_bad_cr << " " << triangles_total << " " << r2.first << " best=" << best_bad_cr << endl;
        // int weight2 = find_edge_disjoint_subgraphs(pair_crossings, in_neighbors, to_do);
        // cout << triangles_total << " " << weight2 << endl;
        // cout << triangles_total << " " << best_bad_cr << endl;

       


        // In-weights lower bound
        vector<int> t(to_do.size(), numeric_limits<int>::max());
        for (const int& v: to_do){
            vector<int> v_in_weights;
            for (const int& w: in_neighbors[v]){
                v_in_weights.push_back(pair_crossings[w][v] - pair_crossings[v][w]);
            }
            sort(v_in_weights.begin(), v_in_weights.end(), [](int a, int b) {
                return a > b;
            });
            int accu = 0;
            for (int i = v_in_weights.size()-1; i >= 0; --i){
                accu += v_in_weights[i];
                t[i] = min(t[i], accu);
            }
            for (int i = v_in_weights.size(); i < to_do.size(); ++i){
                t[i] = 0;
            }
        }

        int lower_bound0 = 0;
        for (int i = 0; i < t.size() && t[i] > 0; ++i){
            lower_bound0 += t[i];
        }
       




         // Cut with min_indegree
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
        

        // Cut with out-degrees
        vector<int> outdegrees;
        for (const int& v: to_do){
            outdegrees.push_back(out_neighbors[v].size());
        }
        sort(outdegrees.begin(), outdegrees.end());
        int rout = 0;
        for (int i = 0; i < to_do.size(); ++i){
            if (outdegrees[i] >= i){
                rout += outdegrees[i] - i;
            }
        }


        // Comparisons of the bounds
        // cout << lower_bound0 << " " << rind << " " << rout << " " << triangles_total << endl;


        // Cuttings
         if (current_bad_cr + lower_bound0 >= best_bad_cr){
            return;
        }

        if (current_bad_cr + rind >= best_bad_cr){
            // cout << "cut with indegree" << endl;
            return;
        }


        if (current_bad_cr + rout >= best_bad_cr){
            // cout << "cut with outdegrees" << endl;
            return;
        }

       
        // Cut with edge disjoint triangles
        if (current_bad_cr + triangles_total >= best_bad_cr){
            return;
        }



        // vector<vector<int>> compo = scc_sub_digraph(digraph.second, digraph.first, to_do);
        vector<pair<vector<int>,bool>> compo = scc_sub_digraph_with_sources(out_neighbors, in_neighbors, to_do);


        if (compo.size() >= 2){

            // Cut with number of components >= 3 (which have each one a cycle)
            int c3 = 0;
            for (int i = 0; i < compo.size(); ++i){
                if (compo[i].first.size() >= 3){
                    c3 ++;
                }
            }
            if (current_bad_cr + c3 >= best_bad_cr) return;

            // Launch computation on each component
            vector<int> sub_order;
            for (int i = 0; i < compo.size(); ++i){



                




                if (compo[i].first.size() == 1) {
                    if (compo[i].second && excluded[compo[i].first[0]]) return;
                    sub_order.push_back(compo[i].first[0]);
                    continue;
                }

                // if (depth <= 500 ) cout << string(depth, '-') << i << "/" << compo.size() << " size " << compo[i].first.size() <<  endl;


                // cout << string(depth, '-') << "component: " << compo[i].second ? "source": "";
                // print(compo[i].first);


                // tout recalculer
                


                // Sort for greedy insertion
                sort(compo[i].first.begin(), compo[i].first.end(), [&adj](int a, int b) {
                    return adj[a].size() > adj[b].size();
                });


                vector<bool> rmask(adj.size(), false);
                for (int j =0; j < compo[i].first.size(); ++j){
                    rmask[compo[i].first[j]] = true;
                }

               

                vector<int> rorder;


                // Recompute upper bound V0
                vector<int> rbest_order;
                for (int i = 0; i < best_order.size(); ++i){
                    if (rmask[best_order[i]]) rbest_order.push_back(best_order[i]);
                }
                int rbest_order_nc = nb_crossings_from_order2(rbest_order, pair_crossings); // compute in greedy algo
                int lb = lower_bound_mask(pair_crossings, compo[i].first);
                int rbest_bad_cr = rbest_order_nc - lb;

                // Recompute upper bound V1
                vector<int> rbest_order2 = order_greedy_sequential_mask3( compo[i].first, pair_crossings);
                int rbest_order_nc2 = nb_crossings_from_order2(rbest_order2, pair_crossings);
                if (rbest_order_nc2 < rbest_order_nc) {
                    rbest_order = rbest_order2;
                    rbest_bad_cr = rbest_order_nc2 - lb;
                } 
                
                

                int rcurrent_bad_cr = 0;

                if (rbest_bad_cr > 0){
                    

                    // Compute sub_digraph: filter neighbors to compo[i]
                    vector<vector<int>> sub_in_neighbors(in_neighbors.size());
                    for (const int& v: compo[i].first){
                        for (const int& w: in_neighbors[v]){
                            if (rmask[w]){
                                sub_in_neighbors[v].push_back(w);
                            }
                        }
                    }
                    vector<vector<int>> sub_out_neighbors(out_neighbors.size());
                    for (const int& v: compo[i].first){
                        for (const int& w: out_neighbors[v]){
                            if (rmask[w]){
                                sub_out_neighbors[v].push_back(w);
                            }
                        }
                    }
                
                   

                    // auto sub_digraph = make_pair(sub_in_neighbors, sub_out_neighbors);

                    vector<bool> next_excluded(adj.size(), false);
                    if (compo[i].second) next_excluded = excluded;

                    // auto start = chrono::high_resolution_clock::now();

                    aux1b(adj, compo[i].first, rorder, rbest_order, rbest_bad_cr, rcurrent_bad_cr, pair_crossings,sub_in_neighbors, sub_out_neighbors, rmask , depth+1, triangles_adj, 0, next_excluded, triangles);

                    // auto end = chrono::high_resolution_clock::now();
                    // chrono::duration<double> elapsed = end - start;
                    // if (elapsed.count() > 5){
                    //     cout << string(depth, '-') << i << "/" << compo.size() << " size " << compo[i].first.size() <<  endl;
                    //     cout << "duration: " << elapsed.count() << endl;
                    //     print_gr_format_sub(adj, compo[i].first);
                    // }

                    // if (depth <= 300) cout << string(depth, '-') << "final best bad cr=" << rbest_bad_cr << endl;

                    
                    
                }

                // Update global current_bad_cr
                current_bad_cr += rbest_bad_cr;

                // Cut with the remaining number of components of size >= 3
                int c3 = 0;
                for (int j = i+1; j < compo.size(); ++j){
                    if (compo[j].first.size() >= 3){
                        c3 ++;
                    }
                }
                if (current_bad_cr + c3 >= best_bad_cr){
                    return;
                }

                // Concatenate rbest_order
                sub_order.insert(sub_order.end(), rbest_order.begin(), rbest_order.end());

            }

            if (best_bad_cr > current_bad_cr){
                
                best_order = order;
                best_order.insert(best_order.end(), sub_order.begin(), sub_order.end());
                best_bad_cr = current_bad_cr;
                // cout << string(depth, '-') << "final better ";
                // print(best_order);
            }

            return;
        } else {
            // cout << compo[0].size() << endl;
        }

        // int triangles_w = find_disjoint_triangles(adj, to_do, pair_crossings);
        // if (current_bad_cr + triangles_w >= best_bad_cr){
        //     return;
        // }

        // triangles_total = find_edge_disjoint_triangles(pair_crossings, in_neighbors, out_neighbors, to_do);
        // if (current_bad_cr + triangles_total >= best_bad_cr){
        //     return;
        // }

        // sort(to_do.begin(), to_do.end(), [&in_neighbors, &out_neighbors](int a, int b) {
        //     auto aindegree = in_neighbors[out_neighbors[a][0]].size();
        //     for (const int& v: out_neighbors[a]){
        //         aindegree = min(aindegree, in_neighbors[v].size());
        //     }
        //     auto bindegree = in_neighbors[out_neighbors[b][0]].size();
        //     for (const int& v: out_neighbors[b]){
        //         bindegree = min(bindegree, in_neighbors[v].size());
        //     }

        //     return aindegree < bindegree;
        // });

        // sort(to_do.begin(), to_do.end(), [&in_neighbors, &out_neighbors, &pair_crossings](int a, int b) {
        //     auto ainw = 0;
        //     for (const int& v: in_neighbors[a]){
        //         ainw += pair_crossings[v][a] - pair_crossings[a][v];
        //     }
        //     auto aoutw = 0;
        //     for (const int& v: out_neighbors[a]){
        //         aoutw += pair_crossings[a][v] - pair_crossings[v][a];
        //     }
        //     auto binw = 0;
        //     for (const int& v: in_neighbors[b]){
        //         binw += pair_crossings[v][b] - pair_crossings[b][v];
        //     }
        //     auto boutw = 0;
        //     for (const int& v: out_neighbors[b]){
        //         boutw += pair_crossings[b][v] - pair_crossings[v][b];
        //     }

        //     return ainw-aoutw < binw-boutw;
        // });

        // vector<int> triangles_adj_w(adj.size());
        // for (const int& x: to_do){
        //     for (const auto& triangle: triangles_adj2[x]){
        //         triangles_adj_w[x] += triangle[2];
        //     }
        // }

        // sort(to_do.begin(), to_do.end(), [&triangles_adj_w, &in_neighbors, &out_neighbors, &pair_crossings](int a, int b) {
        //     auto ainw = 0;
        //     for (const int& v: in_neighbors[a]){
        //         ainw += pair_crossings[v][a] - pair_crossings[a][v];
        //     }
        //     auto binw = 0;
        //     for (const int& v: in_neighbors[b]){
        //         binw += pair_crossings[v][b] - pair_crossings[b][v];
        //     }
        //     return ainw - triangles_adj_w[a] < binw - triangles_adj_w[b];
        // });

        vector<int> in_weights(adj.size());
        for( const int& x: to_do){
            for (const int& v: in_neighbors[x]){
                in_weights[x] += pair_crossings[v][x] - pair_crossings[x][v];
            }
        }

        sort(to_do.begin(), to_do.end(), [&in_weights](int a, int b) {
            return in_weights[a] < in_weights[b];
        });

        // sort(to_do.begin(), to_do.end(), [&in_neighbors, &out_neighbors, &pair_crossings](int a, int b) {
        //     auto ainw = 0;
        //     for (const int& v: in_neighbors[a]){
        //         ainw += pair_crossings[v][a] - pair_crossings[a][v];
        //     }
        //     auto binw = 0;
        //     for (const int& v: in_neighbors[b]){
        //         binw += pair_crossings[v][b] - pair_crossings[b][v];
        //     }
        //     return ainw < binw;
        // });

        // sort(to_do.begin(), to_do.end(), [&in_neighbors](int a, int b) {
        //     return in_neighbors[a].size() < in_neighbors[b].size();
        // });


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
                int w0 = pair_crossings[x][z] - pair_crossings[z][x];
                int w1 = pair_crossings[z][y] - pair_crossings[y][z];
                int w2 = pair_crossings[y][x] - pair_crossings[x][y];
                int mw = min(w1, w2);
                if (w0 > mw) continue;
            }

            // Look for twins of x
            int nb_twins = 0;
            vector<int> twins;
            for (int j = i+1; j < to_do.size(); ++j){
                const int y = to_do[j];
                if (x != y && are_equal( adj[x], adj[y])){
                    // cout << "twin" << endl;
                    // print(adj[x]);
                    // print(adj[y]);
                    excluded[y] = true;
                    nb_twins ++;
                    twins.push_back(y);
                }
            }

           

           
            // int x_bad_cr = 0;
            // for (const int& y: in_neighbors[x]){
            //     int rfirst = pair_crossings[x][y];
            //     int rsecond = pair_crossings[y][x];
            //     x_bad_cr += rsecond- rfirst;
            // }
            int x_bad_cr = in_weights[x];
            int new_current_bad_cr = current_bad_cr + x_bad_cr*(nb_twins +1);


            // if (current_bad_cr + x_bad_cr + triangles_total - triangles_adj_w[x] >= best_bad_cr) break;

            if (current_bad_cr + x_bad_cr >= best_bad_cr) break;

            if (new_current_bad_cr  >= best_bad_cr){
                // cout << string(depth, '-') << "branch cut " << x << endl;
                continue;
            }

            
            // cout << string(depth, '-') << "branch " << x << endl;

            vector<bool> next_excluded(adj.size(), false);
            for (const int& in_neighbor: in_neighbors[x]){
                next_excluded[in_neighbor] = true;
            }
            for (const int& y: to_do){
                if (pair_crossings[x][y] == pair_crossings[y][x] && y > x){
                    next_excluded[y] = true;
                }
            }
            

            

            vector<int> in_neighbors_reinsert;
            for (const int& j: to_do) {
                
                if (pair_crossings[x][j] - pair_crossings[j][x] > 0){
                    auto it = find(in_neighbors[j].begin(), in_neighbors[j].end(), x);
                    in_neighbors[j].erase(it);
                    in_neighbors_reinsert.push_back(j);

                    for (const int& y: twins){
                        
                        if (pair_crossings[y][j] - pair_crossings[j][y] > 0){
                            auto it2 = find(in_neighbors[j].begin(), in_neighbors[j].end(), y);
                            in_neighbors[j].erase(it2);
                        }
                    }
                }
                
            }
            vector<int> out_neighbors_reinsert;
            for (const int& j: to_do) {
                if (pair_crossings[j][x] - pair_crossings[x][j] > 0){
                    auto it = find(out_neighbors[j].begin(), out_neighbors[j].end(), x);
                    out_neighbors[j].erase(it);
                    out_neighbors_reinsert.push_back(j);

                    for (const int& y: twins){
                        
                        if (pair_crossings[j][y] - pair_crossings[y][j] > 0){
                            auto it2 = find(out_neighbors[j].begin(), out_neighbors[j].end(), y);
                            out_neighbors[j].erase(it2);
                        }
                    }
                }
            }

            

            // Remove x and update all the structures
            // to_do.erase(to_do.begin() + i);
            order.push_back(x);
            mask[x] = false;

            vector<int> x_in_neighbors = in_neighbors[x];
            vector<int> x_out_neighbors = out_neighbors[x];
            in_neighbors[x] = {};
            out_neighbors[x] = {};

            vector<int> new_to_do;
            for (size_t j = 0; j < to_do.size(); ++j){
                if (i != j && find(twins.begin(), twins.end(), to_do[j]) == twins.end())
                    new_to_do.push_back(to_do[j]);
            }

            for (const int& y: twins){
                order.push_back(y);
                mask[y] = false;
                in_neighbors[y] = {};
                out_neighbors[y] = {};
            }

            aux1b(adj, new_to_do, order, best_order, best_bad_cr, new_current_bad_cr, pair_crossings, in_neighbors, out_neighbors, mask, depth+1, triangles_adj, 0, next_excluded, triangles);

            for (const int& y: twins){
                order.pop_back();
                mask[y] = true;
                in_neighbors[y] = x_in_neighbors;
                out_neighbors[y] = x_out_neighbors;
            }

            in_neighbors[x] = x_in_neighbors;
            out_neighbors[x] = x_out_neighbors;

            mask[x] = true;
            // to_do.insert(to_do.begin() + i, x);
            order.pop_back();


            for (const int& v: in_neighbors_reinsert){
                in_neighbors[v].push_back(x);
                for (const int& y: twins){
                    in_neighbors[v].push_back(y);
                }
            }
            for (const int& v: out_neighbors_reinsert){
                out_neighbors[v].push_back(x);
                for (const int& y: twins){
                    out_neighbors[v].push_back(y);
                }
            }
        }
    }
}

int solver1b(const vector<vector<int>>& adj, bool verbose) {
    if (verbose){
        cout << "############" << endl;
        cout << "solver1b c" << endl;
        cout << "nb vertices: " << adj.size() << endl;
    }
   

    vector<int> pos = greedy_sequential(adj);
    vector<int> best_order(pos.size());
    for (int i = 0; i < pos.size(); ++i) {
        best_order[pos[i]] = i ;
    }

    // best_order = {263,104,11,165,114,95,266,259,211,171,116,222,47,77,138,50,248,214,252,204,215,48,67,243,154,126,9,258,241,234,209,205,197,188,176,168,125,122,117,110,93,89,78,62,61,59,52,45,40,14,105,79,57,3,237,254,218,207,199,192,190,180,178,257,232,177,159,157,156,149,144,135,131,112,53,36,35,12,187,210,141,238,147,206,239,255,217,19,191,17,23,34,223,44,69,56,174,233,107,43,70,96,83,86,236,102,103,51,229,213,161,132,92,172,245,111,186,198,2,71,118,0,129,18,27,251,94,220,253,196,87,130,200,75,169,160,152,81,216,65,166,115,100,97,91,68,58,37,32,30,26,268,262,260,120,235,226,208,202,201,183,175,164,163,148,146,140,242,98,38,228,247,212,173,21,20,41,264,162,121,256,82,108,28,265,145,55,49,24,261,250,244,225,224,195,189,181,167,134,133,127,124,106,84,80,76,54,46,42,39,31,29,10,6,5,4,249,231,221,219,185,184,155,153,150,139,137,128,123,119,113,109,99,90,88,85,73,72,64,60,33,22,16,15,227,7,267,246,240,230,203,194,193,182,179,170,158,151,143,142,136,101,74,66,63,25,13,8,1};



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
    if (verbose){
        cout << "lower bound: " << lb << "\n";
    }

    int best_order_nc = nb_crossings_from_order(adj, best_order);
    // cout << "greedy nb crossings: " << best_order_nc << "\n";
    int best_bad_cr = best_order_nc - lb;
    //  cout << "greedy bad crossings: " << best_bad_cr << "\n";

     if (verbose){
        int bad_low1 = bad_cr_lower_bound(pair_crossings, digraph.first, digraph.second, to_do);
        cout << "bad cr lower bound: " << bad_low1 << endl;

        auto r = find_edge_disjoint_triangles(pair_crossings, digraph.first, digraph.second, to_do);
        cout << "bad cr lower bound with triangles: " << r.first << endl;

        int badcr_order = find_triangles_from_order(pair_crossings, digraph.first, best_order, to_do);
        cout << "bad cr lower bound from order: " << badcr_order << endl;
    }

    if (best_bad_cr == 0){
        if (verbose){
            cout << "min bad crossings: " << best_bad_cr << "\n";
            cout << "min crossings: " << lb + best_bad_cr << "\n";
        }   
        

        return best_bad_cr;
    }

    // Triangles
    int triangles_total = 0;

    vector<vector<int>> triangles_adj(adj.size());
    
    vector<vector<vector<int>>> triangles(adj.size());
    for (int x = 0; x < adj.size(); ++x){
        for (const int& y: digraph.first[x]){
            int yx = pair_crossings[y][x] - pair_crossings[x][y];
            for (const int& z: digraph.second[x]){
                if (pair_crossings[z][y] - pair_crossings[y][z] > 0){
                    int zy = pair_crossings[z][y] - pair_crossings[y][z];
                    int xz = pair_crossings[x][z] - pair_crossings[z][x];
                    int w = min(yx, min(zy, xz));
                    triangles[x].push_back({y,z,w});
                }
            }
        }
        sort( triangles[x].begin(), triangles[x].end(), [](vector<int> a, vector<int> b){
            return a[2] > b[2];
        });
    }


    



    vector<bool> excluded(adj.size(), false);

    int current_bad_cr = 0;

    if (verbose){
        cout << "initial best bad cr " << best_bad_cr << endl;
    }


    vector<int> order;
    aux1b(adj, to_do, order, best_order, best_bad_cr, current_bad_cr, pair_crossings, digraph.first, digraph.second, mask, 0, triangles_adj, triangles_total, excluded, triangles);

    // Print best order found
    // print(best_order);

    if (verbose){
        cout << "min bad crossings: " << best_bad_cr << "\n";
        cout << "min crossings: " << lb + best_bad_cr << "\n";
        print(best_order);
    }


    // Check solution
    // print(best_order);
    // cout << "best order size: " << best_order.size() << endl;
    // cout << "nb vertices: " << adj.size() << endl;
    // cout << "has duplicates? " << has_duplicates(best_order) << endl;
    // cout << nb_crossings_from_order(adj, best_order) << endl;
    // cout << nb_crossings_from_order(adj, best_order)-lb << endl;

    return lb + best_bad_cr;
}