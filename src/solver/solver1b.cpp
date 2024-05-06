
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <iostream>
#include <algorithm> 
#include <list>

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
    const pair<vector<vector<int>>, 
    const vector<vector<int>>>& digraph, 
    vector<bool>& mask, 
    int deep, 
    vector<vector<int>>& triangles_adj, int triangles_total,
    const vector<bool>& excluded ){


    if (to_do.size() == 0){
        // End of branch
        if (current_bad_cr < best_bad_cr){
            cout << "better " << current_bad_cr << "\n";
            best_order = order;
            best_bad_cr = current_bad_cr;
        }

    } else {

        vector<vector<int>> in_neighbors = digraph.first;
        vector<vector<int>> out_neighbors = digraph.second;

        triangles_total = to_do.size() >= 100 ? 0 : find_edge_disjoint_triangles(pair_crossings, in_neighbors, out_neighbors, to_do);
        // triangles_total = 0;

        // Cut with edge disjoint triangles
        if (current_bad_cr + triangles_total >= best_bad_cr){
            return;
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
        if (current_bad_cr + rind >= best_bad_cr){
            // cout << "cut with indegree" << endl;
            return;
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
        if (current_bad_cr + rout >= best_bad_cr){
            // cout << "cut with outdegrees" << endl;
            return;
        }


        // int min_indegree = in_neighbors[to_do[0]].size();
        // for (const int& v: to_do){
        //     int indegree = in_neighbors[v].size();
        //     if (indegree < min_indegree){
        //         min_indegree = indegree;
        //     }
        // }
        // if (current_bad_cr + min_indegree >= best_bad_cr){
        //     return;
        // }

        // vector<vector<int>> compo = scc_sub_digraph(digraph.second, digraph.first, to_do);
        vector<pair<vector<int>,bool>> compo = scc_sub_digraph_with_sources(digraph.second, digraph.first, to_do);

        // check sources
        // for (int i = 0; i < compo.size(); ++i){
        //     bool is_source = true;
        //     for (int j = 0 ; j < compo.size(); ++j){
        //         if (i == j) continue;
        //         for (const int& v: compo[j].first){
        //             for (const int& w: compo[i].first){
        //                 auto it = find(in_neighbors[w].begin(), in_neighbors[w].end(), v);
        //                 if (it != in_neighbors[w].end()){
        //                     is_source = false;
        //                     break;
        //                 }
        //             }
        //         }
        //     }
        //     if (is_source != compo[i].second){
        //         cout << "BUGG" << i << " " << endl;
        //     }
        // }


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



                vector<bool> rmask(adj.size(), false);
                for (int j =0; j < compo[i].first.size(); ++j){
                    rmask[compo[i].first[j]] = true;
                }

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


                int triangles_weight = find_edge_disjoint_triangles(pair_crossings, sub_in_neighbors, sub_out_neighbors, compo[i].first);
                // triangles_weight = 0;
                // if (deep <= 300 ) cout << string(deep, '-') << i << "/" << compo.size() << " size " << compo[i].first.size() << " triangles= " << triangles_weight <<  endl;

                // if (i == 94){
                //     // for (int x: compo[i].first){
                //     //     cout << x << ": ";
                //     //     print(adj[x]);
                //     // }
                //     print_gr_format_sub(adj, compo[i].first);
                // }

                if (compo[i].first.size() == 1) {
                    if (compo[i].second && excluded[compo[i].first[0]]) return;
                    sub_order.push_back(compo[i].first[0]);
                    continue;
                }

                // cout << string(deep, '-') << "component: " << compo[i].second ? "source": "";
                // print(compo[i].first);


                // tout recalculer
                


                // Sort for greedy insertion
                sort(compo[i].first.begin(), compo[i].first.end(), [&adj](int a, int b) {
                    return adj[a].size() > adj[b].size();
                });


                

                vector<int> rorder;


                // Recompute upper bound V0
                vector<int> rbest_order;
                for (int i = 0; i < best_order.size(); ++i){
                    if (rmask[best_order[i]]) rbest_order.push_back(best_order[i]);
                }
                int rbest_order_nc = nb_crossings_from_order2(rbest_order, pair_crossings); // compute in greedy algo
                int lb = lower_bound_mask(adj, pair_crossings, compo[i].first);
                int rbest_bad_cr = rbest_order_nc - lb;

                // Recompute upper bound V1
                vector<int> rbest_order2 = order_greedy_sequential_mask3(adj, compo[i].first, pair_crossings);
                int rbest_order_nc2 = nb_crossings_from_order2(rbest_order2, pair_crossings);
                if (rbest_order_nc2 < rbest_order_nc) {
                    rbest_order = rbest_order2;
                    rbest_bad_cr = rbest_order_nc2 - lb;
                } 
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

                

                int rcurrent_bad_cr = 0;

                if (rbest_bad_cr > 0){
                    


                   

                    // sort(compo[i].begin(), compo[i].end(), [&sub_in_neighbors, &adj](int a, int b) {
                    //     // if (adj[a][0] == adj[b][0]){
                    //     //     return sub_in_neighbors[a].size() < sub_in_neighbors[b].size();         
                    //     // } else {
                    //     //     return adj[a][0] < adj[b][0];
                    //     // }
                    //     return adj[a][0] < adj[b][0];
                    // });

                    auto sub_digraph = make_pair(sub_in_neighbors, sub_out_neighbors);

                    int compo_triangles_total = 0;
                    // for (int j = 0; j < triangles_adj.size(); ++j){
                    //     if (rmask[j]) compo_triangles_total += triangles_adj[j][2];
                    // }
                    // compo_triangles_total /= 3;

                    // Recompute edge disjoint triangles
                    // compo_triangles_total = find_edge_disjoint_cycles(pair_crossings, sub_in_neighbors, sub_out_neighbors, compo[i].first);

                    vector<bool> next_excluded(adj.size(), false);
                    if (compo[i].second) next_excluded = excluded;

                    // if (deep <= 300) cout << string(deep, '-') << "init best bad cr=" << rbest_bad_cr << endl;
                    if (i == 94){
                        for (int j = 0; j < next_excluded.size(); ++j){
                            if (next_excluded[j]){
                                // cout << "exclude " << j << endl;
                            }
                        }
                    }

                    aux1b(adj, compo[i].first, rorder, rbest_order, rbest_bad_cr, rcurrent_bad_cr, pair_crossings,sub_digraph, rmask , deep+1, triangles_adj, 0, next_excluded);

                    // if (deep <= 300) cout << string(deep, '-') << "final best bad cr=" << rbest_bad_cr << endl;

                    // cout << string(deep, '-');
                    // print(rbest_order);
                    if (nb_crossings_from_order2(rbest_order, pair_crossings)-lb != rbest_bad_cr){
                        // cout << string(deep, '-') << "BUG final best bad cr check=" << nb_crossings_from_order2(rbest_order, pair_crossings)-lb << endl;
                    }
                    
                    
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
                // cout << string(deep, '-') << "final better ";
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

        sort(to_do.begin(), to_do.end(), [&in_neighbors, &out_neighbors](int a, int b) {
            auto aindegree = in_neighbors[out_neighbors[a][0]].size();
            for (const int& v: out_neighbors[a]){
                aindegree = min(aindegree, in_neighbors[v].size());
            }
            auto bindegree = in_neighbors[out_neighbors[b][0]].size();
            for (const int& v: out_neighbors[b]){
                bindegree = min(bindegree, in_neighbors[v].size());
            }

            return aindegree < bindegree;
        });

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

        // sort(to_do.begin(), to_do.end(), [&in_neighbors](int a, int b) {
        //     return in_neighbors[a].size() < in_neighbors[b].size();
        // });


        // Branch on to_do
        // cout << "branch " << to_do.size() << "\n";
        for (int i = 0; i < to_do.size(); ++i){
            int x = to_do[i];
            if (excluded[x]){
                // if ( 91 <= x  && x <= 136){
                //     cout << string(deep, '-') << "exclude " << x << " after "  ;
                //     cout << order.back() << " diff " << pair_crossings[x][order.back()] - pair_crossings[order.back()][x] << endl;
                //     cout << string(deep, '-');
                //     print(order);
                // }
                // cout << string(deep, '-') << "branch exclude " << x << endl;
                continue;
            } 

            // Cut with nb of in_neighbors if sorted increasing
            // if (current_bad_cr + in_neighbors[x].size()  >= best_bad_cr){ // Because triangles crossings are disjoint from the c crossings
            //     break;
            // }


           
            int new_current_bad_cr = current_bad_cr;
            for (const int& y: in_neighbors[x]){
                int rfirst = pair_crossings[x][y];
                int rsecond = pair_crossings[y][x];
                new_current_bad_cr += rsecond- rfirst;
            }

            if (new_current_bad_cr  >= best_bad_cr){
                // cout << string(deep, '-') << "branch cut " << x << endl;
                continue;
            }

            
            // if ( 91 <= x  && x <= 136)
            // cout << string(deep, '-') << "branch " << x << endl;

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
                auto it = find(in_neighbors[j].begin(), in_neighbors[j].end(), x);
                if (it != in_neighbors[j].end()){
                    in_neighbors[j].erase(it);
                    in_neighbors_reinsert.push_back(j);
                }
            }
            vector<int> out_neighbors_reinsert;
            for (const int& j: to_do) {
                auto it = find(out_neighbors[j].begin(), out_neighbors[j].end(), x);
                if (it != out_neighbors[j].end()){
                    out_neighbors[j].erase(it);
                    out_neighbors_reinsert.push_back(j);
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

            vector<int> new_to_do = to_do;
            // for (size_t j = 0; j < to_do.size(); ++j){
            //     if (i != j)
            //         new_to_do.push_back(to_do[j]);
            // }

            aux1b(adj, to_do, order, best_order, best_bad_cr, new_current_bad_cr, pair_crossings, make_pair(in_neighbors, out_neighbors), mask, deep+1, triangles_adj, 0, next_excluded);

            in_neighbors[x] = x_in_neighbors;
            out_neighbors[x] = x_out_neighbors;

            mask[x] = true;
            // to_do.insert(to_do.begin() + i, x);
            order.pop_back();


            for (int v: in_neighbors_reinsert){
                in_neighbors[v].push_back(x);
            }
            for (int v: out_neighbors_reinsert){
                out_neighbors[v].push_back(x);
            }
        }
    }
}

int solver1b(const vector<vector<int>>& adj, bool verbose) {
    if (verbose){
        cout << "############" << endl;
        cout << "solver1b" << endl;
        cout << "nb vertices: " << adj.size() << endl;
    }
   

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
    if (verbose){
        cout << "lower bound: " << lb << "\n";
    }

    int best_order_nc = nb_crossings_from_order(adj, best_order);
    // cout << "greedy nb crossings: " << best_order_nc << "\n";
    int best_bad_cr = best_order_nc - lb;
    //  cout << "greedy bad crossings: " << best_bad_cr << "\n";

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
    

    vector<bool> excluded(adj.size(), false);

    int current_bad_cr = 0;

    if (verbose){
        cout << "initial best bad cr " << best_bad_cr << endl;
    }


    vector<int> order;
    aux1b(adj, to_do, order, best_order, best_bad_cr, current_bad_cr, pair_crossings, digraph, mask, 0, triangles_adj, triangles_total, excluded);

    // Print best order found
    // print(best_order);

    if (verbose){
        cout << "min bad crossings: " << best_bad_cr << "\n";
        cout << "min crossings: " << lb + best_bad_cr << "\n";
        cout << best_order.size() << endl;
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