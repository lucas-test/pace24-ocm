
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
void aux2(const vector<vector<int>>& adj, 
    vector<int>& to_do, 
    vector<int>& order, 
    vector<int>& best_order, 
    int& best_bad_cr, 
    int& current_bad_cr, 
    const vector<vector<int>>& pair_crossings, 
    const pair<vector<vector<int>>,const vector<vector<int>>>& digraph, 
    vector<bool>& mask, 
    int depth, 
    vector<vector<vector<int>>>& triangles_adj, 
    int triangles_total,
    vector<vector<bool>> used_arcs,
    vector<bool> excluded ){

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



        // vector<vector<int>> compo = scc_sub_digraph(digraph.second, digraph.first, to_do);
        vector<pair<vector<int>,bool>> compo = scc_sub_digraph_with_sources(digraph.second, digraph.first, to_do);



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
                
                   

                    auto sub_digraph = make_pair(sub_in_neighbors, sub_out_neighbors);

                    vector<bool> next_excluded(adj.size(), false);
                    if (compo[i].second) next_excluded = excluded;

                    // TODO
                    int compo_triangles_total = 0;
                    for (const int& x: compo[i].first){
                        for (const auto& triangle: triangles_adj[x]){
                            compo_triangles_total += triangle[2];
                        }
                    }
                    compo_triangles_total /= 3;


                    aux2(adj, compo[i].first, rorder, rbest_order, rbest_bad_cr, rcurrent_bad_cr, pair_crossings,sub_digraph, rmask , depth+1, triangles_adj, compo_triangles_total, used_arcs, next_excluded);

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

        sort(to_do.begin(), to_do.end(), [&in_neighbors, &out_neighbors, &pair_crossings](int a, int b) {
            auto ainw = 0;
            for (const int& v: in_neighbors[a]){
                ainw += pair_crossings[v][a] - pair_crossings[a][v];
            }
            auto binw = 0;
            for (const int& v: in_neighbors[b]){
                binw += pair_crossings[v][b] - pair_crossings[b][v];
            }
            return ainw < binw;
        });

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

            // Cut if the the first 3 vertices are not
            if (order.size() >= 2){
                int w0 = pair_crossings[x][order[0]] - pair_crossings[order[0]][x];
                int w1 = pair_crossings[order[0]][order[1]] - pair_crossings[order[1]][order[0]];
                int w2 = pair_crossings[order[1]][x] - pair_crossings[x][order[1]];
                int mw = min(w1, w2);
                if (w0 > mw) continue;
            }

            // Look for twins of x
            int nb_twins = 0;
            vector<int> twins;
            for (const int& y: to_do){
                if (x != y && are_equal( adj[x], adj[y])){
                    // cout << "twin" << endl;
                    // print(adj[x]);
                    // print(adj[y]);
                    excluded[y] = true;
                    nb_twins ++;
                    twins.push_back(y);
                }
            }

           
            int x_bad_cr = 0;
            for (const int& y: in_neighbors[x]){
                int rfirst = pair_crossings[x][y];
                int rsecond = pair_crossings[y][x];
                x_bad_cr += rsecond- rfirst;
            }
            int new_current_bad_cr = current_bad_cr + x_bad_cr*(nb_twins +1);


            if (current_bad_cr + x_bad_cr >= best_bad_cr) break;

            if (new_current_bad_cr  >= best_bad_cr){
                // cout << string(depth, '-') << "branch cut " << x << endl;
                continue;
            }

            
            // cout << string(depth, '-') << "branch " << x << endl;

            vector<bool> next_excluded(adj.size(), false);

            if (nb_twins == 0){
                for (const int& in_neighbor: in_neighbors[x]){
                    next_excluded[in_neighbor] = true;
                }
                for (const int& y: to_do){
                    if (pair_crossings[x][y] == pair_crossings[y][x] && y > x){
                        next_excluded[y] = true;
                    }
                }
            } else {
                for (const int& in_neighbor: in_neighbors[x]){
                    next_excluded[in_neighbor] = true;
                }
            }
            
            // Replace triangles adjacent to x by new triangles
            int new_triangles_total = triangles_total;
            vector<vector<bool>> new_used_arcs = used_arcs;
            vector<vector<vector<int>>> new_triangles_adj(adj.size());
            for (const int& y: to_do){
                if (y != x){
                    for (const auto& triangle: triangles_adj[y]){
                        if (triangle[0] != x && triangle[1] != x){
                            new_triangles_adj[y].push_back(triangle);
                        } else {
                            int z = (triangle[0] == x) ? triangle[1] : triangle[0];
                            vector<int> t = find_triangle_replacement(pair_crossings, y, z, in_neighbors, out_neighbors, triangles_adj, used_arcs );

                            if (t.size() > 0){
                                new_triangles_adj[y].push_back({z, t[0], t[1]});
                                new_triangles_adj[z].push_back({y, t[0], t[1]});
                                new_triangles_adj[t[0]].push_back({y,z,t[1]});
                                new_triangles_total += t[1];
                                new_used_arcs[y][t[0]] = true;
                                new_used_arcs[t[0]][y] = true;
                                new_used_arcs[z][t[0]] = true;
                                new_used_arcs[t[0]][z] = true;
                            }
                        }
                    }
                }
            }
            for (const auto& triangle: triangles_adj[x]){
                int weight = triangle[2];
                new_triangles_total -= weight;
            }
            

            vector<int> in_neighbors_reinsert;
            for (const int& j: to_do) {
                auto it = find(in_neighbors[j].begin(), in_neighbors[j].end(), x);
                if (it != in_neighbors[j].end()){
                    in_neighbors[j].erase(it);
                    in_neighbors_reinsert.push_back(j);

                    for (const int& y: twins){
                        auto it2 = find(in_neighbors[j].begin(), in_neighbors[j].end(), y);
                        if (it2 != in_neighbors[j].end()){
                            in_neighbors[j].erase(it2);
                        }
                    }
                }
                
            }
            vector<int> out_neighbors_reinsert;
            for (const int& j: to_do) {
                auto it = find(out_neighbors[j].begin(), out_neighbors[j].end(), x);
                if (it != out_neighbors[j].end()){
                    out_neighbors[j].erase(it);
                    out_neighbors_reinsert.push_back(j);

                    for (const int& y: twins){
                        auto it2 = find(out_neighbors[j].begin(), out_neighbors[j].end(), y);
                        if (it2 != out_neighbors[j].end()){
                            out_neighbors[j].erase(it2);
                        }
                    }
                }
            }


       
            

            // Remove x and update all the structures
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

           

            aux2(adj, new_to_do, order, best_order, best_bad_cr, new_current_bad_cr, pair_crossings, make_pair(in_neighbors, out_neighbors), mask, depth+1, new_triangles_adj, new_triangles_total, new_used_arcs, next_excluded);

            for (const int& y: twins){
                order.pop_back();
                mask[y] = true;
                in_neighbors[y] = x_in_neighbors;
                out_neighbors[y] = x_out_neighbors;
            }

            in_neighbors[x] = x_in_neighbors;
            out_neighbors[x] = x_out_neighbors;

            mask[x] = true;
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

int solver2(const vector<vector<int>>& adj, bool verbose) {
    if (verbose){
        cout << "############" << endl;
        cout << "solver2" << endl;
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
    vector<bool> mask(adj.size(), true);


    vector<vector<int>> pair_crossings = compute_pair_crossings(adj);

    auto digraph = compute_directed_graph(adj);
    vector<pair<vector<int>,bool>> components = scc_sub_digraph_with_sources(digraph.second, digraph.first, to_do);

    if (verbose){
        cout << "nb strongly co components: " << components.size() << endl;
    }

    int lb = lower_bound_mask(pair_crossings, to_do);
    if (verbose){
        cout << "lower bound: " << lb << "\n";
    }

    int best_order_nc = nb_crossings_from_order(adj, best_order);
    int best_bad_cr = best_order_nc - lb;

    if (best_bad_cr == 0){
        if (verbose){
            cout << "min bad crossings: " << best_bad_cr << "\n";
            cout << "min crossings: " << lb + best_bad_cr << "\n";
        }   
        return best_bad_cr;
    }

    // Triangles

   

    int triangles_total;
    vector<vector<vector<int>>> triangles_adj;
    vector<vector<bool>> used_arcs;

    tie(triangles_total, triangles_adj, used_arcs) = find_edge_disjoint_triangles2(pair_crossings, digraph.first, digraph.second, to_do);


    if (verbose){
        cout << "initial triangles weight: " << triangles_total << endl;
        cout << "initial best bad cr: " << best_bad_cr << endl;
    }

    vector<bool> excluded(adj.size(), false);
    int current_bad_cr = 0;
    vector<int> order;
    aux2(adj, to_do, order, best_order, best_bad_cr, current_bad_cr, pair_crossings, digraph, mask, 0, triangles_adj, triangles_total, used_arcs, excluded);


    // Results
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