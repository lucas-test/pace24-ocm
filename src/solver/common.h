#ifndef SOLVER_H
#define SOLVER_H

#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <list>

using namespace std;

vector<vector<int>> load_file(const string& pathName);
void print_adj(const vector<vector<int>>& adj);
void print(const vector<int> v);
int lower_bound(const vector<vector<int>>& adj);
int lower_bound_mask(const vector<vector<int>>& pair_crossings, const vector<int>& to_do);

vector<int> greedy_sequential(const vector<vector<int>>& adj);
pair<int, int> crossings_between_pair(const vector<int>& adji, const vector<int>& adjj);
int nb_crossings(const vector<vector<int>>& adj, const vector<int>& pos);
int nb_crossings_from_order(const vector<vector<int>>& adj, const vector<int>& order);
int nb_crossings_from_order2(const vector<int>& order, const vector<vector<int>>& pair_crossings);
vector<int> order_greedy_sequential_mask3(const vector<int>& to_do, const vector<vector<int>>& pair_crossings);

void reduce_degree_0(vector<vector<int>>& adj);
vector<vector<int>> generate_random_adj(int n, int n2, double p);
void print_gr_format(vector<vector<int>> adj);
void print_gr_format_sub(vector<vector<int>> adj, vector<int> vertices);



list<vector<int>> find_disjoint_3cycles(const vector<vector<int>>& adj);
list<vector<int>> find_random_disjoint_triangles(const vector<vector<int>>& adj, long m);
vector<int> find_disjoint_triangle(const vector<vector<int>>& pair_crossings, const vector<int>& to_do, int i, int j, const vector<vector<int>>& triangles_adj );
pair<vector<vector<int>>, vector<vector<int>>> compute_directed_graph(const vector<vector<int>>& adj);

vector<vector<int>> scc(const vector<vector<int>>& out_neighbors, const vector<vector<int>>& in_neighbors);
vector<vector<int>> scc_mask(const vector<vector<int>>& out_neighbors, const vector<vector<int>>& in_neighbors, const vector<bool>& mask);
vector<vector<int>> scc_sub_digraph(const vector<vector<int>>& out_neighbors, const vector<vector<int>>& in_neighbors, const vector<int>& sub_vertices);
vector<pair<vector<int>, bool>> scc_sub_digraph_with_sources(const vector<vector<int>>& out_neighbors, const vector<vector<int>>& in_neighbors, const vector<int>& sub_vertices) ;

bool has_duplicates(const vector<int>& numbers);

int find_disjoint_triangles(const vector<vector<int>>& adj, const vector<int>& to_do, const vector<vector<int>>& pair_crossings);
pair<int, vector<vector<vector<int>>>> find_edge_disjoint_triangles(const vector<vector<int>>& pair_crossings, const vector<vector<int>>& in_neighbors, const vector<vector<int>>& out_neighbors, const vector<int>& vertices);
int find_edge_disjoint_cycles(const vector<vector<int>>& pair_crossings, const vector<vector<int>>& in_neighbors, const vector<vector<int>>& out_neighbors, const vector<int>& vertices);
pair<int, vector<vector<vector<int>>>> find_edge_disjoint_triangles_greedy(const vector<vector<int>>& pair_crossings, const vector<vector<int>>& in_neighbors, const vector<vector<int>>& out_neighbors, const vector<int>& vertices);

pair<int, vector<vector<vector<int>>>> find_edge_disjoint_triangles_greedy2(
    const vector<bool>& mask,
    const vector<int>& vertices,
    const vector<vector<vector<int>>>& triangles);

int find_triangles_from_order(const vector<vector<int>>& pair_crossings,
    const vector<vector<int>>& in_neighbors, 
    const vector<int>& order,
    const vector<int>& vertices);

int find_triangles_from_order2(const vector<vector<int>>& pair_crossings,
    const vector<int>& order,
    const vector<int>& vertices);

int find_edge_disjoint_subgraphs(const vector<vector<int>>& pair_crossings, 
    const vector<vector<int>>& in_neighbors,
    const vector<int>& vertices);

tuple<int, vector<vector<vector<int>>>, vector<vector<bool>>> find_edge_disjoint_triangles2(
    const vector<vector<int>>& pair_crossings,
    const vector<vector<int>>& in_neighbors, 
    const vector<vector<int>>& out_neighbors,
     const vector<int>& vertices);


pair<int,int> find_triangle_replacement(const vector<vector<int>>& pair_crossings, vector<bool> excluded, int i, int j, const vector<vector<int>>& in_neighbors, const vector<vector<int>>& out_neighors, const vector<vector<vector<int>>>& triangles_adj, const vector<vector<bool>>& used_arcs );
int bad_cr_lower_bound(const vector<vector<int>>& pair_crossings, const vector<vector<int>>& in_neighbors, const vector<vector<int>>& out_neighbors, const vector<int>& vertices);

vector<vector<int>> compute_pair_crossings(
    const vector<vector<int>>& adj);

vector<vector<int>> compute_twins(const vector<vector<int>>& adj);

void to_dot(
    const vector<vector<int>>& out_neighbors,
    const vector<vector<int>>& pair_crossings
    );

bool are_equal(const vector<int>& vec1, const vector<int>& vec2);

void check_digraph_twins(const vector<int>& vertices, const vector<vector<int>>& pair_crossings);


#endif 