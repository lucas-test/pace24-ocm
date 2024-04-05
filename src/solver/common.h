#ifndef SOLVER_H
#define SOLVER_H

#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <list>

std::vector<std::vector<int>> load_file(const std::string& pathName);
void print_adj(const std::vector<std::vector<int>>& adj);
void print(const std::vector<int> v);
int lower_bound(const std::vector<std::vector<int>>& adj);
std::vector<int> greedy_sequential(const std::vector<std::vector<int>>& adj);
std::pair<int, int> crossings_between_pair(const std::vector<int>& adji, const std::vector<int>& adjj);
int nb_crossings(const std::vector<std::vector<int>>& adj, const std::vector<int>& pos);
int nb_crossings_from_order(const std::vector<std::vector<int>>& adj, const std::vector<int>& order);
void reduce_degree_0(std::vector<std::vector<int>>& adj);
std::vector<std::vector<int>> generate_random_adj(int n, int n2, double p);
void print_gr_format(std::vector<std::vector<int>> adj);
std::list<std::vector<int>> find_disjoint_3cycles(const std::vector<std::vector<int>>& adj);
std::list<std::vector<int>> find_random_disjoint_triangles(const std::vector<std::vector<int>>& adj, long m);



#endif 