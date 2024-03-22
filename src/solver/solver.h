#ifndef SOLVER_H
#define SOLVER_H

#include <fstream>
#include <sstream>
#include <vector>
#include <string>

std::vector<std::vector<int>> load_file(const std::string& pathName);
void print_adj(const std::vector<std::vector<int>>& adj);
int lower_bound(const std::vector<std::vector<int>>& adj);
std::vector<int> greedy_sequential(const std::vector<std::vector<int>>& adj);
std::pair<int, int> crossings_between_pair(const std::vector<int>& adji, const std::vector<int>& adjj);
int nb_crossings(const std::vector<std::vector<int>>& adj, const std::vector<int>& pos);



#endif 