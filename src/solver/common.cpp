#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <iostream>
#include <algorithm> 
#include <cstdlib>
#include <ctime> 
#include <list>
#include <random>

using namespace std;
using std::vector;
using std::list;

/**
 * 
 * Notations
 * 
 * @param adj adjacencies of vertices:
 * adj[i] where i is the index of a below vertex, is the list of the fixed vertices adjacent to i
 * 
 * @param pos pos[i] is the position of vertex i (positions should be different)
 * 
 * Example:
 * adj: {0: [0]} {1: [1]} {2: [2]}
 * pos: 1 0 2
 * 
 * 0 1 2
 *  X  |
 * 1 0 2
 * 
 * Nb crossings here = 1
 * 
 */






/**
 * @brief Load a .gr file
 * 
 * @param path 
 * @return adjacencies 
 */
vector<vector<int>> load_file(const std::string& path) {
    std::string fileContent;
    std::ifstream fileStream(path);

    if (!fileStream.is_open()) {
        std::cout << "Error opening file";
        return {};
    }

    std::string line;
    std::getline(fileStream, line);
    if (line.empty()) {
        return {};
    }

    std::istringstream lineStream(line);
    std::string type, ocr;
    int num1, num2, num3;
    lineStream >> type >> ocr >> num1 >> num2 >> num3;

    vector<vector<int>> adj(num2, vector<int>());

    while (std::getline(fileStream, line)) {
        std::istringstream lineStream(line);
        int firstNumber, secondNumber;
        if (lineStream >> firstNumber >> secondNumber) {
            // adj[secondNumber - num1 - 1].push_back(firstNumber - 1);
            
            // Insert element so that adj[x] is increasing
            auto it = std::upper_bound(adj[secondNumber - num1 - 1].begin(), adj[secondNumber - num1 - 1].end(), firstNumber - 1);
            adj[secondNumber - num1 - 1].insert(it, firstNumber - 1);
        }
    }

    return adj;
}



/**
 * @brief 
 * 
 * @param adj 
 */
void print_adj(const vector<vector<int>>& adj){
    for (size_t i = 0; i < adj.size(); ++i) {
        std::cout << i << ": [";
        for (size_t j = 0; j < adj[i].size(); ++j) {
            std::cout << adj[i][j];
            if (j < adj[i].size() - 1) {
                std::cout << ",";
            }
        }
        std::cout << "]\n";
    }
}


void print(const vector<int> v){
    for (const auto& elem : v) {
        std::cout << elem << " ";
    }
    std::cout << std::endl;
}



/**
 * @brief 
 * 
 * @param adj adjacencies of vertices
 * @param pos pos[i] is the position of vertex i (positions should be different)
 * @return int 
 * @example
 * adj: {0: [0]} {1: [1]} {2: [2]}
 * pos: 1 0 2
 * 
 * 0 1 2
 *  X  |
 * 1 0 2
 * 
 * Nb crossings here = 1
 * 
 */
int nb_crossings(const vector<vector<int>>& adj, const vector<int>& pos) {
    int n = pos.size();
    int nb_crossings = 0;
    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            for (const auto& x : adj[i]) {
                for (const auto& y : adj[j]) {
                    if ((pos[i] - pos[j]) * (x - y) < 0) {
                        ++nb_crossings;
                    }
                }
            }
        }
    }
    return nb_crossings;
}


/**
 * @brief Same as nb_crossings but we compute it from the order.
 * 
 * @param adj 
 * @param order (the inverse permutation of pos)
 * @return nb of crossings.
 */
int nb_crossings_from_order(const vector<vector<int>>& adj, const vector<int>& order) {
    int n = order.size();
    int nb_crossings = 0;
    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            for (const auto& x : adj[order[i]]) {
                for (const auto& y : adj[order[j]]) {
                    if ( y < x ){
                        nb_crossings ++;
                    }
                }
            }
        }
    }
    return nb_crossings;
}







/**
 * @brief 
 * 
 * @param adji indices of the neighbors of i, supposed to be increasing
 * @param adjj indices of the neighbors of j, supposed to be increasing
 * @return [cij, cji] where cij (resp. cji) is the number of crossings between the edges adjacent to i and to j where i < j (resp. j < i)
 */
std::pair<int, int> crossings_between_pair(const std::vector<int>& adji, const std::vector<int>& adjj) {
    int jLeft = 0;
    int jRight = 0;
    int ki = 0;
    int kj = 0;
    while (ki < adji.size()) {
        while (kj < adjj.size() && adjj[kj] <= adji[ki]) {
            jLeft += ki;
            if (adjj[kj] == adji[ki]) {
                jRight += adji.size() - ki - 1;
            } else {
                jRight += adji.size() - ki;
            }
            kj++;
        }
        ki++;
    }
    while (kj < adjj.size()) {
        jLeft += ki;
        jRight += 0;
        kj++;
    }
    return std::make_pair(jLeft, jRight);
}


/**
 * @brief 
 * 
 * @param adj adjacencies
 * @return the sum of the minimum of crossings for each pair
 */
int lower_bound(const vector<vector<int>>& adj){
    int nb_crossings = 0;
    for (int i = 0; i < adj.size(); i ++){
        for (int j = i+1; j < adj.size(); j ++){
            std::pair<int,int> cr = crossings_between_pair( adj[i], adj[j]);
            nb_crossings += std::min(cr.first, cr.second);
        }
    }
    return nb_crossings;
}



/**
 * @brief insert vertex one by one following the adj vector
 * 
 * @param adj 
 * @return positions of the vertices
 */
vector<int> greedy_sequential(const vector<vector<int>>& adj) {
    int n = adj.size();
    vector<int> order;
    order.push_back(0);
    for (int i = 1; i < n; i++) {
        int cr = 0;
        int minCr = 0;
        int minK = 0;
        for (int k = 0; k < i; k++) {
            int j = order[k];
            auto [iLeft, iRight] = crossings_between_pair(adj[j], adj[i]);
            cr += iRight - iLeft;
            if (cr < minCr) {
                minCr = cr;
                minK = k + 1;
            }
        }
        order.insert(order.begin() + minK, i);
    }

    vector<int> pos(n);
    for (int i = 0; i < n; i++) {
        pos[i] = i;
    }
    for (int k = 0; k < n; k++) {
        pos[order[k]] = k;
    }

    return pos;
}




void reduce_degree_0(vector<vector<int>>& adj) {
    for (auto it = adj.begin(); it != adj.end(); ) {
        if (it->empty()) {
            it = adj.erase(it); 
        } else {
            ++it;
        }
    }
}



/**
 * @brief 
 * 
 * @param n size of the below vertices
 * @param n2 size of the fixed vertices
 * @param p proba that ij is an edge for every i in [0,n-1] and j in [0,n2-1]
 * @return vector<vector<int>> with adjacencies sorted
 */
vector<vector<int>> generate_random_adj(int n, int n2, double p){
    vector<vector<int>> adj(n);

    // Seed the random number generator
    // srand(static_cast<unsigned>(time(0)));

    // Generate edges based on probability p
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n2; ++j) {
            // Generate a random number between 0 and 1
            double randomValue = static_cast<double>(rand()) / static_cast<double>(RAND_MAX);

            // If the random number is less than p, add an edge between vertices i and j
            if (randomValue < p) {
                adj[i].push_back(j);
            }
        }
        std::sort(adj[i].begin(), adj[i].end());
    }

    return adj;
}

/**
 * @brief 
 * 
 * @param adj 
 */
void print_gr_format(vector<vector<int>> adj){
    int n1 = 0;
    int n2 = adj.size();
    int m = 0;
    for (int i = 0; i < adj.size(); i++){
        m += adj[i].size();
        for (int j = 0; j < adj[i].size(); j ++){
            if ( adj[i][j] > n1){
                n1 = adj[i][j]+1;
            }
        }
    }
    std::cout << "p ocr " << n1 << " " << n2 << " " << m << "\n";

    for (int i = 0; i < adj.size(); i ++){
        for (int j= 0; j < adj[i].size(); j ++){
            std::cout << adj[i][j]+1 << " " << i+1+n1 << "\n";
        }
    }
}


int pair_diff(const vector<int>& adji, const vector<int>& adjj){
    std::pair<int,int> r = crossings_between_pair(adji, adjj);
    return r.first - r.second;
}




/**
 * @brief Greedy search of disjoint triangles.
 * 
 * @param adj 
 * @return list<vector<int>> lists of triangles [i,j,k,w] where i,j and k are the vertices, w the weight
 */
list<vector<int>> find_disjoint_3cycles(const vector<vector<int>>& adj){
    vector<bool> used(adj.size(), false); // Set to keep track of vertices used in cycles
    // int total = 0; // sum of the minimum of the minimum weight of the cycles
    list<vector<int>> cycles;  // cycle[z] = [i,j,k,w] where i,j and k are the vertices, w the weight

    for (int i = 0; i< adj.size(); ++i){
        if (used[i]) continue;

        for (int j = i+1; j < adj.size(); ++j){
            if (used[i]) break;
            if (used[j]) continue;

            int rij = pair_diff(adj[i], adj[j]);
            for (int k = j+1; k < adj.size(); ++k){
                if (used[i]) break;
                if (used[j]) break;
                if (used[k]) continue;

                int rjk = pair_diff(adj[j], adj[k]);
                int rki = pair_diff(adj[k], adj[i]);

                if (rij > 0 && rjk > 0 && rki > 0){
                    int min = rij < rjk ? rij : rjk;
                    min = min < rki ? min : rki;
                    // std::cout << i << " " << j << " " << k  << " " << min << "\n";
                    // total += min;
                    used[i] = true;
                    used[j] = true;
                    used[k] = true;
                    cycles.push_back( {i,j,k, min});

                } else if (rij < 0 && rjk < 0 && rki < 0){
                    int max = rij > rjk ? rij : rjk;
                    max = max > rki ? max : rki;
                    max = max < 0 ? -max : max;
                    // std::cout << i << " " << j << " " << k << " " << max << "\n";
                    // total += max;
                    used[i] = true;
                    used[j] = true;
                    used[k] = true;
                    cycles.push_back( {i,j,k, max});
                }
            }
        }
    }
    return cycles;
}


/**
 * @brief 0 if not a triange
 * otherwise > 0 the min of the weights of the arcs
 * 
 * @param adj 
 * @param i 
 * @param j 
 * @param k 
 * @return int 
 */
int is_triangle(const vector<vector<int>>& adj, int i, int j, int k){
    int rij = pair_diff(adj[i], adj[j]);
    int rjk = pair_diff(adj[j], adj[k]);
    int rki = pair_diff(adj[k], adj[i]);
    if (rij > 0 && rjk > 0 && rki > 0){
        int min = rij < rjk ? rij : rjk;
        min = min < rki ? min : rki;
        return min;
    } else if (rij < 0 && rjk < 0 && rki < 0){
        int max = rij > rjk ? rij : rjk;
        max = max > rki ? max : rki;
        max = max < 0 ? -max : max;
        return max;
    }
    return 0;
}



/**
 * @brief Take randomly three vertices and check if it is a triangle.
 * Repeat m times.
 * This algorithm is very bad.
 * 
 * 
 * @param adj 
 * @param m 
 * @return list<vector<int>> 
 */
list<vector<int>> find_random_disjoint_triangles(const vector<vector<int>>& adj, long m){
    vector<int> available(adj.size());
    for (int i =0; i < adj.size(); ++i) {
        available[i] = i;
    } 
    list<vector<int>> triangles;  
    random_device rd;
    mt19937 gen(rd());

    for (long i = 0; i < m; ++i) {
        int r1 = gen() % adj.size();
        int r2 = gen() % adj.size();
        int r3 = gen() % adj.size();
        if (r2 == r1 || r3 == r1 || r3 == r2) continue;

        if (r1 > r2) {
            swap(r1, r2);
        }
        if (r1 > r3) {
            swap(r1, r3);
        }
        if (r2 > r3) {
            swap(r2, r3);
        }

        int v1 = available[r1];
        int v2 = available[r2];
        int v3 = available[r3];

        int r = is_triangle(adj, v1, v2, v3);
        // cout << v1 << " " << v2 << " " << v3 << " " << r << endl;
        if (r > 0){
            triangles.push_back({v1, v2, v3, r});
            available.erase(available.begin() + r3);
            available.erase(available.begin() + r2);
            available.erase(available.begin() + r1);
        }
    }
    
    return triangles;
}