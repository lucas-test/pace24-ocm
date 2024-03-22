#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <iostream>
#include <algorithm> 

using std::vector;


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
            adj[secondNumber - num1 - 1].push_back(firstNumber - 1);
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
    int nbCrossings = 0;
    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            for (const auto& x : adj[i]) {
                for (const auto& y : adj[j]) {
                    if ((pos[i] - pos[j]) * (x - y) < 0) {
                        ++nbCrossings;
                    }
                }
            }
        }
    }
    return nbCrossings;
}







/**
 * @brief 
 * 
 * @param adji indices of the neighbors of i
 * @param adjj indices of the neighbors of j
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
    int nbCrossings = 0;
    for (int i = 0; i < adj.size(); i ++){
        for (int j = i+1; j < adj.size(); j ++){
            std::pair<int,int> cr = crossings_between_pair( adj[i], adj[j]);
            nbCrossings += std::min(cr.first, cr.second);
        }
    }
    return nbCrossings;
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
