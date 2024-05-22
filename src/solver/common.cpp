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
#include <stack>
#include <unordered_set>

#include "solver_bruteforce.h"


using namespace std;

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


bool are_equal(const vector<int>& vec1, const vector<int>& vec2) {
    if (vec1.size()!= vec2.size()) {
        return false;
    }
    return equal(vec1.begin(), vec1.end(), vec2.begin());
}



/**
 * @brief Load a .gr file
 * 
 * @param path 
 * @return adjacencies 
 */
pair<vector<vector<int>>, int> load_file(const std::string& path) {
    std::string fileContent;
    std::ifstream fileStream(path);

    if (!fileStream.is_open()) {
        std::cout << "Error opening file";
        return {};
    }

    string line;

    getline(fileStream, line);
    while (!line.empty() && line[0] == 'c') {
        getline(fileStream, line);
    }

    std::istringstream lineStream(line);
    std::string type, ocr;
    int num1, num2, num3;
    lineStream >> type >> ocr >> num1 >> num2 >> num3;

    vector<vector<int>> adj(num2, vector<int>());
    // vector<int> corresp(adj.size());

    while (std::getline(fileStream, line)) {
        if (!line.empty() && line[0] == 'c'){
            continue;
        }

        std::istringstream lineStream(line);
        int firstNumber, secondNumber;
        if (lineStream >> firstNumber >> secondNumber) {
            // adj[secondNumber - num1 - 1].push_back(firstNumber - 1);
            
            // Insert element so that adj[x] is increasing
            auto it = std::upper_bound(adj[secondNumber - num1 - 1].begin(), adj[secondNumber - num1 - 1].end(), firstNumber - 1);
            adj[secondNumber - num1 - 1].insert(it, firstNumber - 1);
            // corresp[secondNumber - num1 - 1] = secondNumber;
        }
    }

    return make_pair(adj, num1);
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
        sort(adj[i].begin(), adj[i].end());
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

void print_gr_format_sub(vector<vector<int>> adj, vector<int> vertices){
    int n1 = 0;
    int n2 = vertices.size();
    int m = 0;
    for (const int& i: vertices){
        m += adj[i].size();
        for (int j = 0; j < adj[i].size(); j ++){
            if ( adj[i][j] > n1){
                n1 = adj[i][j]+1;
            }
        }
    }
    std::cout << "p ocr " << n1 << " " << n2 << " " << m << "\n";

    for (int i = 0; i < vertices.size(); i ++){
        int x = vertices[i];
        for (int j= 0; j < adj[x].size(); j ++){
            std::cout << adj[x][j]+1 << " " << i+1+n1 << "\n";
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


vector<int> find_disjoint_triangle(const vector<vector<int>>& pair_crossings, const vector<int>& to_do, int i, int j, const vector<vector<int>>& triangles_adj ){
    if (pair_crossings[i][j] - pair_crossings[j][i] < 0){
        swap(i,j);
    }
    int ij = pair_crossings[i][j] - pair_crossings[j][i]; // > 0
    
    int bestk = -1;
    int bestw = 0;
    for (int u = 0; u < to_do.size(); ++u){
        int k = to_do[u];
        if (triangles_adj[k][2] > 0) continue;
        if (k == i || k == j) continue;

        int ki = pair_crossings[k][i] - pair_crossings[i][k];
        int jk = pair_crossings[j][k] - pair_crossings[k][j];
        if (ki > 0 && jk > 0){
            int minw = min({ij, jk, ki});
            if (minw > bestw){
                bestw = minw;
                bestk = k;
                if (bestw == ij) break;
            }
        }

    }
    return {bestk, bestw};
}


pair<int,int> find_triangle_replacement(const vector<vector<int>>& pair_crossings, vector<bool> excluded, int i, int j, const vector<vector<int>>& in_neighbors, const vector<vector<int>>& out_neighors, const vector<vector<vector<int>>>& triangles_adj, const vector<vector<bool>>& used_arcs ){
    if (pair_crossings[i][j] - pair_crossings[j][i] < 0){
        swap(i,j);
    }
    int ij = pair_crossings[i][j] - pair_crossings[j][i]; // > 0
    
    int bestk = -1;
    int bestw = 0;
    for (const int& z: in_neighbors[i]){
        if (excluded[z]) continue;
        if (used_arcs[i][z] == false && used_arcs[j][z] == false){
            if (pair_crossings[j][z] - pair_crossings[z][j] > 0){
                int zi = pair_crossings[z][i] - pair_crossings[i][z];
                int jz = pair_crossings[j][z] - pair_crossings[z][j];
                int w = min({ij, jz, zi});
                return make_pair(z, w );
            }
        }
    }
    return make_pair(-1,-1);
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



/**
 * @brief 
 * 
 * @param adj 
 * @return <in_neighbors, out_neighbors>
 */
pair<vector<vector<int>>, vector<vector<int>>> compute_directed_graph(const vector<vector<int>>& adj){
    vector<vector<int>> in_neighbors(adj.size());
    vector<vector<int>> out_neighbors(adj.size());
    for (int i = 0; i < adj.size(); ++i){
        for (int j = 0; j < i; ++j){
            if (i==j) continue;
            int wij = pair_diff(adj[i], adj[j]);
            // cout << i << " "<< j << " " << wij << endl;
            if (wij > 0){
                out_neighbors[i].push_back(j);
                in_neighbors[j].push_back(i);
            } else if (wij < 0){
                in_neighbors[i].push_back(j);
                out_neighbors[j].push_back(i);
            }
        }
    }
    return make_pair(in_neighbors, out_neighbors);
}



void visit (int cur, unordered_set<int>& visited,const vector<vector<int>>& out_neighbors, stack<int>& stack) {
    if (visited.count(cur)) return;
    visited.insert(cur);
    for (const auto& neigh : out_neighbors[cur]) {
        visit(neigh, visited, out_neighbors, stack);
    }
    stack.push(cur);
}

void assign (int cur, const vector<vector<int>>& in_neighbors, unordered_set<int>& assigned, vector<vector<int>>& scc) {
    if (!assigned.count(cur)) {
        assigned.insert(cur);
        vector<int> rootStack;
        if (!scc.empty()) {
            rootStack = scc.back();
            scc.pop_back();
        }
        rootStack.push_back(cur); 
        scc.push_back(rootStack); // scc.back().push_back(cur) // en espérant que scc.back() soit pas une copie
        for (const auto& neigh : in_neighbors[cur]) {
            assign(neigh, in_neighbors, assigned, scc);
        }
    }
}


vector<vector<int>> scc(const vector<vector<int>>& out_neighbors, const vector<vector<int>>& in_neighbors) {
    vector<vector<int>> scc; // Strongly Connected Components
    stack<int> stack;
    unordered_set<int> visited;

    for (int i = 0; i < out_neighbors.size(); ++i) {
        visit(i, visited, out_neighbors, stack);
    }

    unordered_set<int> assigned;

    while (!stack.empty()) {
        int stack_head = stack.top();
        stack.pop();
        if (!assigned.count(stack_head)) {
            scc.push_back(vector<int>()); // The array to stock the new component
            assign(stack_head, in_neighbors, assigned, scc);
        }
    }

    return scc;
    }


/**
 * @brief 
 * 
 * @param cur 
 * @param in_neighbors 
 * @param assigned 
 * @param component 
 * @return true if the component is a source
 * @return false otherwise (there is an arc from a vertex not in the component to this component)
 */
bool assign2 (int cur, const vector<vector<int>>& in_neighbors, vector<bool>& assigned, vector<int>& component) {
    if (assigned[cur] == false) {
        assigned[cur] = true;
        component.push_back(cur);
        bool is_source = true;
        for (const auto& neigh : in_neighbors[cur]) {
            if (assigned[neigh] && find(component.begin(), component.end(), neigh) == component.end()){
                is_source = false;
            } else {
                if (assign2(neigh, in_neighbors, assigned, component) == false){
                    is_source = false;
                }
            }
        }
        return is_source;
    }
    return true;
}


vector<pair<vector<int>, bool>> scc_sub_digraph_with_sources(const vector<vector<int>>& out_neighbors, const vector<vector<int>>& in_neighbors, const vector<int>& sub_vertices) {
    vector<pair<vector<int>,bool>> scc; // Strongly Connected Components
    stack<int> stack;
    unordered_set<int> visited;

    for (const int& i: sub_vertices) {
        visit(i, visited, out_neighbors, stack);
    }

    vector<bool> assigned(out_neighbors.size(), false);

    while (!stack.empty()) {
        int stack_head = stack.top();
        stack.pop();
        if (assigned[stack_head] == false) {
            vector<int> component;
            bool is_source = assign2(stack_head, in_neighbors, assigned, component);
            scc.push_back(make_pair(component, is_source));
        }
    }

    return scc;
}




void visit_mask (int cur, unordered_set<int>& visited,const vector<vector<int>>& out_neighbors, stack<int>& stack, const vector<bool>& mask) {
    if (visited.count(cur)) return;
    visited.insert(cur);
    for (const auto& neigh : out_neighbors[cur]) {
        if (mask[neigh])
            visit_mask(neigh, visited, out_neighbors, stack, mask);
    }
    stack.push(cur);
}

void assign_mask (int cur, const vector<vector<int>>& in_neighbors, unordered_set<int>& assigned, vector<vector<int>>& scc, const vector<bool>& mask) {
    if (!assigned.count(cur)) {
        assigned.insert(cur);
        vector<int> rootStack;
        if (!scc.empty()) {
            rootStack = scc.back();
            scc.pop_back();
        }
        rootStack.push_back(cur);
        scc.push_back(rootStack);
        for (const auto& neigh : in_neighbors[cur]) {
            if (mask[neigh])
                assign_mask(neigh, in_neighbors, assigned, scc, mask);
        }
    }
}


vector<vector<int>> scc_mask(const vector<vector<int>>& out_neighbors, const vector<vector<int>>& in_neighbors, const vector<bool>& mask) {
    vector<vector<int>> scc; // Strongly Connected Components
    stack<int> stack;
    unordered_set<int> visited;

    for (int i = 0; i < out_neighbors.size(); ++i) {
        if (mask[i])  visit_mask(i, visited, out_neighbors, stack, mask);
    }

    unordered_set<int> assigned;

    while (!stack.empty()) {
        int stack_head = stack.top();
        stack.pop();
        if (!assigned.count(stack_head)) {
            scc.push_back(vector<int>()); // The array to stock the new component
            assign_mask(stack_head, in_neighbors, assigned, scc, mask);
        }
    }

    return scc;
}


vector<vector<int>> scc_sub_digraph(const vector<vector<int>>& out_neighbors, const vector<vector<int>>& in_neighbors, const vector<int>& sub_vertices) {
    vector<vector<int>> scc; // Strongly Connected Components
    stack<int> stack;
    unordered_set<int> visited;

    for (int i: sub_vertices) {
        visit(i, visited, out_neighbors, stack);
    }

    unordered_set<int> assigned;

    while (!stack.empty()) {
        int stack_head = stack.top();
        stack.pop();
        if (!assigned.count(stack_head)) {
            scc.push_back(vector<int>()); // The array to stock the new component
            assign(stack_head, in_neighbors, assigned, scc);
        }
    }

    return scc;
}







/**
 * @brief check for duplicates in a vector.
 * Just used to check if the vectors in the algorithm are ok.
 * 
 * @param numbers 
 * @return true if has duplicates
 * @return false 
 */
bool has_duplicates(const vector<int>& numbers) {
    unordered_set<int> seen;
    for (int number : numbers) {
        if (seen.count(number) > 0) {
            return true; // Duplicate found
        }
        seen.insert(number);
    }
    return false; // No duplicates found
}





int lower_bound_mask(const vector<vector<int>>& pair_crossings, const vector<int>& to_do){
    int nb_crossings = 0;
    for (int i = 0; i < to_do.size(); i ++){
        for (int j = i+1; j < to_do.size(); j ++){
            nb_crossings += min(pair_crossings[to_do[i]][to_do[j]], pair_crossings[to_do[j]][to_do[i]]);
        }
    }
    return nb_crossings;
}

int nb_crossings_from_order2(const vector<int>& order, const vector<vector<int>>& pair_crossings) {
    int n = order.size();
    int nb_crossings = 0;
    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            nb_crossings += pair_crossings[order[j]][order[i]];
        }
    }
    return nb_crossings;
}




int find_disjoint_triangles(const vector<vector<int>>& adj, const vector<int>& to_do, const vector<vector<int>>& pair_crossings){
    vector<bool> used(adj.size(), false); // Set to keep track of vertices used in cycles
    int total = 0; // sum of the minimum of the minimum weight of the cycles

    for (int i = 0; i< to_do.size(); ++i){
        int x = to_do[i];
        if (used[x]) continue;

        for (int j = i+1; j < to_do.size(); ++j){
            int y = to_do[j];
            if (used[x]) break;
            if (used[y]) continue;

            int rij = pair_crossings[x][y] - pair_crossings[y][x];
            for (int k = j+1; k < to_do.size(); ++k){
                int z = to_do[k];
                if (used[x]) break;
                if (used[y]) break;
                if (used[z]) continue;

                int rjk = pair_crossings[y][z] - pair_crossings[z][y];
                int rki = pair_crossings[z][x] - pair_crossings[x][z];

                if (rij > 0 && rjk > 0 && rki > 0){
                    int min = rij < rjk ? rij : rjk;
                    min = min < rki ? min : rki;
                    // std::cout << i << " " << j << " " << k  << " " << min << "\n";
                    total += min;
                    used[x] = true;
                    used[y] = true;
                    used[z] = true;
                } else if (rij < 0 && rjk < 0 && rki < 0){
                    int max = rij > rjk ? rij : rjk;
                    max = max > rki ? max : rki;
                    max = max < 0 ? -max : max;
                    total += max;
                    used[x] = true;
                    used[y] = true;
                    used[z] = true;
                }
            }
        }
    }
    return total;
}




vector<int> order_greedy_sequential_mask3(const vector<int>& to_do, const vector<vector<int>>& pair_crossings) {
    vector<int> order;
    for (int iId =0; iId < to_do.size(); ++iId) {
        int i = to_do[iId];
        int cr = 0;
        int minCr = 0;
        int minK = 0;
        for (int k = 0; k < order.size(); ++k) {
            int j = order[k];
            int iLeft = pair_crossings[j][i];
            int iRight = pair_crossings[i][j];
            cr += iRight - iLeft;
            if (cr < minCr) {
                minCr = cr;
                minK = k + 1;
            }
        }
        order.insert(order.begin() + minK, i);
    }
    return order;
}







pair<int, vector<vector<vector<int>>>> find_edge_disjoint_triangles(const vector<vector<int>>& pair_crossings, const vector<vector<int>>& in_neighbors, const vector<vector<int>>& out_neighbors, const vector<int>& vertices){
    vector<vector<bool>> used(pair_crossings.size());
    for (int i= 0; i < used.size(); ++i){
        used[i] = vector<bool>(pair_crossings.size(), false);
    }
    int total = 0;
    // vector<vector<int>> triangles;
    vector<vector<vector<int>>> triangles_adj(pair_crossings.size());

    for (const int& x: vertices){

        vector<int> x_in_neighbors = in_neighbors[x];
        sort(x_in_neighbors.begin(), x_in_neighbors.end(), [&pair_crossings, &x](int a, int b) {
            return pair_crossings[a][x] - pair_crossings[x][a] > pair_crossings[b][x] - pair_crossings[x][b];
        });

        for (const int& y: x_in_neighbors){
            if (used[x][y]) continue;
            int yx = pair_crossings[y][x] - pair_crossings[x][y];

            vector<int> x_out_neighbors = out_neighbors[x];
            sort(x_out_neighbors.begin(), x_out_neighbors.end(), [&pair_crossings, &x](int a, int b) {
                return pair_crossings[x][a] - pair_crossings[a][x] > pair_crossings[x][b] - pair_crossings[b][x];
            });

            for (const int& z: x_out_neighbors){
                if (used[x][y] || used[x][z] || used[y][z]) continue;
                int zy = pair_crossings[z][y] - pair_crossings[y][z];
                if (zy > 0){

                    // bool found = false;
                    // for (const int& w: out_neighbors[z]){
                    //     if (w == y || used[w][z] || used[w][x] || used[x][y] || used[x][z] || used[y][z] ) continue;
                    //     auto it = find(in_neighbors[x].begin(), in_neighbors[x].end(), w);
                    //     if (it != in_neighbors[x].end()){
                    //         // z -> w -> x
                    //         int xz = pair_crossings[x][z] - pair_crossings[z][x];
                    //         int zy = pair_crossings[z][y] - pair_crossings[y][z];
                    //         int zw = pair_crossings[z][w] - pair_crossings[w][z];
                    //         int wx = pair_crossings[w][x] - pair_crossings[x][w];
                    //         int weight = min(xz, min(min(yx, zy), min(zw, wx)));
                    //         if (weight <=  min(min(yx, xz),zy)) continue;
                            
                    //         found = true;

                    //         used[x][y] = true;
                    //         used[y][x] = true;
                    //         used[x][z] = true;
                    //         used[z][x] = true;
                    //         used[z][y] = true;
                    //         used[y][z] = true;
                    //         used[z][w] = true;
                    //         used[w][z] = true;
                    //         used[x][w] = true;
                    //         used[w][x] = true;

                    //         triangles.push_back({x,z,y});
                    //         total += weight;
                    //         break;
                    //     }
                    // }
                    // if (found) break;

                    // triangle y > x -> z -> y
                    used[x][y] = true;
                    used[y][x] = true;
                    used[x][z] = true;
                    used[z][x] = true;
                    used[z][y] = true;
                    used[y][z] = true;
                    

                    
                    int xz = pair_crossings[x][z] - pair_crossings[z][x];
                    int weight = min(min(yx, xz),zy);
                    // triangles.push_back({x,z,y});
                    triangles_adj[x].push_back({z,y,weight});
                    triangles_adj[y].push_back({x,z,weight});
                    triangles_adj[z].push_back({x,y,weight});

                    total += weight;
                    break;
                }

                // Search for 4-cycles
                // if (used[x][y]) break;
                // if (used[x][z]) continue;
                
                // for (const int& w: out_neighbors[z]){
                //     if (used[x][y] || used[x][z] || used[z][w] || used[y][w]
                //     || used[y][x] || used[z][x] ||used[w][z] || used[w][y]) continue;
                //     auto it = find(out_neighbors[w].begin(), out_neighbors[w].end(), y);
                //     if (it != out_neighbors[w].end() ){
                //         used[x][y] = true;
                //         used[y][x] = true;
                //         used[x][z] = true;
                //         used[z][x] = true;
                //         used[z][w] = true;
                //         used[w][z] = true;
                //         used[w][y] = true;
                //         used[y][w] = true;

                //         int yx = pair_crossings[y][x] - pair_crossings[x][y];
                //         int xz = pair_crossings[x][z] - pair_crossings[z][x];
                //         int zw = pair_crossings[z][w] - pair_crossings[w][z];
                //         int wy = pair_crossings[w][y] - pair_crossings[y][w];
                //         int weight = min(min(yx, xz),min(zw, wy));
                //         total += weight;
                //         break;
                //     }
                // }
            }
        }
    }
    return make_pair(total, triangles_adj);
}


int find_edge_disjoint_cycles(const vector<vector<int>>& pair_crossings, const vector<vector<int>>& in_neighbors, const vector<vector<int>>& out_neighbors, const vector<int>& vertices){
    vector<vector<bool>> used(in_neighbors.size());
    for (int i= 0; i < used.size(); ++i){
        used[i] = vector<bool>(in_neighbors.size(), false);
    }
    int total = 0;
    vector<vector<int>> triangles;

    for (const int& x: vertices){
        for (const int& y: in_neighbors[x]){
            if (used[x][y]) continue;
            for (const int& z: out_neighbors[x]){
                if (used[x][y] || used[x][z] || used[y][z]) continue;
                auto it = find(out_neighbors[z].begin(), out_neighbors[z].end(), y);
                if (it != out_neighbors[z].end() ){
                    // triangle (xzy)
                    used[x][y] = true;
                    used[y][x] = true;
                    used[x][z] = true;
                    used[z][x] = true;
                    used[z][y] = true;
                    used[y][z] = true;
                    
                    triangles.push_back({x,z,y});
                    int yx = pair_crossings[y][x] - pair_crossings[x][y];
                    int xz = pair_crossings[x][z] - pair_crossings[z][x];
                    int zy = pair_crossings[z][y] - pair_crossings[y][z];
                    int weight = min(min(yx, xz),zy);
                    total += weight;
                    break;
                }

                // Search for 4-cycles
                if (used[x][y]) break;
                
                if (used[x][z]) continue;
                
                for (const int& w: out_neighbors[z]){
                    if (used[x][y] || used[x][z] || used[z][w] || used[y][w]
                    || used[y][x] || used[z][x] ||used[w][z] || used[w][y]) continue;
                    auto it = find(out_neighbors[w].begin(), out_neighbors[w].end(), y);
                    if (it != out_neighbors[w].end() ){
                        used[x][y] = true;
                        used[y][x] = true;
                        used[x][z] = true;
                        used[z][x] = true;
                        used[z][w] = true;
                        used[w][z] = true;
                        used[w][y] = true;
                        used[y][w] = true;

                        triangles.push_back({x,z,w,y});
                        int yx = pair_crossings[y][x] - pair_crossings[x][y];
                        int xz = pair_crossings[x][z] - pair_crossings[z][x];
                        int zw = pair_crossings[z][w] - pair_crossings[w][z];
                        int wy = pair_crossings[w][y] - pair_crossings[y][w];
                        int weight = min(min(yx, xz),min(zw, wy));
                        total += weight;
                        break;
                    }
                }

            }
        }
    }
    // cout << "########" << endl;
    // for (const vector<int>& triangle: triangles){
    //     print(triangle);
    // }

    // Check if cycles are disjoint
    // for (int i = 0; i < triangles.size(); ++i){
    //     for (int j = 0; j < i; ++j){
    //         for (int k = 0; k < triangles[i].size(); ++k){
    //             for (int l = 0; l < triangles[j].size(); ++l){
    //                 int kk = k+1 == triangles[i].size() ? 0 : k+1;
    //                 int ll = l+1 == triangles[j].size() ? 0 : l+1;
    //                 if (triangles[i][k] == triangles[j][l] && triangles[i][kk] == triangles[j][ll]){
    //                     cout << "bug"<< endl;
    //                 }
    //             }
    //         }
    //     }
    // }
    // cout << "triangles: " << triangles.size() << " " << total << endl;
    return total;
}



int find_edge_disjoint_subgraphs(const vector<vector<int>>& pair_crossings, 
    const vector<vector<int>>& in_neighbors,
    const vector<int>& vertices){

        vector<int> vertices2 = vertices;

        sort(vertices2.begin(), vertices2.end(), [&in_neighbors](int x, int y) {
            return in_neighbors[x].size() < in_neighbors[y].size(); 
        });

        int weight = 0;

        for (int i = 0; i < vertices.size()/6; ++i){
            vector<int> to_insert;
            vector<int> order;
            for (int j = 0; j < 6; ++j){
                to_insert.push_back(vertices2[6*i+j]);
            }
            int lb = lower_bound_mask(pair_crossings, to_insert);
            weight += auxBF(to_insert, order, pair_crossings, false) -lb;
        }
        return weight;
    }



tuple<int, vector<vector<vector<int>>>, vector<vector<bool>>> find_edge_disjoint_triangles2(
    const vector<vector<int>>& pair_crossings,
    const vector<vector<int>>& in_neighbors, 
    const vector<vector<int>>& out_neighbors,
    const vector<int>& vertices){

    vector<vector<bool>> used(pair_crossings.size());
    for (int i= 0; i < used.size(); ++i){
        used[i] = vector<bool>(pair_crossings.size(), false);
    }
    int total = 0;
    vector<vector<vector<int>>> triangles_adj(pair_crossings.size());

    for (const int& x: vertices){

        vector<int> x_in_neighbors = in_neighbors[x];
        sort(x_in_neighbors.begin(), x_in_neighbors.end(), [&pair_crossings, &x](int a, int b) {
            return pair_crossings[a][x] - pair_crossings[x][a] > pair_crossings[b][x] - pair_crossings[x][b];
        });

        for (const int& y: x_in_neighbors){
            if (used[x][y]) continue;
            int yx = pair_crossings[y][x] - pair_crossings[x][y];

            vector<int> x_out_neighbors = out_neighbors[x];
            sort(x_out_neighbors.begin(), x_out_neighbors.end(), [&pair_crossings, &x](int a, int b) {
                return pair_crossings[x][a] - pair_crossings[a][x] > pair_crossings[x][b] - pair_crossings[b][x];
            });

            for (const int& z: x_out_neighbors){
                if (used[x][y] || used[x][z] || used[y][z]) continue;
                auto it = find(out_neighbors[z].begin(), out_neighbors[z].end(), y);
                if (it != out_neighbors[z].end() ){

                    // triangle y > x -> z -> y
                    used[x][y] = true;
                    used[y][x] = true;
                    used[x][z] = true;
                    used[z][x] = true;
                    used[z][y] = true;
                    used[y][z] = true;
                    
                    int xz = pair_crossings[x][z] - pair_crossings[z][x];
                    int zy = pair_crossings[z][y] - pair_crossings[y][z];
                    int weight = min(min(yx, xz),zy);
                    triangles_adj[x].push_back({z,y,weight});
                    triangles_adj[y].push_back({x,z,weight});
                    triangles_adj[z].push_back({x,y,weight});
                    
                    
                    total += weight;
                    break;
                }

              
            }
        }
    }
    return make_tuple(total, triangles_adj, used);
}



vector<vector<int>> compute_pair_crossings(
    const vector<vector<int>>& adj){
        vector<vector<int>> pair_crossings(adj.size());
        for (int i = 0; i < adj.size(); ++i){
            vector<int> crossings_i(adj.size());
            pair_crossings[i] = crossings_i;
            for (int j = 0; j < adj.size(); ++j){
                pair_crossings[i][j] = crossings_between_pair(adj[i], adj[j]).first;
            }
        }
        return pair_crossings;
    }



/**
 * x and y are twins if they have the same adjacencies
*/
vector<vector<int>> compute_twins(const vector<vector<int>>& adj){
    vector<vector<int>> twins(adj.size());
    for (int i = 0; i < adj.size(); ++i){
        for (int j = 0; j < adj.size(); ++j){
            if (are_equal( adj[i], adj[j])){
                twins[i].push_back(j);
            }
        }
    }
    return twins;
}


void to_dot(
    const vector<vector<int>>& out_neighbors,
    const vector<vector<int>>& pair_crossings
    ) {
        cout << "digraph G {" << endl;
        for (int i = 0; i < out_neighbors.size(); ++i){
            for (const int& j: out_neighbors[i]){
                printf("%d -> %d [label=%d];\n", i, j, pair_crossings[i][j] - pair_crossings[j][i]);
            }
        }
        cout << "}" << endl;
    }








pair<int, vector<vector<vector<int>>>> find_edge_disjoint_triangles_greedy(const vector<vector<int>>& pair_crossings, const vector<vector<int>>& in_neighbors, const vector<vector<int>>& out_neighbors, const vector<int>& vertices){
    vector<vector<bool>> used(pair_crossings.size());
    for (int i= 0; i < used.size(); ++i){
        used[i] = vector<bool>(pair_crossings.size(), false);
    }
    int total = 0;
    vector<vector<vector<int>>> triangles_adj(pair_crossings.size());



    for (const int& x: vertices){

        vector<int> x_in_neighbors = in_neighbors[x];
        sort(x_in_neighbors.begin(), x_in_neighbors.end(), [&pair_crossings, &x](int a, int b) {
            return pair_crossings[a][x] - pair_crossings[x][a] > pair_crossings[b][x] - pair_crossings[x][b];
        });

        for (const int& y: x_in_neighbors){
            if (used[x][y]) continue;
            int yx = pair_crossings[y][x] - pair_crossings[x][y];

            int bestz = -1;
            int best_w = 0;

            for (const int& z: out_neighbors[x]){
                if (used[x][y] || used[x][z] || used[y][z]) continue;
                int zy = pair_crossings[z][y] - pair_crossings[y][z];
                if (zy > 0 ){
                    int xz = pair_crossings[x][z] - pair_crossings[z][x];
                    int weight = min(min(yx, xz),zy);

                    if (weight > best_w){
                        bestz = z;
                        best_w = weight;
                    }
                }
            }

            if (bestz >= 0){
                int z = bestz;
                int weight = best_w;
                // triangle y > x -> z -> y
                used[x][y] = true;
                used[y][x] = true;
                used[x][z] = true;
                used[z][x] = true;
                used[z][y] = true;
                used[y][z] = true;
                // triangles.push_back({x,z,y});
                triangles_adj[x].push_back({z,y,weight});
                triangles_adj[y].push_back({x,z,weight});
                triangles_adj[z].push_back({x,y,weight});

                total += weight;
            }
            
        }
    }
    return make_pair(total, triangles_adj);
}

/**
 * Same algo as greedy1 but returns a different structure.
 * vector<vector<int>> triangles_adj[x] contains the array of the triangle_ids adjacent to x
 * vector<int> triangles_weight[tid] is the weight of triangle having id = tid
*/
pair<vector<vector<int>>, vector<int>> find_edge_disjoint_triangles_greedy3(const vector<vector<int>>& pair_crossings, const vector<vector<int>>& in_neighbors, const vector<vector<int>>& out_neighbors, const vector<int>& vertices){
    vector<vector<bool>> used(pair_crossings.size());
    for (int i= 0; i < used.size(); ++i){
        used[i] = vector<bool>(pair_crossings.size(), false);
    }
    vector<vector<int>> triangles_adj(pair_crossings.size());
    vector<int> triangles_weight;
    int triangle_id = 0;



    for (const int& x: vertices){

        vector<int> x_in_neighbors = in_neighbors[x];
        sort(x_in_neighbors.begin(), x_in_neighbors.end(), [&pair_crossings, &x](int a, int b) {
            return pair_crossings[a][x] - pair_crossings[x][a] > pair_crossings[b][x] - pair_crossings[x][b];
        });

        for (const int& y: x_in_neighbors){
            if (used[x][y]) continue;
            int yx = pair_crossings[y][x] - pair_crossings[x][y];

            int bestz = -1;
            int best_w = 0;

            for (const int& z: out_neighbors[x]){
                if (used[x][y] || used[x][z] || used[y][z]) continue;
                int zy = pair_crossings[z][y] - pair_crossings[y][z];
                if (zy > 0 ){
                    int xz = pair_crossings[x][z] - pair_crossings[z][x];
                    int weight = min(min(yx, xz),zy);

                    if (weight > best_w){
                        bestz = z;
                        best_w = weight;
                    }
                }
            }

            if (bestz >= 0){
                int z = bestz;
                int weight = best_w;
                // triangle y > x -> z -> y
                used[x][y] = true;
                used[y][x] = true;
                used[x][z] = true;
                used[z][x] = true;
                used[z][y] = true;
                used[y][z] = true;
                // triangles.push_back({x,z,y});
                // printf("triangle %d: %d %d %d\n", weight, x, y, z);
                triangles_adj[x].push_back(triangle_id);
                triangles_adj[y].push_back(triangle_id);
                triangles_adj[z].push_back(triangle_id);
                triangles_weight.push_back(weight);
                triangle_id ++;

            }
            
        }
    }
    return make_pair(triangles_adj, triangles_weight);
}




pair<int, vector<vector<vector<int>>>> find_edge_disjoint_triangles_greedy2(
    const vector<bool>& mask,
    const vector<int>& vertices,
    const vector<vector<vector<int>>>& triangles){

    vector<vector<bool>> used(mask.size());
    for (int i= 0; i < used.size(); ++i){
        used[i] = vector<bool>(mask.size(), false);
    }
    int total = 0;
    vector<vector<vector<int>>> triangles_adj(mask.size());



    for (const int& x: vertices){

        // int best_triangle_index = -1;
        // int best_weight = 0;

        for (int i = 0; i < triangles[x].size(); ++i){
            auto triangle = triangles[x][i];
            int y = triangle[0];
            int z = triangle[1];
            if (mask[y] && mask[z] &&
            used[x][y] == false && used[x][z] == false && used[y][z] == false){
                int weight = triangle[2];

                 used[x][y] = true;
                used[y][x] = true;
                used[x][z] = true;
                used[z][x] = true;
                used[z][y] = true;
                used[y][z] = true;
                triangles_adj[x].push_back({z,y,weight});
                triangles_adj[y].push_back({x,z,weight});
                triangles_adj[z].push_back({x,y,weight});

                total += weight;

                // if (weight > best_weight){
                //     best_weight = weight;
                //     best_triangle_index = i;
                // }
            }
        }

        // if (best_triangle_index >= 0){
        //     int i = best_triangle_index;
        //     int weight = best_weight;
        //     auto triangle = triangles[x][i];
        //     int y = triangle[0];
        //     int z = triangle[1];
        //     used[x][y] = true;
        //     used[y][x] = true;
        //     used[x][z] = true;
        //     used[z][x] = true;
        //     used[z][y] = true;
        //     used[y][z] = true;
        //     triangles_adj[x].push_back({z,y,weight});
        //     triangles_adj[y].push_back({x,z,weight});
        //     triangles_adj[z].push_back({x,y,weight});

        //     total += weight;
        // }
    }
    return make_pair(total, triangles_adj);
}



int bad_cr_lower_bound(const vector<vector<int>>& pair_crossings, const vector<vector<int>>& in_neighbors, const vector<vector<int>>& out_neighbors, const vector<int>& vertices){
    vector<vector<bool>> used(pair_crossings.size());
    for (int i = 0; i < pair_crossings.size(); ++i){
        used[i] = vector<bool>(pair_crossings.size());
    }


    int total = 0;
    for (int i = 0; i < vertices.size(); ++i){
        int x = vertices[i];

        int best = 0;
        int best_j = -1;

        for (int j= 0; j < vertices.size(); ++j){
            if (i == j) continue;
            int y = vertices[j];

            int wxy = 0;
            if (used[y][x] == false && pair_crossings[y][x] - pair_crossings[x][y] > 0){
                wxy += pair_crossings[y][x] - pair_crossings[x][y];
            }
            for (const int& z: out_neighbors[y]){
                int zx = pair_crossings[z][x] - pair_crossings[x][z];
                if (used[z][x] == false &&  used[y][z] == false && zx > 0){
                    wxy += min(zx, pair_crossings[y][z] - pair_crossings[z][y]);
                }
            }
            int wyx = 0;
            if (used[x][y] == false && pair_crossings[x][y] - pair_crossings[y][x] > 0){
                wyx += pair_crossings[x][y] - pair_crossings[y][x];
            }
            for (const int& z: out_neighbors[x]){
                int zy = pair_crossings[z][y] - pair_crossings[y][z];
                if (used[z][x] == false && used[z][y] == false && zy > 0){
                    wyx += min(zy, pair_crossings[x][z] - pair_crossings[z][x]);
                }
            }
            if (min(wxy, wyx) > best){
                best_j = j;
                best = min(wxy, wyx);
            }
        }

        if (best_j >= 0){
            int j = best_j;
            
            int y = vertices[j];

            if (used[y][x] == false && pair_crossings[y][x] - pair_crossings[x][y] > 0){
                used[y][x] = true;
                used[x][y] = true;
            }
            for (const int& z: out_neighbors[y]){
                int zx = pair_crossings[z][x] - pair_crossings[x][z];
                if (used[z][x] == false &&  used[y][z] == false && zx > 0){
                    used[z][x] = true;
                    used[x][z] = true;
                    used[z][y] = true;
                    used[y][z] = true;
                }
            }
            if (used[x][y] == false && pair_crossings[x][y] - pair_crossings[y][x] > 0){
                used[y][x] = true;
                used[x][y] = true;
            }
            for (const int& z: out_neighbors[x]){
                int zy = pair_crossings[z][y] - pair_crossings[y][z];
                if (used[z][x] == false && used[z][y] == false && zy > 0){
                    used[z][x] = true;
                    used[x][z] = true;
                    used[z][y] = true;
                    used[y][z] = true;
                }
            }
            total += best;
        }
        

        
    }
    return total;
}







int find_triangles_from_order(const vector<vector<int>>& pair_crossings,
const vector<vector<int>>& in_neighbors, 
const vector<int>& order,
const vector<int>& vertices){
    vector<vector<bool>> used(pair_crossings.size());
    for (int i = 0; i < pair_crossings.size(); ++i){
        used[i] = vector<bool>(pair_crossings.size());
    }
    int total = 0;

    for (int i = 0; i < vertices.size(); ++i){
        int x = vertices[i];
        
        for (int j = i+1; j < vertices.size(); ++j){
            int y = vertices[j];
            int yx = pair_crossings[y][x] - pair_crossings[x][y];
            if (yx > 0 && used[x][y] == false){
                for (int k = i+1; k < j; ++k){
                    int z = vertices[k];
                    int xz = pair_crossings[x][z] - pair_crossings[z][x];
                    int zy = pair_crossings[z][y] - pair_crossings[y][z];
                    if (xz > 0 && zy > 0 && used[x][z] == false && used[z][y] == false){
                        int w = min(min(xz, zy), yx);
                        used[x][y] = true;
                        used[y][x] = true;
                        used[x][z] = true;
                        used[z][x] = true;
                        used[z][y] = true;
                        used[y][z] = true;
                        total += w;
                        break;
                    }
                }
            }
        }
    }
    return total;
}





/**
 * This one should be faster in linear time
*/
int find_triangles_from_order2(const vector<vector<int>>& pair_crossings,
    const vector<int>& order,
    const vector<int>& vertices){

    vector<vector<bool>> used(pair_crossings.size());
    for (int i = 0; i < pair_crossings.size(); ++i){
        used[i] = vector<bool>(pair_crossings.size());
    }
    int total = 0;
    int j = 2;
    int i = 0;
    while(i < vertices.size() && j< vertices.size()){
        int x = vertices[i];
        int y = vertices[j];
        int z = vertices[i+1];


        int yx = pair_crossings[y][x] - pair_crossings[x][y];
        if (yx > 0){
            int xz = pair_crossings[x][z] - pair_crossings[z][x];
            int zy = pair_crossings[z][y] - pair_crossings[y][z];
            if (xz > 0 && zy > 0 ){
                int w = min(min(xz, zy), yx);
                total += w;
                i ++;
                if (i+1 >= j){
                    j ++;
                }
            } else {
                i ++;
            }
        } else {
            j ++;
        }
    }
    return total;
}




void check_digraph_twins(const vector<int>& vertices, const vector<vector<int>>& pair_crossings){
    for (const int& x: vertices){
        for (const int& y: vertices){
            if (x < y && pair_crossings[x][y] == pair_crossings[y][x]){
                bool is_twin = true;
                for (const int& z: vertices){
                    if (x != z && y != z){
                        if (pair_crossings[x][z] != pair_crossings[y][z] || pair_crossings[z][x] != pair_crossings[z][y]){
                            is_twin = false;
                            break;
                        }
                    }
                }
                if (is_twin){
                    cout << "twins: " << x << " " << y << endl;
                }
            }
        }
    }
}