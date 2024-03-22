#include <iostream>
#include <vector>
#include <string>
#include <iomanip> 
#include <algorithm>
#include <random>

#include "solver.h"

using std::vector;

int main(int argc, char* argv[]) {


    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " <filename>\n";
        return 1;
    }

    vector<vector<int>> adj = load_file(argv[1]);
    // std::cout << "Adjacencies: " << "\n";
    // print_adj(adj);
    // std::cout << "Lower bound: " << lower_bound(adj) << "\n";

    vector<int> pos = greedy_sequential(adj);
    // std::cout << "Number of crossings: " << nb_crossings(adj, pos) << "\n";

    int lb = lower_bound(adj);
    int nc = nb_crossings(adj, pos);
    std::cout << argv[1] << " " << lb << " " << nc << " ";
    

    std::cout << "\n";

    return 0;
}
