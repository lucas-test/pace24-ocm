#include <iostream>
#include <vector>
#include <string>
#include <iomanip> 
#include <algorithm>
#include <random>
#include <chrono>

#include "common.h"
#include "solver1.h"

using std::vector;

void search_random(void);


int main(int argc, char* argv[]) {
    
    // search_random();


    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " <filename>\n";
        return 1;
    }

    std::cout << argv[1] << std::endl;


    vector<vector<int>> adj = load_file(argv[1]);
    std::cout << adj.size() << "\n";
    // print_adj(adj);

    reduce_degree_0(adj);
    std::cout << adj.size() << "\n";
    print_adj(adj);


    // Solver1 test
    std::cout << "Lower bound: " << lower_bound(adj) << "\n";
    solver1(adj);


    // std::cout << "Adjacencies: " << "\n";
    // print_adj(adj);
    // std::cout << "Lower bound: " << lower_bound(adj) << "\n";

    // vector<int> pos = greedy_sequential(adj);
    // std::cout << "Number of crossings: " << nb_crossings(adj, pos) << "\n";

    // int lb = lower_bound(adj);
    // int nc = nb_crossings(adj, pos);
    // std::cout << argv[1] << " " << lb << " " << nc << " ";
    


    // sort by increasing max neighbor
    // std::sort(adj.begin(), adj.end(), [](const vector<int>& a, const vector<int>& b) {
    //    if (a.empty() || b.empty()) {
    //         return a.size() < b.size();
    //     } else {
    //         return *std::max_element(a.begin(), a.end()) < *std::max_element(b.begin(), b.end());
    //     }
    // });
    // std::vector<int> posGIMX = greedy_sequential(adj);
    // std::cout << nb_crossings(adj, posGIMX) << "\n";
    // std::vector<int> orderGIMX(posGIMX.size());
    // for (int i = 0; i < posGIMX.size(); i++) {
    //     orderGIMX[posGIMX[i]] = i ;
    // }
    // std::cout << nb_crossings_from_order(adj, orderGIMX) << "\n";
    // // print(orderGIMX);

    //  // sort by decreasing max neighbor
    // std::sort(adj.begin(), adj.end(), [](const vector<int>& a, const vector<int>& b) {
    //    if (a.empty() || b.empty()) {
    //         return a.size() < b.size();
    //     } else {
    //         return *std::max_element(a.begin(), a.end()) > *std::max_element(b.begin(), b.end());
    //     }
    // });
    // std::cout << nb_crossings(adj, greedy_sequential(adj)) << " ";

   

    // std::cout << "\n";

    return 0;
}


/**
 * @brief for finding a hard configuration with not that much vertices: it should be done in some seconds by solver1
 * 
 */
void search_random(void){
    srand(static_cast<unsigned>(time(0)));

    for (int i = 0; i < 10 ; i ++){
        vector<vector<int>> adj = generate_random_adj(20, 20, 0.33);
        vector<int> pos = greedy_sequential(adj);
        if (lower_bound(adj) != nb_crossings(adj, pos)){

            auto start = std::chrono::high_resolution_clock::now();
            solver1(adj);
            auto end = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> elapsed = end - start;

            

            if (elapsed.count() > 0.5){
                print_adj(adj);
                print_gr_format(adj);
                std::cout << "Lower bound: " << lower_bound(adj) << "\n";
                std::cout << "Execution time of solver1: " << elapsed.count() << " seconds\n";
            }

           

            
        }
    }
}