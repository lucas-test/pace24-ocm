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
void compare_greedy_insertion_orders(const std::vector<std::vector<int>>& adj);


int main(int argc, char* argv[]) {
    
    // search_random();


    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " <filename>\n";
        return 1;
    }

    std::cout << "\n" << argv[1] << std::endl;


    vector<vector<int>> adj = load_file(argv[1]);

    reduce_degree_0(adj);
    // std::cout << adj.size() << "\n";


    // Solver1 test
    std::cout << "lower bound: " << lower_bound(adj) << "\n";

    // std::cout << "lol: " << nb_crossings(adj, greedy_sequential(adj)) << "\n";

    // std::sort(adj.begin(), adj.end(), [](const std::vector<int>& a, const std::vector<int>& b) {
    //     return a.size() < b.size();
    // });

    //     std::cout << "lol: " << nb_crossings(adj, greedy_sequential(adj)) << "\n";


    // std::sort(adj.begin(), adj.end(), [](const std::vector<int>& a, const std::vector<int>& b) {
    //     return a.size() > b.size();
    // });

    //     std::cout << "lol: " << nb_crossings(adj, greedy_sequential(adj)) << "\n";


    // std::sort(adj.begin(), adj.end(), [](const vector<int>& a, const vector<int>& b) {
    //    if (a.empty() || b.empty()) {
    //         return a.size() < b.size();
    //     } else {
    //         return *std::max_element(a.begin(), a.end()) > *std::max_element(b.begin(), b.end());
    //     }
    // });
    // print_adj(adj);

    solver1(adj);

        // std::cout << "lol: " << nb_crossings(adj, greedy_sequential(adj)) << "\n";



    // // Sort by increasing degree
    // std::sort(adj.begin(), adj.end(), [](const std::vector<int>& a, const std::vector<int>& b) {
    //     return a.size() < b.size();
    // });
    // int greedy_inc_deg = nb_crossings(adj, greedy_sequential(adj));
    // std::cout << "greedy increasing degree: " << greedy_inc_deg << "\n";

    // // Sort by decreasing degree
    // std::sort(adj.begin(), adj.end(), [](const std::vector<int>& a, const std::vector<int>& b) {
    //     return a.size() > b.size();
    // });
    // vector<int> greedy_dec_deg_pos = greedy_sequential(adj);
    // int greedy_dec_deg = nb_crossings(adj, greedy_sequential(adj));
    // std::cout << "greedy decreasing degree: " << greedy_dec_deg << "\n";
    
    // vector<int> greedy_dec_deg_order(greedy_dec_deg_pos.size());
    // for (int i = 0; i < greedy_dec_deg_pos.size(); i++) {
    //     greedy_dec_deg_order[greedy_dec_deg_pos[i]] = i ;
    // }
    // int greedy_dec_deg_order_nc = nb_crossings_from_order(adj, greedy_dec_deg_order);
    // std::cout << "greedy decreasing degree order: " << greedy_dec_deg_order_nc << "\n";


    return 0;
}


void compare_greedy_insertion_orders(const vector<vector<int>>& adj){

    // // Natural
    // int greedy_nat = nb_crossings(adj, greedy_sequential(adj));
    // // std::cout << "greedy: " << greedy_nat << "\n";

    // // Sort by increasing degree
    // std::sort(adj.begin(), adj.end(), [](const std::vector<int>& a, const std::vector<int>& b) {
    //     return a.size() < b.size();
    // });
    // int greedy_inc_deg = nb_crossings(adj, greedy_sequential(adj));
    // // std::cout << "greedy increasing degree: " << greedy_inc_deg << "\n";

    // // Sort by decreasing degree
    // std::sort(adj.begin(), adj.end(), [](const std::vector<int>& a, const std::vector<int>& b) {
    //     return a.size() > b.size();
    // });
    // int greedy_dec_deg = nb_crossings(adj, greedy_sequential(adj));
    // // std::cout << "greedy decreasing degree: " << greedy_dec_deg << "\n";

    // // Sort by increasing max neighbor
    // std::sort(adj.begin(), adj.end(), [](const vector<int>& a, const vector<int>& b) {
    //    if (a.empty() || b.empty()) {
    //         return a.size() < b.size();
    //     } else {
    //         return *std::max_element(a.begin(), a.end()) < *std::max_element(b.begin(), b.end());
    //     }
    // });
    // std::vector<int> posGIMX = greedy_sequential(adj);
    // int greedy_inc_max_neighbor = nb_crossings(adj, greedy_sequential(adj));
    // // std::cout << "greedy increasing max neighbor: " << greedy_inc_max_neighbor  << "\n";
    
    // // std::vector<int> orderGIMX(posGIMX.size());
    // // for (int i = 0; i < posGIMX.size(); i++) {
    // //     orderGIMX[posGIMX[i]] = i ;
    // // }
    // // std::cout << nb_crossings_from_order(adj, orderGIMX) << "\n";
    // // print(orderGIMX);

    //  // Sort by decreasing max neighbor
    // std::sort(adj.begin(), adj.end(), [](const vector<int>& a, const vector<int>& b) {
    //    if (a.empty() || b.empty()) {
    //         return a.size() < b.size();
    //     } else {
    //         return *std::max_element(a.begin(), a.end()) > *std::max_element(b.begin(), b.end());
    //     }
    // });
    // int greedy_dec_max_neighbor = nb_crossings(adj, greedy_sequential(adj));
    // // std::cout <<  "greedy decreasing max neighbor: " << greedy_dec_max_neighbor << "\n";
   
    // std::cout << lower_bound(adj) << " " << greedy_nat << " " << greedy_inc_deg << " " << greedy_dec_deg << " " << greedy_inc_max_neighbor << " " << greedy_dec_max_neighbor << "\n";

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