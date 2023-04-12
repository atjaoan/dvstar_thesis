#pragma once 

#include <execution>

using RI_Kmer = container::RI_Kmer; 

namespace b_tree {

class B_Tree {
    const int B = 16;
    int n; 
    int nblocks; 
    int t = 0; 

    std::vector<std::array<RI_Kmer,16>> btree{};

    public: 

        B_Tree() = default;
        ~B_Tree() = default;

        B_Tree(std::vector<RI_Kmer> container) {
            n = container.size();
            nblocks = (n + B - 1) / B;
            btree.resize(nblocks);
            build(container);
        } 

        int get_nblocks() { return nblocks; }

        int get_B() { return B; }

        std::array<RI_Kmer, 16> &get_block(int k) { return btree[k]; }

        int go(int k, int i) { return k * (B + 1) + i + 1; }

        int go_back(int k, int i) { return (k - i - 1) / (B + 1)}

        void build(std::vector<RI_Kmer> &container, int k = 0) {
            if (k < nblocks) {
                for (int i = 0; i < B; i++) {
                    build(container, go(k, i));
                    btree[k][i] = (t < n ? container[t++] : RI_Kmer{} );
                }
                build(container, go(k, B));
            }
        }

        RI_Kmer &get(int k, int i){
            return btree[k][i]; 
        }

        RI_Kmer &get(std::tuple<int, int> idx){
            return btree[std::get<0>(idx)][std::get<1>(idx)]; 
        }

        std::tuple<int,int> search(const int i_rep) {
            int k = 0; 
            auto res = std::make_tuple(0, 0);
            start: 
            while (k < nblocks) {
                for (int i = 0; i < B; i++) {
                    if (btree[k][i].integer_rep >= i_rep) {
                        res = std::make_tuple(k, i);
                        k = go(k, i);
                        goto start;
                    }
                }
                k = go(k, B);
            }
            return res;
        }

        void prettyPrint(int max = 10){
            int k = 0;
            while (k < nblocks) {
                if (k < max) {
                    for (int i = 0; i < B; i++){
                        std::cout << btree[k][i].integer_rep << " "; 
                    }
                    std::cout << std::endl; 
                    k++; 
                } else {
                    break; 
                }
            }
        }
};
}