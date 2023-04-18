#pragma once 

#include <execution>

using RI_Kmer = container::RI_Kmer; 

namespace b_tree {

class B_Tree {
    const int B = 3;
    int n; 
    int nblocks; 
    int t = 0; 
    int max_level; 

    std::vector<std::array<RI_Kmer,3>> btree{};

    public: 
        B_Tree() = default;
        ~B_Tree() = default;

        B_Tree(std::vector<RI_Kmer> container) {
            n = container.size();
            nblocks = (n + B - 1) / B;
            btree.resize(nblocks);
            build(container);
            max_level = (log(nblocks) + 1) / log(B + 1); 
        } 

        int get_nblocks() { return nblocks; }

        int get_B() { return B; }

        int go(int k, int i) { return k * (B + 1) + i + 1; }

        std::tuple<int,int> go_back(int k) { 
            int old_i = (k - 1) % (B + 1); 
            int old_k = (k - old_i - 1) / (B + 1);
            return std::make_tuple(old_k, old_i);  
        }

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

        std::tuple<int, std::tuple<int,int>> search(const int i_rep) {
            int k = 0; 
            int level = max_level; 
            auto res = std::make_tuple(-1, std::make_tuple(0, 0));
            start: 
            while (k < nblocks) {
                for (int i = 0; i < B; i++) {
                    if (btree[k][i].integer_rep == i_rep) {
                        return std::make_tuple(level, std::make_tuple(k, i)); 
                    } else if (btree[k][i].integer_rep > i_rep) {
                        res = std::make_tuple(level, std::make_tuple(k, i));
                        k = go(k, i);
                        level--; 
                        goto start;
                    }
                }
                k = go(k, B);
                level--; 
            }
            return res;
        }

        std::tuple<bool,RI_Kmer> search(const int i_rep, int search_from_k) {
            int k = search_from_k; 
            auto res = std::make_tuple(false, RI_Kmer{});
            start: 
            while (k < nblocks) {
                for (int i = 0; i < B; i++) {
                    if (btree[k][i].integer_rep >= i_rep) {
                        res = std::make_tuple(true, btree[k][i]);
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