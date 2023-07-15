#ifndef BSAX_ISAX2TREE_H

#include <vector>
#include "RootNode.h"

class bSAX2Tree {
private:
    void DFS(const SAX *search_sax, void *node, bool is_leaf, std::vector<LeafKey> &leaf_ans) const;
    void BFS(ts_type* search_paa, float top_dis, std::vector<LeafKey>& leaf_ans) const;
    void BFSM(ts_type* search_paa, float top_dis, std::vector<LeafNode*>& leaf_ans) const;
    void read_ts_batch(int n, FILE* data_file, std::vector<SAX>& sax_vec, std::vector<TS>& ts_vec);
public:
    bSAX2Tree(const char* _filename): filename(_filename) {}
    void Build();
    void BuildFromSAX(const char* sax_filename);
    void BuildTree();
    void Insert(const SAX* insert_sax, const u_int64_t p);
    void ApproximateSearch(int k, TS *search_ts, std::vector<TS *> &ts_ans, std::vector<float> &top_dis) const;
    void ExactSearch(int k, TS *search_ts, std::vector<TS *> &ts_ans, std::vector<float> &top_dis) const;
    void ExactSearchM(int k, TS *search_ts, std::vector<TS *> &ts_ans, std::vector<float> &top_dis) const;
    void ExactSearchExp3(TS* search_ts, float bsf, std::vector<std::pair<SAX, u_int64_t>>& sax_and_p, std::vector<float>& top_dis) const;
    void materialized(const char* _filename);
    void materialized();

    RootNode* root;
    const char* filename;

};


#endif 
