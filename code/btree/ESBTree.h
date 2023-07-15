//
// Created by seth on 5/19/23.
//

#ifndef BSAX_ESBTREE_H
#define BSAX_ESBTREE_H


#include "InternalNode.h"
#include "LeafKey.h"
#include "LeafNode.h"
#include <vector>
using namespace std;

class ESBTree {
private:
    void read_ts_batch(int n, FILE* data_file, vector<LeafKey>& leaf_keys, long offset, vector<TS>& ts_vec, vector<SAX>& sax_vec);
    void DFS(const SAXT* search_saxt, void* node, bool is_leaf, vector<LeafNode*>& leaf_ans) const;
    void BFS(ts_type* search_paa, float top_dis, vector<LeafNode*>& leaf_ans) const;
public:
    ESBTree(const char* _filename): filename(_filename) {}
    void Build();
    void BuildFromSAX(const char* sax_filename);
    void BuildTree(vector<LeafKey>& leaf_keys);
    void ApproximateSearch(int k, TS* search_ts, vector<TS*>& ts_ans, vector<float>& top_dis) const;
    void ExactSearch(int k, TS* search_ts, vector<TS*>& ts_ans, vector<float>& top_dis) const;
    void ExactSearchM(int k, TS* search_ts, vector<TS*>& ts_ans, vector<float>& top_dis) const;
    void ExactSearchExp3(TS *search_ts, float top_k_dis, vector<pair<SAX, u_int64_t>> &sax_and_p, vector<float> &top_dis) const;
    void materialized(const char* _filename);
    void materialized();

    InternalNode* root;
    LeafNode* leaf_nodes;  // 
    u_int32_t num_leafs;    // 
    const char* filename;


};


#endif //BSAX_ESBTREE_H
