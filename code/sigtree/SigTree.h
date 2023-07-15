



#ifndef BSAX_ISAX2TREE_H
#define BSAX_ISAX2TREE_H

#include <vector>
#include "RootNode.h"

class SigTree {
private:
    void DFS(const SAXT *search_saxt, void *node, bool is_leaf, std::vector<LeafKey> &leaf_ans, u_int8_t card_now) const;
    void BFS(ts_type* search_paa, float top_dis, std::vector<LeafKey>& leaf_ans) const;
    void BFSM(ts_type* search_paa, float top_dis, std::vector<LeafNode*>& leaf_ans) const;
public:
    SigTree(const char* _filename): filename(_filename) {}
    void Build();
    void BuildTree();
    void BuildFromSAX(const char *sax_filename);
    void Insert(const SAXT* insert_saxt, u_int64_t p);
    void ApproximateSearch(int k, TS *search_ts, std::vector<TS *> &ts_ans, std::vector<float> &top_dis) const;
    void ExactSearch(int k, TS *search_ts, std::vector<TS *> &ts_ans, std::vector<float> &top_dis) const;
    void ExactSearchM(int k, TS *search_ts, std::vector<TS *> &ts_ans, std::vector<float> &top_dis) const;
    void ExactSearchExp3(TS* search_ts, float top_k_dis, std::vector<std::pair<SAX, u_int64_t>>& sax_and_p, std::vector<float>& top_dis) const;
    void materialized(const char* _filename);
    void materialized();

    RootNode* root;
    const char* filename;

};


#endif 
