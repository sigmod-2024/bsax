



#ifndef BSAX_LEAFNODE_H
#define BSAX_LEAFNODE_H

#include "LeafKey.h"
#include <vector>

class LeafNode {
public:
    LeafNode(SAXT saxt, u_int8_t card) : saxt_(saxt), card_(card){}
    SAXT saxt_;
    u_int8_t card_;
#if sigtree_leaf_keys_use_vector
    std::vector<LeafKey> leaf_keys;
#else
    u_int32_t len = 0;
    LeafKey leaf_keys[LEAF_MAX_NUM];
#endif
};


#endif 
