



#ifndef BSAX_LEAFNODE_H
#define BSAX_LEAFNODE_H


#include "LeafKey.h"

class LeafNode {
public:
    LeafNode(): len(0) {}

    void add_batch(LeafKey* add_leaf_keys, u_int32_t num);

    u_int32_t len;
    LeafKey leaf_keys[LEAF_MAX_NUM];
    u_int32_t id;   

};


#endif 
