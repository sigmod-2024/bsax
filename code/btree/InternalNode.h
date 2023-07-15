



#ifndef BSAX_INTERNALNODE_H
#define BSAX_INTERNALNODE_H


#include "InternalKey.h"

class InternalNode {
public:
    InternalNode() : len(0) {}

    u_int64_t len;
    InternalKey internal_keys[LEAF_MAX_NUM];
};


#endif 
