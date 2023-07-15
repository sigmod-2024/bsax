#ifndef BSAX_ROOTNODE_H
#define BSAX_ROOTNODE_H

#include "InternalNode.h"

class RootNode {
public:
#if binary_tree_root_full
    
    std::pair<void*, bool> children[(u_int64_t)1 << SEGMENTS];
# else
    void* node;
    bool root_is_leaf = true;
#endif
};


#endif 
