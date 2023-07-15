



#ifndef BSAX_ROOTNODE_H
#define BSAX_ROOTNODE_H

#include "InternalNode.h"
#include <unordered_map>
#include <map>
class RootNode {
public:
    

#if sigtree_delay_create_leafnode

    std::map<u_int32_t, std::pair<void*, bool>> children;
#else
    std::pair<void*, bool> children[(u_int64_t)1 << SEGMENTS];
#endif
};


#endif 
