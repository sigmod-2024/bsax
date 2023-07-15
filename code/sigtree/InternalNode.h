//
// Created by seth on 5/22/23.
//

#ifndef BSAX_INTERNALNODE_H
#define BSAX_INTERNALNODE_H

#include "LeafNode.h"
#include "sax.h"
#include <unordered_map>
#include <map>
class InternalNode {
public:
    InternalNode(SAXT _saxt, u_int8_t _card) : saxt_(_saxt), card_(_card) {}


    SAXT saxt_;
    u_int8_t card_;
    // key: Segmentsax   value: first:is_leaf, second:
#if sigtree_delay_create_leafnode
//    std::unordered_map<u_int32_t, std::pair<void*, bool>> children;
    std::map<u_int32_t, std::pair<void*, bool>> children;
#else
    std::pair<void*, bool> children[(u_int64_t)1 << SEGMENTS];
#endif

};


#endif //BSAX_INTERNALNODE_H
