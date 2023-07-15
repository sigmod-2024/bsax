
#ifndef BSAX_LEAFNODE_H
#define BSAX_LEAFNODE_H

#include "LeafKey.h"

class LeafNode {
public:
    LeafNode() : len(0){
        sax_lb.set_max_value();
        sax_ub.set_min_value();
    }

    SAX sax_lb;
    SAX sax_ub;

    u_int32_t len;
    LeafKey leaf_keys[LEAF_MAX_NUM];

#if exp3_b_binary_use_isax_to_prune
    SAX sax_;
    CARD card_;
#endif
};


#endif 
