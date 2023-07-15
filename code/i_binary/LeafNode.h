



#ifndef BSAX_LEAFNODE_H
#define BSAX_LEAFNODE_H

#include "LeafKey.h"

class LeafNode {
public:
    LeafNode(SAX sax, CARD card) : sax_(sax), card_(card), len(0){
        memset(sum, 0, sizeof(sum));
        memset(square_sum, 0, sizeof(square_sum));

#if exp3_i_binary_use_bsax_to_prune
        sax_lb.set_max_value();
        sax_ub.set_min_value();
#endif
    }
    SAX sax_;
    CARD card_;
    u_int32_t len;
    LeafKey leaf_keys[LEAF_MAX_NUM];

    u_int32_t sum[SEGMENTS];
    u_int32_t square_sum[SEGMENTS];

#if exp3_i_binary_use_bsax_to_prune
    SAX sax_lb;
    SAX sax_ub;
#endif
};


#endif 
