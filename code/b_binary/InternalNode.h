#ifndef BSAX_INTERNALNODE_H
#define BSAX_INTERNALNODE_H

#include "LeafNode.h"
#include "sax.h"

class InternalNode {
public:
    InternalNode(SAX _sax_lb, SAX _sax_ub, LeafNode* _left, LeafNode* _right, u_int8_t _sp, sax_type _split_sax) :
            sax_lb(_sax_lb), sax_ub(_sax_ub), left(_left), right(_right),
            split_segment(_sp), split_sax(_split_sax),
            is_left_leaf(true), is_right_leaf(true) {}

    void* left;
    void* right;
    SAX sax_lb;
    SAX sax_ub;
    u_int8_t split_segment; 
    sax_type split_sax;  

    bool is_left_leaf;   
    bool is_right_leaf;   

#if exp3_b_binary_use_isax_to_prune
    SAX sax_;
    CARD card_;
#endif
};


#endif 
