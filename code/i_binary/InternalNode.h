



#ifndef BSAX_INTERNALNODE_H
#define BSAX_INTERNALNODE_H

#include "LeafNode.h"
#include "sax.h"

class InternalNode {
public:
    InternalNode(SAX _sax, CARD _card, LeafNode* _left, LeafNode* _right, u_int8_t _sp) : sax_(_sax), card_(_card),
    left(_left), right(_right), split_segment(_sp), is_left_leaf(true), is_right_leaf(true) {}

    void* left;
    void* right;
    SAX sax_;
    CARD card_;
    u_int8_t split_segment; 

    bool is_left_leaf;   
    bool is_right_leaf;   

#if exp3_i_binary_use_bsax_to_prune
    SAX sax_lb;
    SAX sax_ub;
#endif
};


#endif 
