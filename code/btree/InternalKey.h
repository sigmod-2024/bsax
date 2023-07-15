



#ifndef BSAX_INTERNALKEY_H
#define BSAX_INTERNALKEY_H


#include "globals.h"

class InternalKey {
public:
    InternalKey(): is_leaf(false), p(nullptr) {}

    SAXT saxt_lb;

#if btree_use_bsax
    SAX sax_lb;
    SAX sax_ub;
#else
    SAX sax_;   
    CARD card_; 
#endif

    void* p;
    bool is_leaf;
};


#endif 
