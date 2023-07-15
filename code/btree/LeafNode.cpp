



#include "LeafNode.h"

void LeafNode::add_batch(LeafKey *add_leaf_keys, u_int32_t num) {
    memcpy(leaf_keys + len, add_leaf_keys, num * sizeof(LeafKey));
    len += num;
}
