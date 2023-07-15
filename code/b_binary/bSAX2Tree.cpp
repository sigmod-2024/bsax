#include <vector>
#include <algorithm>
#include "bSAX2Tree.h"
#include "RootNode.h"
#include "../util/CntRecord.h"
#include "../util/TopKHeap.h"
#include "TimeRecord.h"
#include "unordered_set"

using namespace std;

void bSAX2Tree::read_ts_batch(int n, FILE* data_file, vector<SAX>& sax_vec, vector<TS>& ts_vec) {
    BUILD_READ_TS_START
    for(int i = 0; i < n; i ++) {
        fread(&ts_vec[i], sizeof(TS) , 1, data_file);
    }
    BUILD_READ_TS_END

    BUILD_CONVERT_SAX_START
    for (int i = 0; i < n; i ++ ) {
        sax_from_ts(ts_vec[i].ts, sax_vec[i].sax);
    }
    BUILD_CONVERT_SAX_END
}

void bSAX2Tree::Build() {
    this->BuildTree();

    FILE * data_file = fopen(filename,"r");
    if (!data_file) {
        cout << "" << filename << endl;
        exit(-1);
    }
    vector<SAX> sax_vec(READ_TS_BATCH);
    vector<TS> ts_vec(READ_TS_BATCH);

    BUILD_TIME_START
    for (int i = 0; i < (TOTAL_TS - 1) / READ_TS_BATCH + 1; i ++ ) {
        int read_batch = READ_TS_BATCH;
        if (TOTAL_TS % READ_TS_BATCH && i == (TOTAL_TS - 1) / READ_TS_BATCH) {
            read_batch = TOTAL_TS % READ_TS_BATCH;
        }
        read_ts_batch(read_batch, data_file, sax_vec, ts_vec);
        BUILD_INDEX_START
        for (int j = 0; j < read_batch; j ++ ) {
            this->Insert(&sax_vec[j], i * READ_TS_BATCH + j);
        }
        BUILD_INDEX_END
    }

    BUILD_TIME_END

    fclose (data_file);
}

void bSAX2Tree::BuildFromSAX(const char *sax_filename) {
    this->BuildTree();

    FILE * data_file = fopen(sax_filename,"r");
    if (!data_file) {
        cout << "" << sax_filename << endl;
        exit(-1);
    }

    vector<SAX> sax_vec(TOTAL_TS);
    for(int i = 0; i < TOTAL_TS; i++) {
        fread(&sax_vec[i], sizeof(SAX), 1, data_file);
    }

    BUILD_INDEX_START
    for(int i = 0; i < TOTAL_TS; i++) {
        this->Insert(&sax_vec[i], i);
//        if ((i + 1) % 1000000 == 0) {
//            cout << "" << i + 1 << endl;
//        }
    }
    BUILD_INDEX_END
    fclose (data_file);
}


/**
 * 
 */
void bSAX2Tree::BuildTree() {
    root = new RootNode();
#if binary_tree_root_full
    for (u_int64_t i = 0; i < ((u_int64_t)1 << SEGMENTS); i ++ ) {
        root->children[i] = {new LeafNode(), true};

#if exp3_b_binary_use_isax_to_prune
        SAX sax_;
        CARD card_;
        for (int j = 0; j < SEGMENTS; j ++ ) {
            sax_.sax[j] = ((i >> j) & 1) << (BIT_CARDINALITY - 1);
            card_.card[j] = 1;
        }
        LeafNode* leaf_node = (LeafNode*)root->children[i].first;
        leaf_node->sax_ = sax_;
        leaf_node->card_ = card_;
#endif
    }
#else
    LeafNode* leaf_node = new LeafNode();
    root->node = leaf_node;
#if exp3_b_binary_use_isax_to_prune
    SAX sax;
    CARD card;
    sax.set_min_value();
    card.set_min_card();
    leaf_node->sax_ = sax;
    leaf_node->card_ = card;
#endif
#endif
}


//struct CompLeafKey {
//    u_int8_t s;
//    CompLeafKey(u_int8_t _s) : s(_s) {}
//    bool operator()(const LeafKey& x, const LeafKey& y) const {
//        return x.sax_.sax[s] < y.sax_.sax[s];
//    }
//};
//
//
//LeafNode* SplitNode(LeafNode* leaf_node, void*& pre_node, bool& pre_is_root, const SAX* insert_sax) {
//    int s = -1;
//    float max_dis = -MAXFLOAT;
//
//    for (int i = 0; i < Segments; i ++ ) {
//        float breakpoint_lower = sax_a[leaf_node->sax_lb.sax[i]];   // saxbreakpoint
//        float breakpoint_upper = sax_a[leaf_node->sax_ub.sax[i] + 1];   // saxbreakpoint
//
//        if (max_dis < breakpoint_upper - breakpoint_lower && leaf_node->sax_ub.sax[i] - leaf_node->sax_lb.sax[i] > 0) {
//            max_dis = breakpoint_upper - breakpoint_lower;
//            s = i;
//        }
//    }
//
//    // sax
//    sort(leaf_node->leaf_keys, leaf_node->leaf_keys + leaf_node->len, CompLeafKey(s));
//
//
//    sax_type split_sax = leaf_node->leaf_keys[(leaf_node->len >> 1) - 1].sax_.sax[s];
////    for (int i = 0; i < leaf_node->len; i ++ ) {
////        cout << "(" << i << ":" << (int)leaf_node->leaf_keys[i].sax_.sax[s] << ") ";
////    }
////    cout <<  ":" << (int) split_sax << endl;
//    int l = 0, r = leaf_node->len - 1;
//    while(l < r) {  // 
//        int mid = (l + r + 1) >> 1;
//        if (leaf_node->leaf_keys[mid].sax_.sax[s] <= split_sax) l = mid;
//        else r = mid - 1;
//    }
//    if (l == leaf_node->len - 1) {  // 
//        split_sax -- ; // -1，
//        l = 0, r = leaf_node->len - 1;
//        while(l < r) {
//            int mid = (l + r + 1) >> 1;
//            if (leaf_node->leaf_keys[mid].sax_.sax[s] <= split_sax) l = mid;
//            else r = mid - 1;
//        }
//    }
//
//
//    // 
//    LeafNode* left_leaf_node = new LeafNode();
//    LeafNode* right_leaf_node = new LeafNode();
//
//    memcpy(left_leaf_node->leaf_keys, leaf_node->leaf_keys, sizeof(LeafKey) * (l + 1));
//    left_leaf_node->len = l + 1;
//    for (int i = 0; i < left_leaf_node->len; i ++ ) {
//        for (int j = 0; j < Segments; j ++ ) {
//            left_leaf_node->sax_lb.sax[j] = min(left_leaf_node->sax_lb.sax[j], left_leaf_node->leaf_keys[i].sax_.sax[j]);
//            left_leaf_node->sax_ub.sax[j] = max(left_leaf_node->sax_ub.sax[j], left_leaf_node->leaf_keys[i].sax_.sax[j]);
//        }
//    }
//    memcpy(right_leaf_node->leaf_keys, leaf_node->leaf_keys + (l + 1), sizeof(LeafKey) * (leaf_node->len - l - 1));
//    right_leaf_node->len = leaf_node->len - l - 1;
//    for (int i = 0; i < right_leaf_node->len; i ++ ) {
//        for (int j = 0; j < Segments; j ++ ) {
//            right_leaf_node->sax_lb.sax[j] = min(right_leaf_node->sax_lb.sax[j], right_leaf_node->leaf_keys[i].sax_.sax[j]);
//            right_leaf_node->sax_ub.sax[j] = max(right_leaf_node->sax_ub.sax[j], right_leaf_node->leaf_keys[i].sax_.sax[j]);
//        }
//    }
//
//
////    cout << "s:" << s << " sax:" << (int)leaf_node->sax_lb.sax[s]  << " sax:" <<  (int)leaf_node->sax_ub.sax[s] << endl;
////    cout << "breakpoint:" << max_dis << " :" << (int) split_sax << endl;
////    std::cout << ":" << left_leaf_node->len << " :" << right_leaf_node->len << std::endl;
////    std::cout << ":" << (int)left_leaf_node->sax_lb.sax[s] << " " << (int)left_leaf_node->sax_ub.sax[s] <<
////    " :" << (int)right_leaf_node->sax_lb.sax[s] << " " << (int)right_leaf_node->sax_ub.sax[s] << std::endl;
////    if (left_leaf_node->len == Leaf_maxnum) {
////        for (int i = 0; i < leaf_node->len; i ++ ) {
////            cout << (int)leaf_node->leaf_keys[i].sax_.sax[s] << " ";
////        }
////        exit(0);
////    }
////    if (right_leaf_node->len == Leaf_maxnum) {
////        for (int i = 0; i < right_leaf_node->len; i ++ ) {
////            cout << (int)right_leaf_node->leaf_keys[i].sax_.sax[s] << " ";
////        }
////        exit(0);
////    }
//
//
//    InternalNode* new_internal_node = new InternalNode(leaf_node->sax_lb, leaf_node->sax_ub, left_leaf_node, right_leaf_node, s, split_sax);
//    if (pre_is_root) {  // root，rootmapleaf_nodeinternal_node
//        RootNode* root_node = (RootNode*) pre_node;
//        u_int64_t sax_first_bit = 0;
//        saxt_type first_bit = 1 << (Bit_cardinality - 1);
//        for (int i = 0; i < Segments; i ++ ) {
//            sax_first_bit |= ((insert_sax->sax[i] & first_bit) >> (Bit_cardinality - 1)) << i;
//        }
//        root_node->mp[sax_first_bit].first = new_internal_node;
//        root_node->mp[sax_first_bit].second = false;
//    }
//    else {
//        InternalNode* internal_node = (InternalNode*) pre_node;
//        if (internal_node->left == leaf_node) { // leaf_nodepre_node
//            internal_node->is_left_leaf = false;
//            internal_node->left = new_internal_node;
//        }
//        else {
//            internal_node->is_right_leaf = false;
//            internal_node->right = new_internal_node;
//        }
//    }
//    delete leaf_node;
//
//    pre_node = new_internal_node;
//    pre_is_root = false;
//    if (insert_sax->sax[s] <= split_sax) { // 
//        return left_leaf_node;
//    }
//    else {  // 
//        return right_leaf_node;
//    }
//}


sax_type quick_select_median(int l, int r, int k, LeafKey* leaf_keys, u_int8_t s) {
    if (l >= r) {
        return leaf_keys[l].sax_.sax[s];
    }
    int i = l - 1;
    int j = r + 1;
    sax_type mid_sax = leaf_keys[(l + r) >> 1].sax_.sax[s];
    while(i < j) {
        while(leaf_keys[ ++ i].sax_.sax[s] < mid_sax);
        while(leaf_keys[ -- j].sax_.sax[s] > mid_sax);
        if (i < j) {
            swap(leaf_keys[i], leaf_keys[j]);
        }
    }
    int num = j - l + 1;
    if (num >= k) {
        return quick_select_median(l, j, k, leaf_keys, s);
    }
    else {
        return quick_select_median(j + 1, r, k - num, leaf_keys, s);
    }
}

LeafNode* SplitNode(LeafNode* leaf_node, void*& pre_node, bool& pre_is_root, const SAX* insert_sax) {
    int s = -1;
#if b_binary_use_breakpoint_to_split
    float max_dis = -MAXFLOAT;
    for (int i = 0; i < SEGMENTS; i ++ ) {

        float breakpoint_lower = sax_a[leaf_node->sax_lb.sax[i]];   // saxbreakpoint
        float breakpoint_upper = sax_a[leaf_node->sax_ub.sax[i] + 1];   // saxbreakpoint

        if (max_dis < breakpoint_upper - breakpoint_lower && leaf_node->sax_ub.sax[i] - leaf_node->sax_lb.sax[i] > 0) {
            max_dis = breakpoint_upper - breakpoint_lower;
            s = i;
        }
    }
#else
    int max_dis = INT32_MIN;
    for (int i = 0; i < SEGMENTS; i ++ ) {
        if (max_dis < leaf_node->sax_ub.sax[i] - leaf_node->sax_lb.sax[i]) {
            max_dis = leaf_node->sax_ub.sax[i] - leaf_node->sax_lb.sax[i];
            s = i;
        }
    }
#endif
    // sax
    sax_type split_sax = quick_select_median(0, leaf_node->len - 1, leaf_node->len / 2, leaf_node->leaf_keys, s);

    // 
    LeafNode* left_leaf_node = new LeafNode();
    LeafNode* right_leaf_node = new LeafNode();


    for (int i = 0; i < leaf_node->len; i ++ ) {
        LeafKey* leaf_key = &leaf_node->leaf_keys[i];
        if (leaf_key->sax_.sax[s] <= split_sax) {   // leaf_key
            left_leaf_node->leaf_keys[left_leaf_node->len].sax_ = leaf_key->sax_;
            left_leaf_node->leaf_keys[left_leaf_node->len].p = leaf_key->p;
            left_leaf_node->len ++ ;
            for (int j = 0; j < SEGMENTS; j ++ ) {
                left_leaf_node->sax_lb.sax[j] = min(left_leaf_node->sax_lb.sax[j], leaf_key->sax_.sax[j]);
                left_leaf_node->sax_ub.sax[j] = max(left_leaf_node->sax_ub.sax[j], leaf_key->sax_.sax[j]);
            }
        }
        else {  // leaf_key
            right_leaf_node->leaf_keys[right_leaf_node->len].sax_ = leaf_key->sax_;
            right_leaf_node->leaf_keys[right_leaf_node->len].p = leaf_key->p;
            right_leaf_node->len ++ ;
            for (int j = 0; j < SEGMENTS; j ++ ) {
                right_leaf_node->sax_lb.sax[j] = min(right_leaf_node->sax_lb.sax[j], leaf_key->sax_.sax[j]);
                right_leaf_node->sax_ub.sax[j] = max(right_leaf_node->sax_ub.sax[j], leaf_key->sax_.sax[j]);
            }
        }
    }



    if (left_leaf_node->len >= leaf_node->len) {    // ,
        if (split_sax == 0) {
            cout << "0，，" << endl;
            exit(233);
        }
        split_sax --;   // -1,
        left_leaf_node->len = 0;    right_leaf_node->len = 0;
        left_leaf_node->sax_lb.set_max_value(); left_leaf_node->sax_ub.set_min_value();
        right_leaf_node->sax_lb.set_max_value(); right_leaf_node->sax_ub.set_min_value();
        for (int i = 0; i < leaf_node->len; i ++ ) {
            LeafKey* leaf_key = &leaf_node->leaf_keys[i];
            if (leaf_key->sax_.sax[s] <= split_sax) {   // leaf_key
                left_leaf_node->leaf_keys[left_leaf_node->len].sax_ = leaf_key->sax_;
                left_leaf_node->leaf_keys[left_leaf_node->len].p = leaf_key->p;
                left_leaf_node->len ++ ;
                for (int j = 0; j < SEGMENTS; j ++ ) {
                    left_leaf_node->sax_lb.sax[j] = min(left_leaf_node->sax_lb.sax[j], leaf_key->sax_.sax[j]);
                    left_leaf_node->sax_ub.sax[j] = max(left_leaf_node->sax_ub.sax[j], leaf_key->sax_.sax[j]);
                }
            }
            else {  // leaf_key
                right_leaf_node->leaf_keys[right_leaf_node->len].sax_ = leaf_key->sax_;
                right_leaf_node->leaf_keys[right_leaf_node->len].p = leaf_key->p;
                right_leaf_node->len ++ ;
                for (int j = 0; j < SEGMENTS; j ++ ) {
                    right_leaf_node->sax_lb.sax[j] = min(right_leaf_node->sax_lb.sax[j], leaf_key->sax_.sax[j]);
                    right_leaf_node->sax_ub.sax[j] = max(right_leaf_node->sax_ub.sax[j], leaf_key->sax_.sax[j]);
                }
            }
        }
    }
//    cout << "s:" << s << " sax:" << (int)leaf_node->sax_lb.sax[s]  << " sax:" <<  (int)leaf_node->sax_ub.sax[s] << endl;
//    cout << "breakpoint:" << max_dis << " :" << (int) split_sax << endl;

//    std::cout << ":" << left_leaf_node->len << " :" << right_leaf_node->len << std::endl;

//    if (left_leaf_node->len == Leaf_maxnum) {
//        for (int i = 0; i < leaf_node->len; i ++ ) {
//            cout << (int)leaf_node->leaf_keys[i].sax_.sax[s] << " ";
//        }
//        exit(0);
//    }
//    if (right_leaf_node->len == Leaf_maxnum) {
//        for (int i = 0; i < right_leaf_node->len; i ++ ) {
//            cout << (int)right_leaf_node->leaf_keys[i].sax_.sax[s] << " ";
//        }
//        exit(0);
//    }

//    cout<< (int)split_sax << endl;
//
//    cout<< (int)right_leaf_node->sax_lb.sax[s] << " " << (int)right_leaf_node->sax_ub.sax[s] << endl;
//
//    assert(left_leaf_node->len < 100);
//    assert(right_leaf_node->len < 100);

    InternalNode* new_internal_node = new InternalNode(leaf_node->sax_lb, leaf_node->sax_ub, left_leaf_node, right_leaf_node, s, split_sax);
#if exp3_b_binary_use_isax_to_prune
    new_internal_node->sax_ = leaf_node->sax_;
    new_internal_node->card_ = leaf_node->card_;

    SAX node_isax;
    CARD node_card;
    for (int j = 0; j < left_leaf_node->len; j ++ ) {
        SAX key_sax = left_leaf_node->leaf_keys[j].sax_;
        // nodeisax
        if (j == 0) {
            node_isax = key_sax;
            node_card.set_max_card();
        }
        else {
            for (int k = 0; k < SEGMENTS; k ++ ) {
                for (int l = 1; l <= node_card.card[k]; l ++ ) {
                    saxt_type bit = 1 << (BIT_CARDINALITY - l);
                    if ((node_isax.sax[k] & bit) != (key_sax.sax[k] & bit)) {
                        node_card.card[k] = l - 1;
                        node_isax.sax[k] &= bit;
                    }
                }
            }
        }
    }
    left_leaf_node->sax_ = node_isax;
    left_leaf_node->card_ = node_card;
    for (int j = 0; j < right_leaf_node->len; j ++ ) {
        SAX key_sax = right_leaf_node->leaf_keys[j].sax_;
        // nodeisax
        if (j == 0) {
            node_isax = key_sax;
            node_card.set_max_card();
        }
        else {
            for (int k = 0; k < SEGMENTS; k ++ ) {
                for (int l = 1; l <= node_card.card[k]; l ++ ) {
                    saxt_type bit = 1 << (BIT_CARDINALITY - l);
                    if ((node_isax.sax[k] & bit) != (key_sax.sax[k] & bit)) {
                        node_card.card[k] = l - 1;
                        node_isax.sax[k] &= bit;
                    }
                }
            }
        }
    }
    right_leaf_node->sax_ = node_isax;
    right_leaf_node->card_ = node_card;
#endif

    if (pre_is_root) {  // root，rootmapleaf_nodeinternal_node
#if binary_tree_root_full
        RootNode* root_node = (RootNode*) pre_node;
        u_int64_t sax_first_bit = 0;
        saxt_type first_bit = 1 << (BIT_CARDINALITY - 1);
        for (int i = 0; i < SEGMENTS; i ++ ) {
            sax_first_bit |= ((insert_sax->sax[i] & first_bit) >> (BIT_CARDINALITY - 1)) << i;
        }
        root_node->children[sax_first_bit].first = new_internal_node;
        root_node->children[sax_first_bit].second = false;
#else
        RootNode* root_node = (RootNode*) pre_node;
        root_node->node = new_internal_node;
        root_node->root_is_leaf = false;
#endif
    }
    else {
        InternalNode* internal_node = (InternalNode*) pre_node;
        if (internal_node->left == leaf_node) { // leaf_nodepre_node
            internal_node->is_left_leaf = false;
            internal_node->left = new_internal_node;
        }
        else if (internal_node->right== leaf_node){
            internal_node->is_right_leaf = false;
            internal_node->right = new_internal_node;
        }
        else {
            cout << "pre node " << endl;
            exit(-1);
        }
    }
    delete leaf_node;

    pre_node = new_internal_node;
    pre_is_root = false;

//    cout << "split finish" << endl;

    if(left_leaf_node->len >= LEAF_MAX_NUM || right_leaf_node->len >= LEAF_MAX_NUM) {
        cout << "" << endl;
//        std::cout << ":" << left_leaf_node->len << " :" << right_leaf_node->len << std::endl;
//        for (int i = 0; i < right_leaf_node->len; i ++ ) {
//            sax_print(right_leaf_node->leaf_keys[i].sax_);
//        }
//        cout << endl;
        exit(233);
    }

    if (insert_sax->sax[s] <= split_sax) { // 
        return left_leaf_node;
    }
    else {  // 
        return right_leaf_node;
    }

}

void InsertNode(const SAX* insert_sax, const u_int64_t p, void* pre_node, bool pre_is_root,
                void* node, bool is_leaf) {

    if (!is_leaf) {

        InternalNode* internal_node = (InternalNode*) node;

        u_int8_t s = internal_node->split_segment;
        if (insert_sax->sax[s] <= internal_node->split_sax) {    // 
            InsertNode(insert_sax, p, internal_node, false, internal_node->left, internal_node->is_left_leaf);
        }
        else {  // 
            InsertNode(insert_sax, p, internal_node, false, internal_node->right, internal_node->is_right_leaf);
        }
    }
    else {

        LeafNode* leaf_node = (LeafNode*) node;

        if (leaf_node->len >= LEAF_MAX_NUM) {
//            cout << "start split" << endl;
            leaf_node = SplitNode(leaf_node, pre_node, pre_is_root, insert_sax);
//            cout << "split finish" << endl;
        }

        assert(leaf_node->len < LEAF_MAX_NUM);
        LeafKey* leaf_key = &leaf_node->leaf_keys[leaf_node->len];
//        memcpy(&leaf_key->sax_, insert_sax, sizeof(SAX));
        leaf_key->sax_ = *insert_sax;
        leaf_key->p = p;
        leaf_node->len ++ ;
        for (int i = 0; i < SEGMENTS; i ++ ) {
            leaf_node->sax_lb.sax[i] = min(leaf_node->sax_lb.sax[i], insert_sax->sax[i]);
            leaf_node->sax_ub.sax[i] = max(leaf_node->sax_ub.sax[i], insert_sax->sax[i]);
        }
    }

}


void bSAX2Tree::Insert(const SAX* insert_sax, const u_int64_t p) {
#if binary_tree_root_full
    // mapkey
    u_int64_t sax_first_bit = 0;
    saxt_type first_bit = 1 << (BIT_CARDINALITY - 1);
    for (int i = 0; i < SEGMENTS; i ++ ) {
        sax_first_bit |= ((insert_sax->sax[i] & first_bit) >> (BIT_CARDINALITY - 1)) << i;
    }

    std::pair<void*, bool> node_pair = root->children[sax_first_bit];
    InsertNode(insert_sax, p, root, true, node_pair.first, node_pair.second);
#else
    InsertNode(insert_sax, p, root, true, root->node, root->root_is_leaf);
#endif
}

/**
 * 
 */
void GetTopKAns(const char* filename, int k, vector<LeafKey>& leaf_ans, TS* search_ts, ts_type* search_paa,
                vector<TS*>& ts_ans, vector<float>& top_dis) {
    vector<DisP> dis_p(leaf_ans.size());

    for (int i = 0; i < leaf_ans.size(); i ++) {
        float low_dis = min_dist_paa_to_sax(search_paa, leaf_ans[i].sax_);
        dis_p[i] = {low_dis, leaf_ans[i].p};
    }

    FILE * file;
    file = fopen(filename, "r");
    if (!file) {
        cout << "" << filename << endl;
        exit(-1);
    }
    TopKHeap* heap = new TopKHeap(k);
    TS* ts_arr = new TS[K+1];
    TS* read_ts = ts_arr;


#if sort_strategy == 0
    sort(dis_p.begin(), dis_p.end(), PCmp);
    for(int i = 0; i < dis_p.size(); i ++) {

        if (!heap->check_approximate(dis_p[i].dis)) continue;


        long offset = dis_p[i].p * TS_LENGTH * sizeof(ts_type);
        fseek(file, offset, SEEK_SET);
        size_t bytes_read = fread(read_ts->ts, sizeof(TS), 1, file);
        COUNT_APPROXIMATE_READ_TS(1)
        float true_dis = ts_euclidean_distance(search_ts->ts, read_ts->ts);
        heap->push_ans_approximate(true_dis, dis_p[i].p);
    }
#endif


#if sort_strategy >= 1
    sort(dis_p.begin(), dis_p.end(), DisCmp);
    for(int i = 0; i < dis_p.size(); i ++) {

        if (!heap->check_approximate(dis_p[i].dis)) break;


        long offset = dis_p[i].p * TS_LENGTH * sizeof(ts_type);
        fseek(file, offset, SEEK_SET);
        size_t bytes_read = fread(read_ts->ts, sizeof(TS), 1, file);
        COUNT_APPROXIMATE_READ_TS(1)
        float true_dis = ts_euclidean_distance(search_ts->ts, read_ts->ts);
        heap->push_ans_approximate(true_dis, dis_p[i].p);
    }

#endif

#if sort_strategy == 3
    sort(dis_p.begin(), dis_p.end(), DisCmp);

    int b_size = dis_p.size() / sort_batch_num;

    int num = 0;
    if (b_size != 0) num = dis_p.size() / b_size;


    int last = dis_p.size() - num * b_size;

    bool is_break = false;

    // 
    for (int i = 0; i < num; i ++) {
        cout<<""<<i<<""<<endl;
        sort(dis_p.begin() + b_size * i, dis_p.begin() + b_size * (i + 1), PCmp);
        for (int j = b_size * i; j < b_size * (i + 1); j ++) {
            if (heap->check_approximate(dis_p[j].dis)) {
                is_break = true;
                continue;
            }
            long offset = dis_p[j].p * TS_LENGTH * sizeof(ts_type);
            fseek(file, offset, SEEK_SET);
            size_t bytes_read = fread(read_ts->ts, sizeof(TS), 1, file);
            COUNT_EXACT_READ_TS(1)
            float true_dis = ts_euclidean_distance(search_ts->ts, read_ts->ts);
            heap->push_ans_approximate(true_dis, &read_ts);
        }
        if (is_break) break;
    }

    if (!is_break && last > 0) {
        sort(dis_p.begin() + b_size * num, dis_p.end(), PCmp);
        for (int j = b_size * num; j < dis_p.size(); j ++) {
            if (heap->check_approximate(dis_p[j].dis)) continue;


            long offset = dis_p[j].p * TS_LENGTH * sizeof(ts_type);
            fseek(file, offset, SEEK_SET);
            size_t bytes_read = fread(read_ts->ts, sizeof(TS), 1, file);
            COUNT_APPROXIMATE_READ_TS(1)
            float true_dis = ts_euclidean_distance(search_ts->ts, read_ts->ts);
            heap->push_ans_approximate(true_dis, &read_ts);
        }
    }
#endif

    while(!heap->pq.empty()) {
        top_dis.emplace_back(heap->pq.top().first);
        ts_ans.emplace_back(heap->pq.top().second);
        heap->pq.pop();
    }




    delete heap;
    delete[] ts_arr;
}
/**
 * ，bsf
 */
void GetTopKAns(const char* filename, int k, vector<LeafKey>& leaf_ans, TS* search_ts, ts_type* search_paa,
                vector<TS*>& ts_ans, vector<float>& top_dis, vector<float> bsfs, vector<TS*> bsf_p) {

    TopKHeap* heap = new TopKHeap(k);
    unordered_set<TS*> found_p;

    for(int i=0;i<bsfs.size();i++) {
        heap->pq.push({bsfs[i], bsf_p[i]});
        found_p.insert(bsf_p[i]);
    }
    float bsf = heap->pq.top().first;
    cout<<"："<<bsf<<endl;



    vector<DisP> dis_p;
    COMPUTE_MIN_DIS_START
    for (auto & leaf_an : leaf_ans) {
        if(found_p.count((TS*)leaf_an.p)) continue;
        float low_dis = min_dist_paa_to_sax(search_paa, leaf_an.sax_);
        if (low_dis < bsf) dis_p.emplace_back(low_dis, leaf_an.p);
    }
    COMPUTE_MIN_DIS_END
    cout<<"dis_p"<<dis_p.size()<<endl;
    FILE * file;
    file = fopen(filename, "r");
    if (!file) {
        cout << "" << filename << endl;
        exit(-1);
    }

    TS* read_ts = new TS;


#if sort_strategy == 0
    sort(dis_p.begin(), dis_p.end(), PCmp);

    for(auto & i : dis_p) {

//        if(i % 10000 == 0) cout<< (u_int64_t)dis_p[i].p <<endl;

        if (i.dis <= bsf) {

//        TS* read_ts = new TS;
        fseek(file, i.p * TS_LENGTH * sizeof(ts_type), SEEK_SET);
        fread(read_ts->ts, sizeof(ts_type), TS_LENGTH, file);
        COUNT_EXACT_READ_TS(1)
        float true_dis = ts_euclidean_distance(search_ts->ts, read_ts->ts);
        heap->push_ans_exact(true_dis, i.p, bsf);

        }
    }

#endif

#if sort_strategy == 1
    sort(dis_p.begin(), dis_p.end(), DisCmp);

    for(int i = 0; i < dis_p.size(); i ++) {
        if (dis_p[i].dis > bsf) break;


        long offset = dis_p[i].p * TS_LENGTH * sizeof(ts_type);
        fseek(file, offset, SEEK_SET);
        size_t bytes_read = fread(read_ts->ts, sizeof(TS), 1, file);
        COUNT_EXACT_READ_TS(1)
        float true_dis = ts_euclidean_distance(search_ts->ts, read_ts->ts);
        heap->push_ans_exact(true_dis, dis_p[i].p, bsf);
    }
#endif


#if sort_strategy == 2
    sort(dis_p.begin(), dis_p.end(), DisCmp);


    //
    int b_size = dis_p.size() / sort_batch_num;


    int num = 0;
    if (b_size != 0) num = dis_p.size() / b_size;


    int last = dis_p.size() - num * b_size;

    bool is_break = false;

    // 
    for (int i = 0; i < num; i ++) {
        cout<<""<<i<<""<<endl;
        sort(dis_p.begin() + b_size * i, dis_p.begin() + b_size * (i + 1), PCmp);
        for (int j = b_size * i; j < b_size * (i + 1); j ++) {
            if (dis_p[j].dis > bsf) {
                is_break = true;
                continue;
            }


            long offset = dis_p[j].p * TS_LENGTH * sizeof(ts_type);
            fseek(file, offset, SEEK_SET);
            size_t bytes_read = fread(read_ts->ts, sizeof(TS), 1, file);
            COUNT_EXACT_READ_TS(1)
            float true_dis = ts_euclidean_distance(search_ts->ts, read_ts->ts);
            heap->push_ans_exact(true_dis, dis_p[j].p, bsf);
        }
        if (is_break) break;
    }

    if (!is_break && last > 0) {
        sort(dis_p.begin() + b_size * num, dis_p.end(), PCmp);
        for (int j = b_size * num; j < dis_p.size(); j ++) {
            if (dis_p[j].dis > bsf) continue;


            long offset = dis_p[j].p * TS_LENGTH * sizeof(ts_type);
            fseek(file, offset, SEEK_SET);
            size_t bytes_read = fread(read_ts->ts, sizeof(TS), 1, file);
            COUNT_EXACT_READ_TS(1)
            float true_dis = ts_euclidean_distance(search_ts->ts, read_ts->ts);
            heap->push_ans_exact(true_dis, dis_p[j].p, bsf);
        }
    }

#endif
//    cout<<"ppp:"<<endl;
//    vector<u_int64_t> ppp;
//    for(int i=0;i<100;i++) {
//        ppp.push_back((u_int64_t)heap->pq.top().second);
//        heap->pq.pop();
//    }

//    sort(ppp.begin(), ppp.end());
//    for(auto item: ppp) {
//        cout<<item<<endl;
//    }

    cout<<"p"<<(u_int64_t)heap->pq.top().second<<endl;

    delete heap;
//    delete[] ts_arr;
    delete read_ts;
}

void GetTopKAnsM(const char* filename, int k, vector<LeafNode*>& leaf_ans, TS* search_ts, ts_type* search_paa,
                vector<TS*>& ts_ans, vector<float>& top_dis, vector<float> bsfs, vector<TS*> bsf_p) {

    TopKHeap* heap = new TopKHeap(k);
    unordered_set<TS*> found_p;

    for(int i=0;i<bsfs.size();i++) {
        heap->pq.push({bsfs[i], bsf_p[i]});
        found_p.insert(bsf_p[i]);
    }
    float bsf = heap->pq.top().first;
    cout<<"："<<bsf<<endl;


    COMPUTE_MIN_DIS_START
    vector<DisP> dis_p;
    vector<DisP> dis_p_tmp;
    vector<DisNode> dis_node;
    for(auto & leaf_an : leaf_ans) {

        DisNode tmp_node(MAXFLOAT, dis_p_tmp.size(), 0);
        COUNT_EXACT_ANS(leaf_an->len)

        for(int i=0;i<leaf_an->len;i++) {
            LeafKey& leafKey = leaf_an->leaf_keys[i];
            if(found_p.count((TS*)leafKey.p)) continue;
            float low_dis = min_dist_paa_to_sax(search_paa, leafKey.sax_);
            if (low_dis < bsf) {
                dis_p_tmp.emplace_back(low_dis, leafKey.p);
                tmp_node.size++;
                tmp_node.dis = min(low_dis, tmp_node.dis);
            }
        }

        if (tmp_node.size != 0) dis_node.push_back(tmp_node);
    }

    dis_p.resize(dis_p_tmp.size());
    sort(dis_node.begin(), dis_node.end(), DisNodeCmp);

    u_int64_t now_size = 0;
    for(auto & a_node: dis_node) {
        memcpy(&dis_p[now_size], &dis_p_tmp[a_node.begin], a_node.size * sizeof(DisP));
        now_size += a_node.size;
    }

    COMPUTE_MIN_DIS_END

    cout<<"dis_p"<<dis_p.size()<<endl;
    FILE * file;
    file = fopen(filename, "r");
    if (!file) {
        cout << "" << filename << endl;
        exit(-1);
    }

    TS* read_ts = new TS;


#if sort_strategy == 1

    for(auto & i : dis_p) {

//        if(i % 10000 == 0) cout<< (u_int64_t)dis_p[i].p <<endl;

        if (i.dis <= bsf) {

//        TS* read_ts = new TS;
            fseek(file, i.p * TS_LENGTH * sizeof(ts_type), SEEK_SET);
            fread(read_ts->ts, sizeof(ts_type), TS_LENGTH, file);
            COUNT_EXACT_READ_TS(1)
            float true_dis = ts_euclidean_distance(search_ts->ts, read_ts->ts);
            heap->push_ans_exact(true_dis, i.p, bsf);

        }
    }
#elif sort_strategy == 2

    int b_size = dis_p.size() / sort_batch_num;


    int num = 0;
    if (b_size != 0) num = dis_p.size() / b_size;


    int last = dis_p.size() - num * b_size;



    // 
    for (int i = 0; i < num; i ++) {
        cout<<""<<i<<""<<endl;
        sort(dis_p.begin() + b_size * i, dis_p.begin() + b_size * (i + 1), PCmp);
        for (int j = b_size * i; j < b_size * (i + 1); j ++) {
            if (dis_p[j].dis > bsf) {
                continue;
            }


            long offset = dis_p[j].p * TS_LENGTH * sizeof(ts_type);
            fseek(file, offset, SEEK_SET);
            size_t bytes_read = fread(read_ts->ts, sizeof(TS), 1, file);
            COUNT_EXACT_READ_TS(1)
            float true_dis = ts_euclidean_distance(search_ts->ts, read_ts->ts);
            heap->push_ans_exact(true_dis, dis_p[j].p, bsf);
        }
    }

    if (last > 0) {
        sort(dis_p.begin() + b_size * num, dis_p.end(), PCmp);
        for (int j = b_size * num; j < dis_p.size(); j ++) {
            if (dis_p[j].dis > bsf) continue;


            long offset = dis_p[j].p * TS_LENGTH * sizeof(ts_type);
            fseek(file, offset, SEEK_SET);
            size_t bytes_read = fread(read_ts->ts, sizeof(TS), 1, file);
            COUNT_EXACT_READ_TS(1)
            float true_dis = ts_euclidean_distance(search_ts->ts, read_ts->ts);
            heap->push_ans_exact(true_dis, dis_p[j].p, bsf);
        }
    }

#endif

    delete heap;
//    delete[] ts_arr;
    delete read_ts;
}

/**
 * 
 * @param search_sax
 * @param node
 * @param is_leaf
 * @param leaf_ans
 * @param card_now
 */
void bSAX2Tree::DFS(const SAX* search_sax, void* node, bool is_leaf, vector<LeafKey>& leaf_ans) const {
    if (!is_leaf) {
        InternalNode* internal_node = (InternalNode*) node;
        u_int8_t s = internal_node->split_segment;
//        cout << ":" << (int)internal_node->split_segment << " sax:" << (int) internal_node->split_sax << endl;
        if (search_sax->sax[s] <= internal_node->split_sax) {    // 
            DFS(search_sax, internal_node->left, internal_node->is_left_leaf, leaf_ans);
            if (leaf_ans.size() < num_approximate_search_key) { // 
                DFS(search_sax, internal_node->right, internal_node->is_right_leaf, leaf_ans);
            }
        }
        else {  // 
            DFS(search_sax, internal_node->right, internal_node->is_right_leaf, leaf_ans);
            if (leaf_ans.size() < num_approximate_search_key) { // 
                DFS(search_sax, internal_node->left, internal_node->is_left_leaf, leaf_ans);
            }
        }
    }
    else {
        LeafNode* leaf_node = (LeafNode*) node;

        for (int i = 0; i < leaf_node->len; i ++ ) {
            leaf_ans.emplace_back(leaf_node->leaf_keys[i]);
            if (leaf_ans.size() >= num_approximate_search_key) {
                break;
            }
        }
    }
}

void bSAX2Tree::ApproximateSearch(int k, TS* search_ts, vector<TS*>& ts_ans, vector<float>& top_dis) const {
    APPROXIMATE_TIME_START
    SAX search_sax;
    sax_from_ts(search_ts->ts, search_sax.sax);

//    sax_print(search_sax);
//    sax_print_bit(search_sax);

    vector<LeafKey> leaf_ans;
    APPROXIMATE_INDEX_START
#if binary_tree_root_full
    // mapkey
    u_int64_t sax_first_bit = 0;
    saxt_type first_bit = 1 << (BIT_CARDINALITY - 1);
    for (int i = 0; i < SEGMENTS; i ++ ) {
        sax_first_bit |= ((search_sax.sax[i] & first_bit) >> (BIT_CARDINALITY - 1)) << i;
    }

    std::pair<void*, bool> node_pair = root->children[sax_first_bit];

    DFS(&search_sax, node_pair.first, node_pair.second, leaf_ans);
#else
    DFS(&search_sax, root->node, root->root_is_leaf, leaf_ans);
#endif
    APPROXIMATE_INDEX_END

    COUNT_APPROXIMATE_ANS(leaf_ans.size())

    ts_type search_paa[SEGMENTS];
    paa_from_ts(search_ts->ts, search_paa);

    APPROXIMATE_GET_ANS_START
    GetTopKAns(filename, k, leaf_ans, search_ts, search_paa, ts_ans, top_dis);
    APPROXIMATE_GET_ANS_END

    APPROXIMATE_TIME_END
}


/**
 * ，leaf_ans
 * @param search_paa
 * @param root
 * @param leaf_ans
 * @param top_dis top_dis
 */
void bSAX2Tree::BFS(ts_type* search_paa, float top_dis, vector<LeafKey>& leaf_ans) const {
    queue<pair<void*, bool>> q; // node, is_leaf
#if binary_tree_root_full
    for (auto& node_pair: root->children) {
        q.emplace(node_pair);
    }
#else
    q.push({root->node, root->root_is_leaf});
#endif
    while(!q.empty()) {
        pair<void*, bool> node_pair = q.front();
        q.pop();
        if (!node_pair.second) {    // 
            InternalNode* internal_node = (InternalNode*) node_pair.first;
#if exp3_b_binary_use_isax_to_prune
            float dis = min_dist_paa_to_isax(search_paa, internal_node->sax_, internal_node->card_);
#else
            float dis = min_dist_paa_to_bsax(search_paa, internal_node->sax_lb, internal_node->sax_ub);
#endif
            COUNT_NODE_COMPUTE_MIN_DIS(1)
//            cout << "dis:" << dis << " top dis:" << top_dis << endl;

            if (dis >= top_dis) continue;

            pair<void*, bool> left_p = {internal_node->left, internal_node->is_left_leaf};
            q.emplace(left_p);
            pair<void*, bool> right_p = {internal_node->right, internal_node->is_right_leaf};
            q.emplace(right_p);
        }
        else {  // 
            LeafNode* leaf_node = (LeafNode*) node_pair.first;
#if exp3_b_binary_use_isax_to_prune
            float dis = min_dist_paa_to_isax(search_paa, leaf_node->sax_, leaf_node->card_);
#else
            float dis = min_dist_paa_to_bsax(search_paa, leaf_node->sax_lb, leaf_node->sax_ub);
#endif
            COUNT_NODE_COMPUTE_MIN_DIS(1)
            if (dis >= top_dis) continue;
            for (int i = 0; i < leaf_node->len; i ++ ) {
                leaf_ans.emplace_back(leaf_node->leaf_keys[i]);
            }
        }
    }
}

void bSAX2Tree::BFSM(ts_type* search_paa, float top_dis, vector<LeafNode*>& leaf_ans) const {
    queue<pair<void*, bool>> q; // node, is_leaf
#if binary_tree_root_full
    for (auto& node_pair: root->children) {
        q.emplace(node_pair);
    }
#else
    q.push({root->node, root->root_is_leaf});
#endif
    while(!q.empty()) {
        pair<void*, bool> node_pair = q.front();
        q.pop();
        if (!node_pair.second) {    // 
            InternalNode* internal_node = (InternalNode*) node_pair.first;
#if exp3_b_binary_use_isax_to_prune
            float dis = min_dist_paa_to_isax(search_paa, internal_node->sax_, internal_node->card_);
#else
            float dis = min_dist_paa_to_bsax(search_paa, internal_node->sax_lb, internal_node->sax_ub);
#endif
            COUNT_NODE_COMPUTE_MIN_DIS(1)
//            cout << "dis:" << dis << " top dis:" << top_dis << endl;

            if (dis >= top_dis) continue;

            pair<void*, bool> left_p = {internal_node->left, internal_node->is_left_leaf};
            q.emplace(left_p);
            pair<void*, bool> right_p = {internal_node->right, internal_node->is_right_leaf};
            q.emplace(right_p);
        }
        else {  // 
            LeafNode* leaf_node = (LeafNode*) node_pair.first;
#if exp3_b_binary_use_isax_to_prune
            float dis = min_dist_paa_to_isax(search_paa, leaf_node->sax_, leaf_node->card_);
#else
            float dis = min_dist_paa_to_bsax(search_paa, leaf_node->sax_lb, leaf_node->sax_ub);
#endif
            COUNT_NODE_COMPUTE_MIN_DIS(1)
            if (dis >= top_dis) continue;
            leaf_ans.emplace_back(leaf_node);
        }
    }
}

void bSAX2Tree::ExactSearch(int k, TS* search_ts, vector<TS*>& ts_ans, vector<float>& top_dis) const {
    EXACT_TIME_START
    vector<TS*> approximate_ts_ans;
    vector<float> approximate_top_dis;
    ApproximateSearch(k, search_ts, approximate_ts_ans, approximate_top_dis);

    float bsf = approximate_top_dis[0];


    ts_type search_paa[SEGMENTS];
    paa_from_ts(search_ts->ts, search_paa);

    vector<LeafKey> leaf_ans;
    EXACT_INDEX_START
    BFS(search_paa, bsf, leaf_ans);
    EXACT_INDEX_END

    COUNT_EXACT_ANS(leaf_ans.size())

    EXACT_GET_ANS_START
    GetTopKAns(filename, k, leaf_ans, search_ts, search_paa, ts_ans, top_dis, approximate_top_dis, approximate_ts_ans);
    EXACT_GET_ANS_END
    EXACT_TIME_END


    cout << "bsf:" << bsf << endl;
    PRINT_COUNT
    COUNT_CLEAR
}

void bSAX2Tree::ExactSearchM(int k, TS* search_ts, vector<TS*>& ts_ans, vector<float>& top_dis) const {
    EXACT_TIME_START
    vector<TS*> approximate_ts_ans;
    vector<float> approximate_top_dis;
    ApproximateSearch(k, search_ts, approximate_ts_ans, approximate_top_dis);

    float bsf = approximate_top_dis[0];


    ts_type search_paa[SEGMENTS];
    paa_from_ts(search_ts->ts, search_paa);

    vector<LeafNode*> leaf_ans;
    EXACT_INDEX_START
    BFSM(search_paa, bsf, leaf_ans);
    EXACT_INDEX_END



    EXACT_GET_ANS_START
    GetTopKAnsM(filename, k, leaf_ans, search_ts, search_paa, ts_ans, top_dis, approximate_top_dis, approximate_ts_ans);
    EXACT_GET_ANS_END
    EXACT_TIME_END


    cout << "bsf:" << bsf << endl;
    PRINT_COUNT
    COUNT_CLEAR
}

//bool cmpp(pair<SAX, u_int64_t>& a, pair<SAX, u_int64_t>& b) {
//    return a.second < b.second;
//}


void bSAX2Tree::ExactSearchExp3(TS* search_ts, float bsf, vector<pair<SAX, u_int64_t>>& sax_and_p, vector<float>& top_dis) const {
    ts_type search_paa[SEGMENTS];
    paa_from_ts(search_ts->ts, search_paa);

    vector<LeafKey> leaf_ans;

    EXACT_INDEX_START
    BFS(search_paa, bsf, leaf_ans);
    EXACT_INDEX_END

    for (auto& leaf_key: leaf_ans) {
        sax_and_p.emplace_back(leaf_key.sax_, leaf_key.p);
    }
    COUNT_EXACT_ANS(sax_and_p.size())

    cout << "bsf:" << bsf << endl;
    PRINT_COUNT
    COUNT_CLEAR


//    sort(sax_and_p.begin(), sax_and_p.end(), cmpp);
//    for(int i=0;i<sax_and_p.size();i++) {
//        if(i%10000 == 0)
//        cout<<(u_int64_t)sax_and_p[i].second<<endl;
//    }

}

void bSAX2Tree::materialized(const char *_filename) {
    FILE * m_data_file = fopen(_filename,"w");
    if (!m_data_file) {
        cout << "" << _filename << endl;
        exit(-1);
    }

    FILE * data_file = fopen(filename,"r");
    if (!data_file) {
        cout << "" << filename << endl;
        exit(-1);
    }

    u_int64_t tmp_ts_cap = 100000;
    TS* tmp_ts = new TS[tmp_ts_cap];
#if is_reorder_m
    TS* tmp_ts1 = new TS[tmp_ts_cap];
#endif
    u_int64_t tmp_ts_size = 0;
    u_int64_t new_p = 0;
    u_int64_t read_batch_size = (u_int64_t) build_batch_size * 1000 *1000 * 1024 / 4 / TS_LENGTH;
    TS* tmp_read_ts = new TS[read_batch_size];

    for(int j=0;j< read_batch_size && j < TOTAL_TS;j++) {
        fread(tmp_read_ts + j, sizeof(TS), 1, data_file);
    }

    queue<pair<void*, bool>> q; // node, is_leaf
#if binary_tree_root_full
    for (auto& node_pair: root->children) {
        q.emplace(node_pair);
    }
#else
    q.push({root->node, root->root_is_leaf});
#endif
    while(!q.empty()) {
        pair<void*, bool> node_pair = q.front();
        q.pop();
        if (!node_pair.second) {    // 
            InternalNode* internal_node = (InternalNode*) node_pair.first;
            pair<void*, bool> left_p = {internal_node->left, internal_node->is_left_leaf};
            q.emplace(left_p);
            pair<void*, bool> right_p = {internal_node->right, internal_node->is_right_leaf};
            q.emplace(right_p);
        }
        else {  // 
            LeafNode* leaf_node = (LeafNode*) node_pair.first;
            for (int i = 0; i < leaf_node->len; i ++ ) {
                if (leaf_node->leaf_keys[i].p < read_batch_size) {
#if is_reorder_m
                    if (new_p < tmp_ts_cap) {
                        tmp_ts1[tmp_ts_size] = tmp_read_ts[leaf_node->leaf_keys[i].p];
                    } else {
                        tmp_ts[tmp_ts_size] = tmp_read_ts[leaf_node->leaf_keys[i].p];
                    }
#else
                    tmp_ts[tmp_ts_size] = tmp_read_ts[leaf_node->leaf_keys[i].p];
#endif
                }
                tmp_ts_size++;
                if(tmp_ts_size == tmp_ts_cap) {
#if is_reorder_m
                    if (new_p < tmp_ts_cap) {
                    } else if (new_p == 4*tmp_ts_cap-1) {
                        fwrite(tmp_ts, sizeof(TS), tmp_ts_cap, m_data_file);
                        fwrite(tmp_ts1, sizeof(TS), tmp_ts_cap, m_data_file);
                    }
                    else {
                        fwrite(tmp_ts, sizeof(TS), tmp_ts_cap, m_data_file);
                    }
#else
                    cout<<new_p<<endl;
                    fwrite(tmp_ts, sizeof(TS), 1000000, m_data_file);
#endif
                    tmp_ts_size = 0;
                }
                new_p++;
            }
        }
    }

    if(tmp_ts_size != 0) {
        fwrite(tmp_ts, sizeof(TS), tmp_ts_size, m_data_file);
    }
    delete[] tmp_ts;


    for (u_int64_t kk = read_batch_size; kk < TOTAL_TS; kk += read_batch_size) {
        for(int j=0;j< read_batch_size && kk + j < TOTAL_TS;j++) {
            fread(tmp_read_ts + j, sizeof(TS), 1, data_file);
        }
        cout<<0<<endl;
        new_p = 0;
#if binary_tree_root_full
        for (auto& node_pair: root->children) {
        q.emplace(node_pair);
    }
#else
        q.push({root->node, root->root_is_leaf});
#endif
        while(!q.empty()) {
            pair<void*, bool> node_pair = q.front();
            q.pop();
            if (!node_pair.second) {    // 
                InternalNode* internal_node = (InternalNode*) node_pair.first;
                pair<void*, bool> left_p = {internal_node->left, internal_node->is_left_leaf};
                q.emplace(left_p);
                pair<void*, bool> right_p = {internal_node->right, internal_node->is_right_leaf};
                q.emplace(right_p);
            }
            else {  // 
                LeafNode* leaf_node = (LeafNode*) node_pair.first;
                for (int i = 0; i < leaf_node->len; i ++ ) {
                    if (leaf_node->leaf_keys[i].p < kk + read_batch_size && leaf_node->leaf_keys[i].p >= kk) {
                        fseek(m_data_file, new_p * TS_LENGTH * sizeof(ts_type), SEEK_SET);
                        fwrite(&tmp_read_ts[leaf_node->leaf_keys[i].p-kk], sizeof(TS), 1, m_data_file);
                    }
                    new_p++;
                }
            }
        }
    }


}

void bSAX2Tree::materialized() {
    u_int64_t new_p = 0;

    queue<pair<void*, bool>> q; // node, is_leaf
#if binary_tree_root_full
    for (auto& node_pair: root->children) {
        q.emplace(node_pair);
    }
#else
    q.push({root->node, root->root_is_leaf});
#endif
    while(!q.empty()) {
        pair<void*, bool> node_pair = q.front();
        q.pop();
        if (!node_pair.second) {    // 
            InternalNode* internal_node = (InternalNode*) node_pair.first;
            pair<void*, bool> left_p = {internal_node->left, internal_node->is_left_leaf};
            q.emplace(left_p);
            pair<void*, bool> right_p = {internal_node->right, internal_node->is_right_leaf};
            q.emplace(right_p);
        }
        else {  // 
            LeafNode* leaf_node = (LeafNode*) node_pair.first;
            for (int i = 0; i < leaf_node->len; i ++ ) {
#if is_reorder_m
                if (new_p < 100000) {
                    leaf_node->leaf_keys[i].p = new_p + 3 * 100000;
                }
                else if (new_p < 400000){
                    leaf_node->leaf_keys[i].p = new_p - 100000;
                }
                else {
                    leaf_node->leaf_keys[i].p = new_p;
                }
                new_p++;
#else
                leaf_node->leaf_keys[i].p = new_p++;
#endif
            }
        }
    }
}


