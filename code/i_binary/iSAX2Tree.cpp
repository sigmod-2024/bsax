


#include <vector>
#include <algorithm>
#include "iSAX2Tree.h"
#include "RootNode.h"
#include "../util/TopKHeap.h"
#include "CntRecord.h"
#include "TimeRecord.h"
#include "unordered_set"
using namespace std;

void read_ts_batch(int n, FILE* data_file, vector<SAX>& sax_vec, vector<TS>& ts_vec) {
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

void iSAX2Tree::Build() {
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

void iSAX2Tree::BuildFromSAX(const char *sax_filename) {
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



    }
    BUILD_INDEX_END

    fclose (data_file);
}



void iSAX2Tree::BuildTree() {
    root = new RootNode();
#if binary_tree_root_full
    for (u_int64_t i = 0; i < ((u_int64_t)1 << SEGMENTS); i ++ ) {
        SAX sax_;
        CARD card_;
        for (int j = 0; j < SEGMENTS; j ++ ) {
            sax_.sax[j] = ((i >> j) & 1) << (BIT_CARDINALITY - 1);
            card_.card[j] = 1;
        }
        root->children[i] = {new LeafNode(sax_, card_), true};
    }
#else
    SAX sax;
    CARD card;
    sax.set_min_value();
    card.set_min_card();
    root->node = new LeafNode(sax, card);
#endif
}

LeafNode* SplitNode(LeafNode* leaf_node, void*& pre_node, bool& pre_is_root, const SAX* insert_sax) {
    int min_s = -1, min_s_bak = -1;
    float min_dis = MAXFLOAT, min_dis_bak = MAXFLOAT;

    for (int i = 0; i < SEGMENTS; i ++ ) {
        if (leaf_node->card_.card[i] >= BIT_CARDINALITY) continue;  
        int mean = leaf_node->sum[i] / leaf_node->len;
        int stdev = std::sqrt((leaf_node->square_sum[i] / leaf_node->len) - mean * mean);





        float breakpoint = sax_a[leaf_node->sax_.sax[i] | (1 << (BIT_CARDINALITY - leaf_node->card_.card[i] - 1))];   

        if (min_dis_bak > dist_breakpoint_to_sax(breakpoint, mean)) {   
            min_dis_bak = dist_breakpoint_to_sax(breakpoint, mean);
            min_s_bak = i;
        }
        float avg_minus_var = sax_a[std::max(mean - 3 * stdev, 0)]; 
        float avg_add_var = sax_a[std::min(mean + 3 * stdev + 1, 256)];
        if (avg_minus_var <= breakpoint && avg_add_var >= breakpoint) { 
            if (min_dis > dist_breakpoint_to_sax(breakpoint, mean)) {
                min_dis = dist_breakpoint_to_sax(breakpoint, mean);
                min_s = i;
            }
        }
    }
    int s;
    if (min_s == -1) {
        s = min_s_bak;
    }
    else {
        s = min_s;
    }

    
    CARD new_card = leaf_node->card_;
    new_card.card[s] ++ ;
    SAX left_sax = leaf_node->sax_;
    SAX right_sax = leaf_node->sax_;

    
    sax_type bit = (1 << (BIT_CARDINALITY - new_card.card[s]));
    right_sax.sax[s] |= bit;

    LeafNode* left_leaf_node = new LeafNode(left_sax, new_card);
    LeafNode* right_leaf_node = new LeafNode(right_sax, new_card);

    for (int i = 0; i < leaf_node->len; i ++ ) {
        LeafKey* new_leaf_key;

        if (leaf_node->leaf_keys[i].sax_.sax[s] & bit) { 
            new_leaf_key = &right_leaf_node->leaf_keys[right_leaf_node->len];
            right_leaf_node->len ++ ;
        }
        else {  
            new_leaf_key = &left_leaf_node->leaf_keys[left_leaf_node->len];
            left_leaf_node->len ++ ;
        }
        new_leaf_key->sax_ = leaf_node->leaf_keys[i].sax_;
        new_leaf_key->p = leaf_node->leaf_keys[i].p;
    }

#if exp3_i_binary_use_bsax_to_prune
    for (int i = 0; i < left_leaf_node->len; i ++ ) {
        for (int j = 0; j < SEGMENTS; j ++ ) {
            left_leaf_node->sax_lb.sax[j] = min(left_leaf_node->sax_lb.sax[j], left_leaf_node->leaf_keys[i].sax_.sax[j]);
            left_leaf_node->sax_ub.sax[j] = max(left_leaf_node->sax_ub.sax[j], left_leaf_node->leaf_keys[i].sax_.sax[j]);
        }
    }
    for (int i = 0; i < right_leaf_node->len; i ++ ) {
        for (int j = 0; j < SEGMENTS; j ++ ) {
            right_leaf_node->sax_lb.sax[j] = min(right_leaf_node->sax_lb.sax[j], right_leaf_node->leaf_keys[i].sax_.sax[j]);
            right_leaf_node->sax_ub.sax[j] = max(right_leaf_node->sax_ub.sax[j], right_leaf_node->leaf_keys[i].sax_.sax[j]);
        }
    }
#endif



    InternalNode* new_internal_node = new InternalNode(leaf_node->sax_, leaf_node->card_, left_leaf_node, right_leaf_node, s);
#if exp3_i_binary_use_bsax_to_prune
    new_internal_node->sax_lb = leaf_node->sax_lb;
    new_internal_node->sax_ub = leaf_node->sax_ub;
#endif
    if (pre_is_root) {  
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
        if (internal_node->left == leaf_node) { 
            internal_node->is_left_leaf = false;
            internal_node->left = new_internal_node;
        }
        else {
            internal_node->is_right_leaf = false;
            internal_node->right = new_internal_node;
        }
    }
    delete leaf_node;

    pre_node = new_internal_node;
    pre_is_root = false;
    if (insert_sax->sax[s] & bit) { 
        return right_leaf_node;
    }
    else {  
        return left_leaf_node;
    }
}

void InsertNode(const SAX* insert_sax, const u_int64_t p, void* pre_node, bool pre_is_root,
                void* node, bool is_leaf, CARD* card_now) {

    if (!is_leaf) {
        InternalNode* internal_node = (InternalNode*) node;
        u_int8_t s = internal_node->split_segment;
        card_now->card[s] ++ ;

        
        sax_type bit = insert_sax->sax[s] & (1 << (BIT_CARDINALITY - card_now->card[s]));
        if (bit >> (BIT_CARDINALITY - card_now->card[s])) { 
            InsertNode(insert_sax, p, internal_node, false, internal_node->right, internal_node->is_right_leaf, card_now);
        }
        else {  
            InsertNode(insert_sax, p, internal_node, false, internal_node->left, internal_node->is_left_leaf, card_now);
        }
    }
    else {
        LeafNode* leaf_node = (LeafNode*) node;

        while (leaf_node->len >= LEAF_MAX_NUM) {
            leaf_node = SplitNode(leaf_node, pre_node, pre_is_root, insert_sax);
        }

        LeafKey* leaf_key = &leaf_node->leaf_keys[leaf_node->len];
        leaf_key->sax_ = *insert_sax;
        leaf_key->p = p;
        leaf_node->len ++ ;
        for (int i = 0; i < SEGMENTS; i ++ ) {
            leaf_node->sum[i] += insert_sax->sax[i];
            leaf_node->square_sum[i] += insert_sax->sax[i] * insert_sax->sax[i];

#if exp3_i_binary_use_bsax_to_prune
            leaf_node->sax_lb.sax[i] = min(leaf_node->sax_lb.sax[i], insert_sax->sax[i]);
            leaf_node->sax_ub.sax[i] = max(leaf_node->sax_ub.sax[i], insert_sax->sax[i]);
#endif
        }
    }
}


void iSAX2Tree::Insert(const SAX* insert_sax, const u_int64_t p) {
    CARD card_now;

#if binary_tree_root_full
    
    u_int64_t sax_first_bit = 0;
    saxt_type first_bit = 1 << (BIT_CARDINALITY - 1);
    for (int i = 0; i < SEGMENTS; i ++ ) {
        sax_first_bit |= ((insert_sax->sax[i] & first_bit) >> (BIT_CARDINALITY - 1)) << i;
        card_now.card[i] = 1;
    }

    std::pair<void*, bool> node_pair = root->children[sax_first_bit];
    InsertNode(insert_sax, p, root, true, node_pair.first, node_pair.second, &card_now);
#else
    card_now.set_min_card();
    InsertNode(insert_sax, p, root, true, root->node, root->root_is_leaf, &card_now);
#endif
}


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



        if (i.dis <= bsf) {


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


    
    int b_size = dis_p.size() / sort_batch_num;


    int num = 0;
    if (b_size != 0) num = dis_p.size() / b_size;


    int last = dis_p.size() - num * b_size;

    bool is_break = false;

    
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













    delete heap;

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



        if (i.dis <= bsf) {


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

    delete read_ts;
}


void iSAX2Tree::DFS(const SAX* search_sax, void* node, bool is_leaf, vector<LeafKey>& leaf_ans, CARD* card_now) const {
    if (!is_leaf) {
        InternalNode* internal_node = (InternalNode*) node;
        u_int8_t s = internal_node->split_segment;



        card_now->card[s] ++ ;

        
        sax_type bit = search_sax->sax[s] & (1 << (BIT_CARDINALITY - card_now->card[s]));




        if (bit >> (BIT_CARDINALITY - card_now->card[s])) { 
            DFS(search_sax, internal_node->right, internal_node->is_right_leaf, leaf_ans, card_now);
            if (leaf_ans.size() < num_approximate_search_key) { 
                DFS(search_sax, internal_node->left, internal_node->is_left_leaf, leaf_ans, card_now);
            }
        }
        else {  
            DFS(search_sax, internal_node->left, internal_node->is_left_leaf, leaf_ans, card_now);
            if (leaf_ans.size() < num_approximate_search_key) { 
                DFS(search_sax, internal_node->right, internal_node->is_right_leaf, leaf_ans, card_now);
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

void iSAX2Tree::ApproximateSearch(int k, TS* search_ts, vector<TS*>& ts_ans, vector<float>& top_dis) const {
    APPROXIMATE_TIME_START
    SAX search_sax;
    sax_from_ts(search_ts->ts, search_sax.sax);

    vector<LeafKey> leaf_ans;
    CARD card_now;

    APPROXIMATE_INDEX_START
#if binary_tree_root_full
    
    u_int64_t sax_first_bit = 0;
    saxt_type first_bit = 1 << (BIT_CARDINALITY - 1);
    for (int i = 0; i < SEGMENTS; i ++ ) {
        sax_first_bit |= ((search_sax.sax[i] & first_bit) >> (BIT_CARDINALITY - 1)) << i;
        card_now.card[i] = 1;
    }
    std::pair<void*, bool> node_pair = root->children[sax_first_bit];

    DFS(&search_sax, node_pair.first, node_pair.second, leaf_ans, &card_now);
#else
    card_now.set_min_card();
    DFS(&search_sax, root->node, root->root_is_leaf, leaf_ans, &card_now);
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



void iSAX2Tree::BFS(ts_type* search_paa, float top_dis, vector<LeafKey>& leaf_ans) const {
    queue<pair<void*, bool>> q; 
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
        if (!node_pair.second) {    
            InternalNode* internal_node = (InternalNode*) node_pair.first;
#if exp3_i_binary_use_bsax_to_prune
            float dis = min_dist_paa_to_bsax(search_paa, internal_node->sax_lb, internal_node->sax_ub);
#else
            float dis = min_dist_paa_to_isax(search_paa, internal_node->sax_, internal_node->card_);
#endif
            COUNT_NODE_COMPUTE_MIN_DIS(1)

            if (dis >= top_dis) continue;

            q.push({internal_node->left, internal_node->is_left_leaf});
            q.push( {internal_node->right, internal_node->is_right_leaf});
        }
        else {  
            LeafNode* leaf_node = (LeafNode*) node_pair.first;
#if exp3_i_binary_use_bsax_to_prune
            float dis = min_dist_paa_to_bsax(search_paa, leaf_node->sax_lb, leaf_node->sax_ub);
#else
            float dis = min_dist_paa_to_isax(search_paa, leaf_node->sax_, leaf_node->card_);
#endif
            COUNT_NODE_COMPUTE_MIN_DIS(1)
            if (dis >= top_dis) continue;
            for (int i = 0; i < leaf_node->len; i ++ ) {
                leaf_ans.emplace_back(leaf_node->leaf_keys[i]);
            }
        }
    }
}

void iSAX2Tree::BFSM(ts_type* search_paa, float top_dis, vector<LeafNode*>& leaf_ans) const {
    queue<pair<void*, bool>> q; 
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
        if (!node_pair.second) {    
            InternalNode* internal_node = (InternalNode*) node_pair.first;
#if exp3_i_binary_use_bsax_to_prune
            float dis = min_dist_paa_to_bsax(search_paa, internal_node->sax_lb, internal_node->sax_ub);
#else
            float dis = min_dist_paa_to_isax(search_paa, internal_node->sax_, internal_node->card_);
#endif
            COUNT_NODE_COMPUTE_MIN_DIS(1)

            if (dis >= top_dis) continue;

            q.push({internal_node->left, internal_node->is_left_leaf});
            q.push( {internal_node->right, internal_node->is_right_leaf});
        }
        else {  
            LeafNode* leaf_node = (LeafNode*) node_pair.first;
#if exp3_i_binary_use_bsax_to_prune
            float dis = min_dist_paa_to_bsax(search_paa, leaf_node->sax_lb, leaf_node->sax_ub);
#else
            float dis = min_dist_paa_to_isax(search_paa, leaf_node->sax_, leaf_node->card_);
#endif
            COUNT_NODE_COMPUTE_MIN_DIS(1)
            if (dis >= top_dis) continue;
            leaf_ans.emplace_back(leaf_node);
        }
    }
}

void iSAX2Tree::ExactSearch(int k, TS* search_ts, vector<TS*>& ts_ans, vector<float>& top_dis) const {
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

void iSAX2Tree::ExactSearchM(int k, TS* search_ts, vector<TS*>& ts_ans, vector<float>& top_dis) const {
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


void iSAX2Tree::ExactSearchExp3(TS* search_ts, float bsf, vector<pair<SAX, u_int64_t>>& sax_and_p, vector<float>& top_dis) const {
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

}

void iSAX2Tree::materialized(const char *_filename) {
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

    queue<pair<void*, bool>> q; 
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
        if (!node_pair.second) {    
            InternalNode* internal_node = (InternalNode*) node_pair.first;
            pair<void*, bool> left_p = {internal_node->left, internal_node->is_left_leaf};
            q.emplace(left_p);
            pair<void*, bool> right_p = {internal_node->right, internal_node->is_right_leaf};
            q.emplace(right_p);
        }
        else {  
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
                    } else if (new_p == 5*tmp_ts_cap-1) {
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
            if (!node_pair.second) {    
                InternalNode* internal_node = (InternalNode*) node_pair.first;
                pair<void*, bool> left_p = {internal_node->left, internal_node->is_left_leaf};
                q.emplace(left_p);
                pair<void*, bool> right_p = {internal_node->right, internal_node->is_right_leaf};
                q.emplace(right_p);
            }
            else {  
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

void iSAX2Tree::materialized() {
    u_int64_t new_p = 0;

    queue<pair<void*, bool>> q; 
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
        if (!node_pair.second) {    
            InternalNode* internal_node = (InternalNode*) node_pair.first;
            pair<void*, bool> left_p = {internal_node->left, internal_node->is_left_leaf};
            q.emplace(left_p);
            pair<void*, bool> right_p = {internal_node->right, internal_node->is_right_leaf};
            q.emplace(right_p);
        }
        else {  
            LeafNode* leaf_node = (LeafNode*) node_pair.first;
            for (int i = 0; i < leaf_node->len; i ++ ) {
#if is_reorder_m
                if (new_p < 100000) {
                    leaf_node->leaf_keys[i].p = new_p + 4 * 100000;
                }
                else if (new_p < 500000){
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

