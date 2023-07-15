



#include "ESBTree.h"
#include "LeafNode.h"
#include "sax.h"
#include "../util/TopKHeap.h"
#include "CntRecord.h"
#include "TimeRecord.h"
#include <vector>
#include <queue>
#include <algorithm>
#include "unordered_set"
using namespace std;

void ESBTree:: read_ts_batch(int n, FILE* data_file, vector<LeafKey>& leaf_keys, long offset, vector<TS>& ts_vec, vector<SAX>& sax_vec) {
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

    for (int i = 0; i < n; i ++ ) {
        saxt_from_sax(sax_vec[i].sax, (saxt_type *) &leaf_keys[offset + i].sax_);   
        leaf_keys[offset + i].p = offset + i;
    }


}

void ESBTree::Build() {
    FILE * data_file = fopen(filename,"r");
    if (!data_file) {
        cout << "" << filename << endl;
        exit(-1);
    }
    vector<LeafKey> leaf_keys(TOTAL_TS);
    vector<TS> ts_vec(READ_TS_BATCH);
    vector<SAX> sax_vec(READ_TS_BATCH);

    BUILD_TIME_START
    for (int i = 0; i < (TOTAL_TS - 1) / READ_TS_BATCH + 1; i ++ ) {
        int read_batch = READ_TS_BATCH;
        if (TOTAL_TS % READ_TS_BATCH && i == (TOTAL_TS - 1) / READ_TS_BATCH) {
            read_batch = TOTAL_TS % READ_TS_BATCH;
        }
        read_ts_batch(read_batch, data_file, leaf_keys, i * READ_TS_BATCH, ts_vec, sax_vec);
    }
    fclose (data_file);


    this->BuildTree(leaf_keys);
    BUILD_TIME_END
}

void ESBTree::BuildFromSAX(const char *sax_filename) {
    FILE * data_file = fopen(sax_filename,"r");
    if (!data_file) {
        cout << "" << sax_filename << endl;
        exit(-1);
    }
    vector<LeafKey> leaf_keys(TOTAL_TS);
    vector<SAX> sax_vec(TOTAL_TS);
    for(int i = 0; i < TOTAL_TS; i++) {
        fread(&sax_vec[i], sizeof(SAX), 1, data_file);
    }

    BUILD_INDEX_START
    for(int i = 0; i < TOTAL_TS; i++ ) {
        saxt_from_sax(sax_vec[i].sax, (saxt_type *) &leaf_keys[i].sax_);

        leaf_keys[i].p = i;
    }
    this->BuildTree(leaf_keys);
    BUILD_INDEX_END

    fclose (data_file);
}



void BuildInternalNode(LeafNode* input_leaf_nodes, u_int64_t num_input_nodes,
                       InternalNode* output_internal_nodes) {
    int cnt_internal_node = 0;
#if btree_use_bsax
    SAX sax_lb;
    SAX sax_ub;
#else
    SAX new_sax;
    CARD new_card;
#endif
    SAXT saxt_lb;
    for (int i = 0; i < num_input_nodes; i ++ ) {
        LeafNode* leaf_node = &input_leaf_nodes[i];
#if btree_use_bsax
        sax_lb.set_max_value();
        sax_ub.set_min_value();
        for (int j = 0; j < leaf_node->len; j ++ ) {
            SAX sax_now;
            sax_from_saxt((saxt_type*) &leaf_node->leaf_keys[j].sax_,    
                          sax_now.sax);
            for (int k = 0; k < SEGMENTS; k ++ ) {
                sax_lb.sax[k] = min(sax_now.sax[k], sax_lb.sax[k]);
                sax_ub.sax[k] = max(sax_now.sax[k], sax_ub.sax[k]);
            }
            if (j == 0) {
                saxt_lb = *(SAXT*)& leaf_node->leaf_keys[j].sax_;
            }
            leaf_node->leaf_keys[j].sax_ = sax_now;
        }
#else
        for (int j = 0; j < leaf_node->len; j ++ ) {
            SAX sax_now;
            sax_from_saxt((saxt_type *) &leaf_node->leaf_keys[j].sax_,    
                               sax_now.sax);
            if (j == 0) {
                new_sax = sax_now;
                new_card.set_max_card();
                saxt_lb = *(SAXT*)& leaf_node->leaf_keys[j].sax_;
            }
            else {
                for (int k = 0; k < SEGMENTS; k ++ ) {
                    for (int l = 1; l <= new_card.card[k]; l ++ ) {
                        saxt_type bit = 1 << (BIT_CARDINALITY - l);
                        if ((new_sax.sax[k] & bit) != (sax_now.sax[k] & bit)) {
                            new_card.card[k] = l - 1;
                            new_sax.sax[k] &= bit;
                        }
                    }
                }
            }
            leaf_node->leaf_keys[j].sax_ = sax_now;
        }
#endif
        InternalNode* internal_node = &output_internal_nodes[cnt_internal_node];
        InternalKey* internal_key = &(internal_node->internal_keys[internal_node->len]);
        internal_key->p = leaf_node;
        internal_key->is_leaf = true;

        internal_key->saxt_lb = saxt_lb;
#if btree_use_bsax
        internal_key->sax_lb = sax_lb;
        internal_key->sax_ub = sax_ub;
#else
        internal_key->sax_ = new_sax;
        internal_key->card_ = new_card;
#endif
        internal_node->len ++;
        if (internal_node->len >= LEAF_MAX_NUM) {
            cnt_internal_node ++;
        }
    }
}


void BuildInternalNode(InternalNode* input_internal_nodes, u_int64_t num_input_nodes,
                       InternalNode* output_internal_nodes) {
    int cnt_internal_node = 0;
#if btree_use_bsax
    SAX sax_lb;
    SAX sax_ub;
#else
    SAX new_sax;
    CARD new_card;
#endif
    for (int i = 0; i < num_input_nodes; i ++ ) {
        InternalNode* input_internal_node = &input_internal_nodes[i];

#if btree_use_bsax
        sax_lb.set_max_value();
        sax_ub.set_min_value();
        for (int j = 0; j < input_internal_node->len; j ++ ) {
            InternalKey* internal_key = &input_internal_node->internal_keys[j];
            for (int k = 0; k < SEGMENTS; k ++ ) {
                sax_lb.sax[k] = min(internal_key->sax_lb.sax[k], sax_lb.sax[k]);
                sax_ub.sax[k] = max(internal_key->sax_ub.sax[k], sax_ub.sax[k]);
            }
        }
#else
        for (int j = 0; j < input_internal_node->len; j ++ ) {
            InternalKey* internal_key = &input_internal_node->internal_keys[j];

            if (j == 0) {
                new_sax = internal_key->sax_;
                new_card.set_max_card();
            }
            else {
                for (int k = 0; k < SEGMENTS; k ++ ) {
                    new_card.card[k] = min(new_card.card[k], internal_key->card_.card[k]);
                    for (int l = 1; l <= new_card.card[k]; l ++ ) {
                        saxt_type bit = 1 << (BIT_CARDINALITY - l);
                        if ((new_sax.sax[k] & bit) != (internal_key->sax_.sax[k] & bit)) {
                            new_card.card[k] = l - 1;
                            new_sax.sax[k] &= bit;
                        }
                    }
                }
            }
        }
#endif

        InternalNode* internal_node = &output_internal_nodes[cnt_internal_node];
        InternalKey* internal_key = &(internal_node->internal_keys[internal_node->len]);
        internal_key->p = input_internal_node;
        internal_key->is_leaf = false;

        internal_key->saxt_lb = input_internal_node->internal_keys[0].saxt_lb;

#if btree_use_bsax
        internal_key->sax_lb = sax_lb;
        internal_key->sax_ub = sax_ub;
#else
        internal_key->sax_ = new_sax;
        internal_key->card_ = new_card;
#endif
        internal_node->len ++;
        if (internal_node->len >= LEAF_MAX_NUM) {
            cnt_internal_node ++;
        }
    }
}

bool comp_saxt(LeafKey& x, LeafKey& y) {    
    return *(SAXT*)&x.sax_ < *(SAXT*)&y.sax_;
}

void ESBTree::BuildTree(vector<LeafKey>& leaf_keys) {
    u_int64_t n = leaf_keys.size();

    sort(leaf_keys.begin(), leaf_keys.end(), comp_saxt);    

    num_leafs = (n - 1) / LEAF_MAX_NUM + 1;
    leaf_nodes = new LeafNode[num_leafs];
    for (int i = 0; i < n / LEAF_MAX_NUM; i ++ ) {
        leaf_nodes[i].add_batch(&leaf_keys[i * LEAF_MAX_NUM], LEAF_MAX_NUM);
        leaf_nodes[i].id = i;
    }
    leaf_nodes[num_leafs - 1].add_batch(&leaf_keys[n / LEAF_MAX_NUM * LEAF_MAX_NUM], n - n / LEAF_MAX_NUM * LEAF_MAX_NUM);


    u_int64_t num_nodes = (num_leafs - 1) / LEAF_MAX_NUM + 1;
    InternalNode* internal_nodes = new InternalNode[num_nodes];
    BuildInternalNode(leaf_nodes, num_leafs, internal_nodes);

    while(num_nodes > 1) {
        u_int64_t num_output_nodes = (num_nodes - 1) / LEAF_MAX_NUM + 1;
        InternalNode* new_internal_nodes = new InternalNode[num_output_nodes];
        BuildInternalNode(internal_nodes, num_nodes, new_internal_nodes);
        internal_nodes = new_internal_nodes;
        num_nodes = num_output_nodes;
    }

    root = internal_nodes;
}







































void ESBTree::DFS(const SAXT* search_saxt, void* node, bool is_leaf, vector<LeafNode*>& leaf_ans) const {
    if (!is_leaf) {

        InternalNode* internal_node = (InternalNode*) node;
        int l = 0, r = internal_node->len - 1;
        while(l < r) {
            int mid = (l + r + 1) >> 1;
            if (internal_node->internal_keys[mid].saxt_lb <= *search_saxt) {
                l = mid;
            }
            else {
                r = mid - 1;
            }
        }

        if (!internal_node->internal_keys[l].is_leaf) {
            DFS(search_saxt, internal_node->internal_keys[l].p, false, leaf_ans);
        }
        else {
            DFS(search_saxt, internal_node->internal_keys[l].p, true, leaf_ans);
        }
    }
    else {

        LeafNode* leaf_node = (LeafNode*) node;
        int l = max(0, (int)leaf_node->id - (num_approximate_search_nodes - 1) / 2);
        int cnt = 0;
        for (int i = l; i < num_leafs && cnt < num_approximate_search_nodes; i ++ ) {
            leaf_ans.emplace_back(&leaf_nodes[i]);
            cnt ++;
        }
        if (cnt < num_approximate_search_nodes) {   
            for (int i = l - 1; i >= 0 && cnt < num_approximate_search_nodes; i -- ) {
                leaf_ans.emplace_back(&leaf_nodes[i]);
                cnt ++;
            }
        }
    }
}


void ESBTree::BFS(ts_type* search_paa, float top_dis, vector<LeafNode*>& leaf_ans) const {
    queue<void*> q;
    bool flag = true, next_level_is_leaf = false;

    q.emplace(root);
    while(!q.empty()) {
        int curr_level_size = q.size();
        for (int i = 0; i < curr_level_size; i ++ ) {
            if (flag) {
                InternalNode* node = (InternalNode *)q.front();
                for (int j = 0; j < node->len; j ++ ) {
                    InternalKey internal_key = node->internal_keys[j];
                    if (internal_key.is_leaf) next_level_is_leaf = true;
#if btree_use_bsax
                    float dis = min_dist_paa_to_bsax(search_paa, internal_key.sax_lb, internal_key.sax_ub);
#else
                    float dis = min_dist_paa_to_isax(search_paa, internal_key.sax_, internal_key.card_);
#endif
                    COUNT_NODE_COMPUTE_MIN_DIS(1)

                    if (dis >= top_dis) continue;

                    q.emplace((InternalNode *) internal_key.p);
                }
            }
            else {
                LeafNode* node = (LeafNode *) q.front();
                leaf_ans.emplace_back(node);
            }
            q.pop();
        }
        if (next_level_is_leaf) flag = false;
    }
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

void ESBTree::ApproximateSearch(int k, TS* search_ts, vector<TS*>& ts_ans, vector<float>& top_dis) const {
    APPROXIMATE_TIME_START
    SAXT search_saxt;
    saxt_from_ts(search_ts->ts, search_saxt.saxt);

    vector<LeafNode*> leaf_ans;
    ts_type search_paa[SEGMENTS];
    paa_from_ts(search_ts->ts, search_paa);

    APPROXIMATE_INDEX_START
    DFS(&search_saxt, root, false, leaf_ans);

    vector<LeafKey> leaf_ans_real;
    for (int i = 0; i < leaf_ans.size(); i ++) {
        for (int j = 0; j < leaf_ans[i]->len; j ++) {
            leaf_ans_real.emplace_back(leaf_ans[i]->leaf_keys[j]);
        }
    }
    APPROXIMATE_INDEX_END


    COUNT_APPROXIMATE_ANS(leaf_ans_real.size())

    APPROXIMATE_GET_ANS_START
    GetTopKAns(filename, k, leaf_ans_real, search_ts, search_paa, ts_ans, top_dis);
    APPROXIMATE_GET_ANS_END

    APPROXIMATE_TIME_END
}

void ESBTree::ExactSearch(int k, TS* search_ts, vector<TS*>& ts_ans, vector<float>& top_dis) const {
    EXACT_TIME_START
    vector<TS*> approximate_ts_ans;
    vector<float> approximate_top_dis;
    ApproximateSearch(k, search_ts, approximate_ts_ans, approximate_top_dis);
    float bsf = approximate_top_dis[0];

    for (float dis : approximate_top_dis) {
        cout << dis << " ";
    }
    cout << endl;

    ts_type search_paa[SEGMENTS];
    paa_from_ts(search_ts->ts, search_paa);

    vector<LeafNode*> leaf_ans;

    EXACT_INDEX_START
    BFS(search_paa, bsf, leaf_ans);
    vector<LeafKey> leaf_ans_real;
    for (int i = 0; i < leaf_ans.size(); i ++) {
        for (int j = 0; j < leaf_ans[i]->len; j ++) {
            leaf_ans_real.emplace_back(leaf_ans[i]->leaf_keys[j]);
        }
    }
    EXACT_INDEX_END


    COUNT_EXACT_ANS(leaf_ans_real.size())

    EXACT_GET_ANS_START
    GetTopKAns(filename, k, leaf_ans_real, search_ts, search_paa, ts_ans, top_dis, approximate_top_dis, approximate_ts_ans);
    EXACT_GET_ANS_END

    EXACT_TIME_END

    cout << "bsf:" << bsf << endl;
    PRINT_COUNT
    COUNT_CLEAR
}

void ESBTree::ExactSearchM(int k, TS* search_ts, vector<TS*>& ts_ans, vector<float>& top_dis) const {
    EXACT_TIME_START
    vector<TS*> approximate_ts_ans;
    vector<float> approximate_top_dis;
    ApproximateSearch(k, search_ts, approximate_ts_ans, approximate_top_dis);
    float bsf = approximate_top_dis[0];

    for (float dis : approximate_top_dis) {
        cout << dis << " ";
    }
    cout << endl;

    ts_type search_paa[SEGMENTS];
    paa_from_ts(search_ts->ts, search_paa);

    vector<LeafNode*> leaf_ans;

    EXACT_INDEX_START
    BFS(search_paa, bsf, leaf_ans);


    EXACT_GET_ANS_START
    GetTopKAnsM(filename, k, leaf_ans, search_ts, search_paa, ts_ans, top_dis, approximate_top_dis, approximate_ts_ans);
    EXACT_GET_ANS_END

    EXACT_TIME_END

    cout << "bsf:" << bsf << endl;
    PRINT_COUNT
    COUNT_CLEAR
}



void ESBTree::ExactSearchExp3(TS* search_ts, float top_k_dis, vector<pair<SAX, u_int64_t>>& sax_and_p, vector<float>& top_dis) const {
    ts_type search_paa[SEGMENTS];
    paa_from_ts(search_ts->ts, search_paa);
    vector<LeafNode*> leaf_ans;

    EXACT_INDEX_START
    BFS(search_paa, top_k_dis, leaf_ans);
    EXACT_INDEX_END

    for (auto leaf_node : leaf_ans) {
        for (auto& leaf_key: leaf_node->leaf_keys) {
            sax_and_p.emplace_back(leaf_key.sax_, leaf_key.p);
        }
    }
    COUNT_EXACT_ANS(sax_and_p.size())

    cout << "topk:" << top_k_dis << endl;
    PRINT_COUNT
    COUNT_CLEAR
}

void ESBTree::materialized(const char *_filename) {
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
    TS* tmp_ts = new TS[1000000];
    u_int64_t tmp_ts_size = 0;
    u_int64_t new_p = 0;
    u_int64_t read_batch_size = (u_int64_t) build_batch_size * 1000 *1000 * 1024 / 4 / TS_LENGTH;
    TS* tmp_read_ts = new TS[read_batch_size];

    for(int j=0;j< read_batch_size && j < TOTAL_TS;j++) {
        fread(tmp_read_ts + j, sizeof(TS), 1, data_file);
    }

    queue<void*> q;
    bool flag = true, next_level_is_leaf = false;

    q.emplace(root);
    while(!q.empty()) {
        int curr_level_size = q.size();
        for (int i = 0; i < curr_level_size; i ++ ) {
            if (flag) {
                InternalNode* node = (InternalNode *)q.front();
                for (int j = 0; j < node->len; j ++ ) {
                    InternalKey internal_key = node->internal_keys[j];
                    if (internal_key.is_leaf) next_level_is_leaf = true;

                    q.emplace((InternalNode *) internal_key.p);
                }
            }
            else {
                LeafNode* node = (LeafNode *) q.front();
                for (int j = 0; j < node->len; j ++ ) {
                    if (node->leaf_keys[j].p < read_batch_size) {
                        tmp_ts[tmp_ts_size] = tmp_read_ts[node->leaf_keys[j].p];
                    }
                    tmp_ts_size++;
                    if(tmp_ts_size == 1000000) {
                        cout<<new_p<<endl;
                        fwrite(tmp_ts, sizeof(TS), 1000000, m_data_file);
                        tmp_ts_size = 0;
                    }
                    new_p++;
                }
            }
            q.pop();
        }
        if (next_level_is_leaf) flag = false;
    }
    if(tmp_ts_size != 0) {
        fwrite(tmp_ts, sizeof(TS), tmp_ts_size, m_data_file);
    }
    delete[] tmp_ts;

    for (u_int64_t kk = read_batch_size; kk < TOTAL_TS; kk += read_batch_size) {

        for(int j=0;j< read_batch_size && kk + j < TOTAL_TS;j++) {
            fread(tmp_read_ts + j, sizeof(TS), 1, data_file);
        }
        cout<<"new_batch"<<endl;
        new_p = 0;
        flag = true, next_level_is_leaf = false;

        q.emplace(root);
        while(!q.empty()) {
            int curr_level_size = q.size();
            for (int i = 0; i < curr_level_size; i ++ ) {
                if (flag) {
                    InternalNode* node = (InternalNode *)q.front();
                    for (int j = 0; j < node->len; j ++ ) {
                        InternalKey internal_key = node->internal_keys[j];
                        if (internal_key.is_leaf) next_level_is_leaf = true;

                        q.emplace((InternalNode *) internal_key.p);
                    }
                }
                else {
                    LeafNode* node = (LeafNode *) q.front();
                    for (int j = 0; j < node->len; j ++ ) {
                        if (node->leaf_keys[j].p < kk + read_batch_size && node->leaf_keys[j].p >= kk) {

                            fseek(m_data_file, new_p * TS_LENGTH * sizeof(ts_type), SEEK_SET);
                            fwrite(&tmp_read_ts[node->leaf_keys[j].p - kk], sizeof(TS), 1, m_data_file);
                        }
                        new_p++;
                    }
                }
                q.pop();
            }
            if (next_level_is_leaf) flag = false;
        }
    }


}

void ESBTree::materialized() {
    u_int64_t new_p = 0;

    queue<void*> q;
    bool flag = true, next_level_is_leaf = false;

    q.emplace(root);
    while(!q.empty()) {
        int curr_level_size = q.size();
        for (int i = 0; i < curr_level_size; i ++ ) {
            if (flag) {
                InternalNode* node = (InternalNode *)q.front();
                for (int j = 0; j < node->len; j ++ ) {
                    InternalKey internal_key = node->internal_keys[j];
                    if (internal_key.is_leaf) next_level_is_leaf = true;

                    q.emplace((InternalNode *) internal_key.p);
                }
            }
            else {
                LeafNode* node = (LeafNode *) q.front();
                for (int j = 0; j < node->len; j ++ ) {
                    node->leaf_keys[j].p = new_p++;
                }

            }
            q.pop();
        }
        if (next_level_is_leaf) flag = false;
    }
}


