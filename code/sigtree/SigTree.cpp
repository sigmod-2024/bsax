


#include <vector>
#include <algorithm>
#include "SigTree.h"
#include "RootNode.h"
#include "../util/TopKHeap.h"
#include "CntRecord.h"
#include "TimeRecord.h"
#include "unordered_set"
using namespace std;



void ConvertLeafKey(RootNode *root) {
    queue<pair<void *, bool>> q; 
    for (auto &node_pair: root->children) {
#if sigtree_delay_create_leafnode
        q.emplace(node_pair.second);
#else
        q.emplace(node_pair);
#endif
    }
    while (!q.empty()) {
        pair<void *, bool> node_pair = q.front();
        q.pop();
        if (!node_pair.second) {    
            InternalNode *internal_node = (InternalNode *) node_pair.first;
            for (auto& child_pair: internal_node->children) {
#if sigtree_delay_create_leafnode
                q.emplace(child_pair.second);
#else
                q.emplace(child_pair);
#endif
            }
        } else {  
            LeafNode *leaf_node = (LeafNode *) node_pair.first;

#if sigtree_leaf_keys_use_vector
            for (int i = 0; i < leaf_node->leaf_keys.size(); i++) {
#else
                for (int i = 0; i < leaf_node->len; i ++ ) {
#endif
                LeafKey *leaf_key = &leaf_node->leaf_keys[i];

                SAX sax_;
                sax_from_saxt((saxt_type *) &leaf_key->sax_.sax, sax_.sax);
                leaf_key->sax_ = sax_;  
            }
        }
    }
}

void read_ts_batch(int n, FILE *data_file, vector<SAXT> &saxt_vec, vector<TS> &ts_vec, vector<SAX> &sax_vec) {
    BUILD_READ_TS_START
    for (int i = 0; i < n; i++) {
        fread(&ts_vec[i], sizeof(TS), 1, data_file);
    }
    BUILD_READ_TS_END
    BUILD_CONVERT_SAX_START
    for (int i = 0; i < n; i++) {
        sax_from_ts(ts_vec[i].ts, sax_vec[i].sax);
    }
    BUILD_CONVERT_SAX_END
    for (int i = 0; i < n; i++) {
        saxt_from_sax(sax_vec[i].sax, saxt_vec[i].saxt);
    }
}

void SigTree::Build() {
    this->BuildTree();

    FILE *data_file = fopen(filename, "r");
    if (!data_file) {
        cout << "" << filename << endl;
        exit(-1);
    }
    vector<SAXT> saxt_vec(READ_TS_BATCH);
    vector<TS> ts_vec(READ_TS_BATCH);
    vector<SAX> sax_vec(READ_TS_BATCH);
    BUILD_TIME_START
    for (int i = 0; i < (TOTAL_TS - 1) / READ_TS_BATCH + 1; i++) {
        int read_batch = READ_TS_BATCH;
        if (TOTAL_TS % READ_TS_BATCH && i == (TOTAL_TS - 1) / READ_TS_BATCH) {
            read_batch = TOTAL_TS % READ_TS_BATCH;
        }
        read_ts_batch(read_batch, data_file, saxt_vec, ts_vec, sax_vec);
        for (int j = 0; j < read_batch; j++) {
            this->Insert(&saxt_vec[j], i * READ_TS_BATCH + j);
        }
    }
    ConvertLeafKey(root);
    BUILD_TIME_END

    fclose(data_file);
}

void SigTree::BuildFromSAX(const char *sax_filename) {
    this->BuildTree();
    FILE *data_file = fopen(sax_filename, "r");
    if (!data_file) {
        cout << "" << sax_filename << endl;
        exit(-1);
    }
    vector<SAXT> saxt_vec(TOTAL_TS);
    vector<SAX> sax_vec(TOTAL_TS);
    for (int i = 0; i < TOTAL_TS; i++) {
        fread(&sax_vec[i], sizeof(SAX), 1, data_file);
    }

    BUILD_INDEX_START
    for (int i = 0; i < TOTAL_TS; i++) {
        saxt_from_sax(sax_vec[i].sax, saxt_vec[i].saxt);
        this->Insert(&saxt_vec[i], i);
        if ((i + 1) % 1000000 == 0) {
            cout << i + 1 << endl;
        }
    }
    ConvertLeafKey(root);
    BUILD_INDEX_END

    fclose(data_file);
}


void SigTree::BuildTree() {
    root = new RootNode();
    for (u_int64_t i = 0; i < (1 << SEGMENTS); i++) {
        SAXT saxt_;
        std::memset(saxt_.saxt, 0, sizeof(saxt_));
        saxt_.saxt[BIT_CARDINALITY - 1] = i;
#if sigtree_delay_create_leafnode
#else
        root->children[i] = {new LeafNode(saxt_, 1), true};
#endif
    }
}

LeafNode *SplitNode(LeafNode *leaf_node, void *&pre_node, bool &pre_is_root, const SAXT *insert_saxt) {

    if (leaf_node->card_ >= BIT_CARDINALITY) {
        printf("\n");
        exit(-1);
    }


    InternalNode *new_internal_node = new InternalNode(leaf_node->saxt_, leaf_node->card_);
    
    u_int8_t new_card = leaf_node->card_ + 1;
    SAXT new_saxt = leaf_node->saxt_;
    for (int i = 0; i < (1 << SEGMENTS); i++) {
        new_saxt.saxt[BIT_CARDINALITY - new_card] = i;
#if sigtree_delay_create_leafnode
#else
        LeafNode *new_leaf_node = new LeafNode(new_saxt, new_card);
        new_internal_node->children[i] = {new_leaf_node, true};
#endif
    }


#if sigtree_leaf_keys_use_vector
    for (int i = 0; i < leaf_node->leaf_keys.size(); i++) {
#else
        for (int i = 0; i < leaf_node->len; i ++ ) {
#endif
        SAXT *leaf_keys_saxt = (SAXT *) &leaf_node->leaf_keys[i].sax_;   


#if sigtree_delay_create_leafnode
        if (new_internal_node->children.count(leaf_keys_saxt->saxt[BIT_CARDINALITY - new_card]) == 0) {
            new_saxt.saxt[BIT_CARDINALITY - new_card] = leaf_keys_saxt->saxt[BIT_CARDINALITY -  new_card];
            new_internal_node->children[leaf_keys_saxt->saxt[BIT_CARDINALITY - new_card]] = {new LeafNode(new_saxt, new_card), true};
        }
#endif
        LeafNode *new_leaf_node = (LeafNode *) new_internal_node->children[leaf_keys_saxt->saxt[BIT_CARDINALITY - new_card]].first;


#if sigtree_leaf_keys_use_vector
        new_leaf_node->leaf_keys.emplace_back(leaf_node->leaf_keys[i]);
#else
        new_leaf_node->leaf_keys[new_leaf_node->len] = leaf_node->leaf_keys[i];
        new_leaf_node->len ++;
#endif
    }


    
    if (pre_is_root) {  
        RootNode *root_node = (RootNode *) pre_node;
        saxt_type sax_first_bit = insert_saxt->saxt[BIT_CARDINALITY - 1];
        root_node->children[sax_first_bit] = {new_internal_node, false};
    } else {
        InternalNode *internal_node = (InternalNode *) pre_node;
        saxt_type sax_card_now_bit = leaf_node->saxt_.saxt[BIT_CARDINALITY - leaf_node->card_];
        internal_node->children[sax_card_now_bit] = {new_internal_node, false};
    };
    delete leaf_node;

    pre_node = new_internal_node;
    pre_is_root = false;





#if sigtree_delay_create_leafnode
    if (new_internal_node->children.count(insert_saxt->saxt[BIT_CARDINALITY - new_card]) == 0) {
        new_saxt.saxt[BIT_CARDINALITY - new_card] = insert_saxt->saxt[BIT_CARDINALITY -  new_card];
        new_internal_node->children[insert_saxt->saxt[BIT_CARDINALITY - new_card]] = {new LeafNode(new_saxt, new_card), true};
    }
#endif

    return (LeafNode *) new_internal_node->children[insert_saxt->saxt[BIT_CARDINALITY - new_card]].first;
}

void InsertNode(const SAXT *insert_saxt, const u_int64_t p, void *pre_node, bool pre_is_root,
                void *node, bool is_leaf, u_int8_t card_now) {

    if (!is_leaf) {
        InternalNode *internal_node = (InternalNode *) node;

        card_now ++;
        
#if sigtree_delay_create_leafnode
        if (internal_node->children.count(insert_saxt->saxt[BIT_CARDINALITY - card_now]) == 0) {
            SAXT saxt_now = internal_node->saxt_;
            saxt_now.saxt[BIT_CARDINALITY - card_now] = insert_saxt->saxt[BIT_CARDINALITY - card_now];
            internal_node->children[insert_saxt->saxt[BIT_CARDINALITY - card_now]] = {new LeafNode(saxt_now, card_now), true};
        }
#endif

        std::pair<void *, bool> child_pair = internal_node->children[insert_saxt->saxt[BIT_CARDINALITY - card_now]];

        InsertNode(insert_saxt, p, internal_node, false, child_pair.first, child_pair.second, card_now);
    } else {
        LeafNode *leaf_node = (LeafNode *) node;

#if sigtree_leaf_keys_use_vector
        while (leaf_node->leaf_keys.size() >= LEAF_MAX_NUM) {
#else
        if (leaf_node->len >= LEAF_MAX_NUM) {
#endif
            leaf_node = SplitNode(leaf_node, pre_node, pre_is_root, insert_saxt);
        }


#if sigtree_leaf_keys_use_vector
        LeafKey leaf_key;
        leaf_key.sax_ = *(SAX *) &insert_saxt->saxt;   
        leaf_key.p = p;
        leaf_node->leaf_keys.emplace_back(leaf_key);

#else
        LeafKey* leaf_key = &leaf_node->leaf_keys[leaf_node->len];
        leaf_node->len ++ ;
        leaf_key->sax_ = *(SAX*)&insert_saxt->saxt;   
        leaf_key->p = p;
#endif
    }


}


void SigTree::Insert(const SAXT *insert_saxt, u_int64_t p) {

#if sigtree_delay_create_leafnode
    if (root->children.count(insert_saxt->saxt[BIT_CARDINALITY - 1]) == 0) {
        SAXT saxt_;
        std::memset(saxt_.saxt, 0, sizeof(saxt_));
        saxt_.saxt[BIT_CARDINALITY - 1] = insert_saxt->saxt[BIT_CARDINALITY - 1];
        root->children[insert_saxt->saxt[BIT_CARDINALITY - 1]] = {new LeafNode(saxt_, 1), true};
    }
#endif

    std::pair<void *, bool> node_pair = root->children[insert_saxt->saxt[BIT_CARDINALITY - 1]];
    InsertNode(insert_saxt, p, root, true, node_pair.first, node_pair.second, 1);
}


void GetTopKAns(const char *filename, int k, vector<LeafKey> &leaf_ans, TS *search_ts, ts_type *search_paa,
                vector<TS *> &ts_ans, vector<float> &top_dis) {
    vector<DisP> dis_p(leaf_ans.size());

    for (int i = 0; i < leaf_ans.size(); i++) {
        float low_dis = min_dist_paa_to_sax(search_paa, leaf_ans[i].sax_);
        dis_p[i] = {low_dis, leaf_ans[i].p};
    }

    FILE *file;
    file = fopen(filename, "r");
    if (!file) {
        cout << "" << filename << endl;
        exit(-1);
    }
    TopKHeap *heap = new TopKHeap(k);
    TS *ts_arr = new TS[K + 1];
    TS *read_ts = ts_arr;


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
    for (int i = 0; i < dis_p.size(); i++) {

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

    while (!heap->pq.empty()) {
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
        COUNT_EXACT_ANS(leaf_an->leaf_keys.size())

        for(int i=0;i<leaf_an->leaf_keys.size();i++) {
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


void
SigTree::DFS(const SAXT *search_saxt, void *node, bool is_leaf, vector<LeafKey> &leaf_ans, u_int8_t card_now) const {
    if (leaf_ans.size() >= num_approximate_search_key) return;
    if (!is_leaf) {
        InternalNode *internal_node = (InternalNode *) node;
#if sigtree_delay_create_leafnode
        if (internal_node->children.count(search_saxt->saxt[BIT_CARDINALITY - card_now - 1])) {
            std::pair<void *, bool> child_pair = internal_node->children[search_saxt->saxt[BIT_CARDINALITY - card_now - 1]];
            DFS(search_saxt, child_pair.first, child_pair.second, leaf_ans, card_now + 1);
        }









        
        if (leaf_ans.size() < num_approximate_search_key) {
            auto l = internal_node->children.lower_bound(search_saxt->saxt[BIT_CARDINALITY - card_now - 1]);
            auto r = l;
            l --;
            if (r != root->children.end()) r ++;
            while(l != internal_node->children.end() || r != internal_node->children.end() && leaf_ans.size() < num_approximate_search_key) {
                if (r != internal_node->children.end()) {
                    std::pair<void *, bool> child_pair = r->second;
                    DFS(search_saxt, child_pair.first, child_pair.second, leaf_ans, card_now + 1);
                    r ++;
                }
                if (l != internal_node->children.end()) {
                    std::pair<void *, bool> child_pair = l->second;
                    DFS(search_saxt, child_pair.first, child_pair.second, leaf_ans, card_now + 1);
                    if (l == internal_node->children.begin()) {
                        l = internal_node->children.end();
                    }
                    else {
                        l --;
                    }
                }
            }
        }
#else
        std::pair<void *, bool> child_pair = internal_node->children[search_saxt->saxt[BIT_CARDINALITY - card_now - 1]];
        DFS(search_saxt, child_pair.first, child_pair.second, leaf_ans, card_now + 1);
        
        if (leaf_ans.size() < num_approximate_search_key) {
            long int l = search_saxt->saxt[BIT_CARDINALITY - card_now - 1] - 1, r =
                    search_saxt->saxt[BIT_CARDINALITY - card_now - 1] + 1;
            while (l >= 0 || r < ((u_int64_t)1 << SEGMENTS) && leaf_ans.size() < num_approximate_search_key) {
                if (r < ((u_int64_t)1 << SEGMENTS)) {
                    child_pair = internal_node->children[r];
                    DFS(search_saxt, child_pair.first, child_pair.second, leaf_ans, card_now + 1);
                    r++;
                }
                if (l >= 0) {
                    child_pair = internal_node->children[l];
                    DFS(search_saxt, child_pair.first, child_pair.second, leaf_ans, card_now + 1);
                    l--;
                }
            }
        }
#endif
    } else {
        LeafNode *leaf_node = (LeafNode *) node;
#if sigtree_leaf_keys_use_vector
        for (int i = 0; i < leaf_node->leaf_keys.size(); i++) {
#else
            for (int i = 0; i < leaf_node->len; i ++ ) {
#endif
            leaf_ans.emplace_back(leaf_node->leaf_keys[i]);
            if (leaf_ans.size() >= num_approximate_search_key) {
                break;
            }
        }
    }
}

void SigTree::ApproximateSearch(int k, TS *search_ts, vector<TS *> &ts_ans, vector<float> &top_dis) const {
    APPROXIMATE_TIME_START
    SAXT search_saxt;
    saxt_from_ts(search_ts->ts, search_saxt.saxt);
    vector<LeafKey> leaf_ans;

    APPROXIMATE_INDEX_START


#if sigtree_delay_create_leafnode
    std::pair<void *, bool> node_pair = root->children.lower_bound(search_saxt.saxt[BIT_CARDINALITY - 1])->second;
    DFS(&search_saxt, node_pair.first, node_pair.second, leaf_ans, 1);
    










    if (leaf_ans.size() < num_approximate_search_key) {
        auto l = root->children.lower_bound(search_saxt.saxt[BIT_CARDINALITY - 1]);
        auto r = l;
        l --;
        if (r != root->children.end()) r ++;
        while(l != root->children.end() || r != root->children.end() && leaf_ans.size() < num_approximate_search_key) {
            if (r != root->children.end()) {
                std::pair<void *, bool> child_pair = r->second;
                DFS(&search_saxt, child_pair.first, child_pair.second, leaf_ans, 1);
                r ++;
            }
            if (l != root->children.end()) {
                std::pair<void *, bool> child_pair = l->second;
                DFS(&search_saxt, child_pair.first, child_pair.second, leaf_ans, 1);
                if (l == root->children.begin()) {
                    l = root->children.end();
                }
                else {
                    l --;
                }
            }
        }
    }
#else
    std::pair<void *, bool> node_pair = root->children[search_saxt.saxt[BIT_CARDINALITY - 1]];
    DFS(&search_saxt, node_pair.first, node_pair.second, leaf_ans, 1);
    
    if (leaf_ans.size() < num_approximate_search_key) {
        long int l = search_saxt.saxt[BIT_CARDINALITY - 1] - 1, r = search_saxt.saxt[BIT_CARDINALITY - 1] + 1;
        while (l >= 0 || r < ((u_int64_t)1 << SEGMENTS) && leaf_ans.size() < num_approximate_search_key) {
            if (r < ((u_int64_t)1 << SEGMENTS)) {
                node_pair = root->children[r];
                DFS(&search_saxt, node_pair.first, node_pair.second, leaf_ans, 1);
                r++;
            }
            if (l >= 0) {
                node_pair = root->children[l];
                DFS(&search_saxt, node_pair.first, node_pair.second, leaf_ans, 1);
                l--;
            }
        }
    }
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



void SigTree::BFS(ts_type *search_paa, float top_dis, vector<LeafKey> &leaf_ans) const {
    queue<pair<void *, bool>> q; 

    cout<<"root_size:"<<root->children.size()<<endl;

    for (auto &node_pair: root->children) {
#if sigtree_delay_create_leafnode
        q.emplace(node_pair.second);
#else
        q.emplace(node_pair);
#endif
    }

    while (!q.empty()) {
        pair<void *, bool> node_pair = q.front();
        q.pop();
        if (!node_pair.second) {    
            InternalNode *internal_node = (InternalNode *) node_pair.first;
            float dis = min_dist_paa_to_saxt(search_paa, internal_node->saxt_, internal_node->card_);
            COUNT_NODE_COMPUTE_MIN_DIS(1)

            if (dis >= top_dis) continue;

            for (auto &child_pair: internal_node->children) {
#if sigtree_delay_create_leafnode
                q.emplace(child_pair.second);
#else
                q.emplace(child_pair);
#endif
            }
        } else {  
            LeafNode *leaf_node = (LeafNode *) node_pair.first;

            
            if (leaf_node->leaf_keys.size() <= 20) {
                for(auto & leaf_key : leaf_node->leaf_keys) {
                    leaf_ans.emplace_back(leaf_key);
                }
                continue;
            }

            float dis = min_dist_paa_to_saxt(search_paa, leaf_node->saxt_, leaf_node->card_);
            COUNT_NODE_COMPUTE_MIN_DIS(1)
            if (dis >= top_dis) continue;

#if sigtree_leaf_keys_use_vector
            for (int i = 0; i < leaf_node->leaf_keys.size(); i++) {
#else
                for (int i = 0; i < leaf_node->len; i ++ ) {
#endif
                leaf_ans.emplace_back(leaf_node->leaf_keys[i]);
            }
        }
    }
}

void SigTree::BFSM(ts_type *search_paa, float top_dis, vector<LeafNode*> &leaf_ans) const {
    queue<pair<void *, bool>> q; 

    cout<<"root_size:"<<root->children.size()<<endl;

    for (auto &node_pair: root->children) {
#if sigtree_delay_create_leafnode
        q.emplace(node_pair.second);
#else
        q.emplace(node_pair);
#endif
    }

    while (!q.empty()) {
        pair<void *, bool> node_pair = q.front();
        q.pop();
        if (!node_pair.second) {    
            InternalNode *internal_node = (InternalNode *) node_pair.first;
            float dis = min_dist_paa_to_saxt(search_paa, internal_node->saxt_, internal_node->card_);
            COUNT_NODE_COMPUTE_MIN_DIS(1)

            if (dis >= top_dis) continue;

            for (auto &child_pair: internal_node->children) {
#if sigtree_delay_create_leafnode
                q.emplace(child_pair.second);
#else
                q.emplace(child_pair);
#endif
            }
        } else {  
            LeafNode *leaf_node = (LeafNode *) node_pair.first;

            
            if (leaf_node->leaf_keys.size() <= 20) {
                leaf_ans.emplace_back(leaf_node);
                continue;
            }

            float dis = min_dist_paa_to_saxt(search_paa, leaf_node->saxt_, leaf_node->card_);
            COUNT_NODE_COMPUTE_MIN_DIS(1)
            if (dis >= top_dis) continue;
            leaf_ans.emplace_back(leaf_node);


        }
    }
}

void SigTree::ExactSearch(int k, TS *search_ts, vector<TS *> &ts_ans, vector<float> &top_dis) const {
    EXACT_TIME_START
    vector<TS *> approximate_ts_ans;
    vector<float> approximate_top_dis;
    ApproximateSearch(k, search_ts, approximate_ts_ans, approximate_top_dis);

    for (float dis : approximate_top_dis) {
        cout << dis << " ";
    }
    cout << endl;

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

void SigTree::ExactSearchM(int k, TS *search_ts, vector<TS *> &ts_ans, vector<float> &top_dis) const {
    EXACT_TIME_START
    vector<TS *> approximate_ts_ans;
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

void SigTree::materialized(const char *_filename) {
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

    queue<pair<void *, bool>> q; 


    for (auto &node_pair: root->children) {
#if sigtree_delay_create_leafnode
        q.emplace(node_pair.second);
#else
        q.emplace(node_pair);
#endif
    }

    while (!q.empty()) {
        pair<void *, bool> node_pair = q.front();
        q.pop();
        if (!node_pair.second) {    
            InternalNode *internal_node = (InternalNode *) node_pair.first;

            for (auto &child_pair: internal_node->children) {
#if sigtree_delay_create_leafnode
                q.emplace(child_pair.second);
#else
                q.emplace(child_pair);
#endif
            }
        } else {  
            LeafNode *leaf_node = (LeafNode *) node_pair.first;



#if sigtree_leaf_keys_use_vector
            for (int i = 0; i < leaf_node->leaf_keys.size(); i++) {
#else
                for (int i = 0; i < leaf_node->len; i ++ ) {
#endif
                if (leaf_node->leaf_keys[i].p < read_batch_size) {
                    tmp_ts[tmp_ts_size] = tmp_read_ts[leaf_node->leaf_keys[i].p];
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
        for (auto &node_pair: root->children) {
#if sigtree_delay_create_leafnode
            q.emplace(node_pair.second);
#else
            q.emplace(node_pair);
#endif
        }

        while (!q.empty()) {
            pair<void *, bool> node_pair = q.front();
            q.pop();
            if (!node_pair.second) {    
                InternalNode *internal_node = (InternalNode *) node_pair.first;

                for (auto &child_pair: internal_node->children) {
#if sigtree_delay_create_leafnode
                    q.emplace(child_pair.second);
#else
                    q.emplace(child_pair);
#endif
                }
            } else {  
                LeafNode *leaf_node = (LeafNode *) node_pair.first;



#if sigtree_leaf_keys_use_vector
                for (int i = 0; i < leaf_node->leaf_keys.size(); i++) {
#else
                    for (int i = 0; i < leaf_node->len; i ++ ) {
#endif
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

void SigTree::materialized() {
    u_int64_t new_p = 0;

    queue<pair<void *, bool>> q; 


    for (auto &node_pair: root->children) {
#if sigtree_delay_create_leafnode
        q.emplace(node_pair.second);
#else
        q.emplace(node_pair);
#endif
    }

    while (!q.empty()) {
        pair<void *, bool> node_pair = q.front();
        q.pop();
        if (!node_pair.second) {    
            InternalNode *internal_node = (InternalNode *) node_pair.first;

            for (auto &child_pair: internal_node->children) {
#if sigtree_delay_create_leafnode
                q.emplace(child_pair.second);
#else
                q.emplace(child_pair);
#endif
            }
        } else {  
            LeafNode *leaf_node = (LeafNode *) node_pair.first;



#if sigtree_leaf_keys_use_vector
            for (int i = 0; i < leaf_node->leaf_keys.size(); i++) {
#else
                for (int i = 0; i < leaf_node->len; i ++ ) {
#endif
                leaf_node->leaf_keys[i].p = new_p++;
            }
        }
    }
}

















