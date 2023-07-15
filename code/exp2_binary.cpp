#include <iostream>
#include <vector>
#include <queue>
#include <algorithm>
#include "globals.h"
#include "sax.h"

#if exp2_binary_use_isax
#include "i_binary/iSAX2Tree.h"
#else
#include "b_binary/bSAX2Tree.h"
#endif
using namespace std;

void compute_dis(ts_type* search_paa, LeafNode* leaf_node, vector<float>& bsax_ratios, vector<float>& isax_ratios) {
    float min_key_dis = MAXFLOAT;

    CARD key_card;  key_card.set_max_card();

    SAX node_sax_lb, node_sax_ub;
    node_sax_lb.set_max_value();
    node_sax_ub.set_min_value();

    SAX node_isax;
    CARD node_card;

    for (int j = 0; j < leaf_node->len; j ++ ) {
        SAX key_sax = leaf_node->leaf_keys[j].sax_;

        float key_dis = min_dist_paa_to_isax(search_paa, key_sax, key_card);
        min_key_dis = min(min_key_dis, key_dis);


        for (int k = 0; k < SEGMENTS; k ++ ) {
            node_sax_lb.sax[k] = min(key_sax.sax[k], node_sax_lb.sax[k]);
            node_sax_ub.sax[k] = max(key_sax.sax[k], node_sax_ub.sax[k]);
        }


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
    float node_bsax_dis = min_dist_paa_to_bsax(search_paa, node_sax_lb, node_sax_ub);
    float node_isax_dis = min_dist_paa_to_isax(search_paa, node_isax, node_card);


    if(min_key_dis == 0) {
        bsax_ratios.emplace_back(1);
        isax_ratios.emplace_back(1);
    } else {
        float bsax_rations_2 = node_bsax_dis / min_key_dis;
        float isax_ratios_2 = node_isax_dis / min_key_dis;
        if (bsax_rations_2 > 1e-6) bsax_ratios.emplace_back(pow(bsax_rations_2, 0.5));
        else bsax_ratios.emplace_back(0);
        if (isax_ratios_2 > 1e-6) isax_ratios.emplace_back(pow(isax_ratios_2, 0.5));
        else isax_ratios.emplace_back(0);
    }

}

void bfs(ts_type* search_paa, void* _tree, vector<float>& bsax_ratios, vector<float>& isax_ratios) {
#if exp2_binary_use_isax
    iSAX2Tree* tree = (iSAX2Tree*)_tree;
#else
    bSAX2Tree* tree = (bSAX2Tree*)_tree;
#endif
    queue<pair<void*, bool>> q;
#if binary_tree_root_full
    for (auto& node_pair: tree->root->children) {
        q.emplace(node_pair);
    }
#else
    q.push({tree->root->node, tree->root->root_is_leaf});
#endif
    while(!q.empty()) {
        pair<void*, bool> node_pair = q.front();
        q.pop();
        if (!node_pair.second) {
            InternalNode* internal_node = (InternalNode*) node_pair.first;
            q.push({internal_node->left, internal_node->is_left_leaf});
            q.push( {internal_node->right, internal_node->is_right_leaf});
        }
        else {
            LeafNode* leaf_node = (LeafNode*) node_pair.first;
            compute_dis(search_paa, leaf_node, bsax_ratios, isax_ratios);
        }
    }
}




int main() {
//    const char *filename = "../../dataset/ts.bin";
//    const char *sax_filename = "../../dataset/ts_sax.bin";
//    const char *query_filename = "../../dataset/ts_query.bin";
//    const char *bsax_answer_filename = "../../dataset_answers/ts_ibinary_bsax.bin";
//    const char *isax_answer_filename = "../../dataset_answers/ts_ibinary_isax.bin";

//    const char *filename = "../../dataset/glove.100d.bin";
//    const char *sax_filename = "../../dataset/glove_sax.100d.bin";
//    const char *query_filename = "../../dataset/glove_query.100d.bin";
//    const char *bsax_answer_filename = "../../dataset_answers/glove_bbinary_bsax.bin";
//    const char *isax_answer_filename = "../../dataset_answers/glove_bbinary_isax.bin";

    const char *filename = "../../dataset/deep1b_99974688.bin";
    const char *sax_filename = "../../dataset/deep1b_sax_99974688.bin";
    const char *query_filename = "../../dataset/deep1b_query_var005.bin";
    const char *bsax_answer_filename = "../../dataset_answers/deep1b_ibinary_bsax.bin";
    const char *isax_answer_filename = "../../dataset_answers/deep1b_ibinary_isax.bin";

//    const char *filename = "../../dataset/shift1e8.bin";
//    const char *sax_filename = "../../dataset/shift1e8_sax.bin";
//    const char *query_filename = "../../dataset/shift1b_query_var005.bin";
//    const char *bsax_answer_filename = "../../dataset_answers/shift1e8_ibinary_bsax.bin";
//    const char *isax_answer_filename = "../../dataset_answers/shift1e8_ibinary_isax.bin";

#if exp2_binary_use_isax
    iSAX2Tree tree(filename);
    tree.BuildFromSAX(sax_filename);
#else
    bSAX2Tree tree(filename);
    tree.BuildFromSAX(sax_filename);
#endif


    cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;

    TS* search_ts = new TS();
    FILE* data_file = fopen(query_filename,"r");
    FILE* bsax_answer_file = fopen(bsax_answer_filename,"w");
    FILE* isax_answer_file = fopen(isax_answer_filename,"w");
    if (!data_file) {
        cout << "can't find" << filename << endl;
        exit(-1);
    }

    for(int i = 0; i<100; i++) {
        fread(search_ts->ts, sizeof(TS) , 1, data_file);
        ts_type search_paa[SEGMENTS];
        paa_from_ts(search_ts->ts, search_paa);

        vector<float> bsax_ratios, isax_ratios;

        bfs(search_paa, (void *)&tree, bsax_ratios, isax_ratios);
        cout<<i<<" "<<bsax_ratios.size()<<endl;
        for(float & bsax_ratio : bsax_ratios) {
            fwrite(&bsax_ratio, 4, 1, bsax_answer_file);
        }
        for(float & isax_ratio : isax_ratios) {
            fwrite(&isax_ratio, 4, 1, isax_answer_file);
        }
    }



    fclose(data_file);
    fclose(bsax_answer_file);
    fclose(isax_answer_file);
    delete search_ts;

    return 0;
}
