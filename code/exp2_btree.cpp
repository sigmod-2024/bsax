#include <iostream>
#include <vector>
#include <algorithm>
#include "globals.h"
#include "sax.h"
#include "btree/ESBTree.h"
#include "btree/LeafNode.h"

using namespace std;


int main() {
//    const char *filename = "../../dataset/ts.bin";
//    const char *sax_filename = "../../dataset/ts_sax.bin";
//    const char *query_filename = "../../dataset/ts_query.bin";
//    const char *bsax_answer_filename = "../../dataset_answers/ts_btree_bsax.bin";
//    const char *isax_answer_filename = "../../dataset_answers/ts_btree_isax.bin";

//    const char *filename = "../../dataset/glove.100d.bin";
//    const char *sax_filename = "../../dataset/glove_sax.100d.bin";
//    const char *query_filename = "../../dataset/glove_query.100d.bin";
//    const char *bsax_answer_filename = "../../dataset_answers/glove_btree_bsax.bin";
//    const char *isax_answer_filename = "../../dataset_answers/glove_btree_isax.bin";

    const char *filename = "../../dataset/deep1b_99974688.bin";
    const char *sax_filename = "../../dataset/deep1b_sax_99974688.bin";
    const char *query_filename = "../../dataset/deep1b_query_var005.bin";
    const char *bsax_answer_filename = "../../dataset_answers/deep1b_btree_bsax.bin";
    const char *isax_answer_filename = "../../dataset_answers/deep1b_btree_isax.bin";


//    const char *filename = "../../dataset/shift1e8.bin";
//    const char *sax_filename = "../../dataset/shift1e8_sax.bin";
//    const char *query_filename = "../../dataset/shift1b_query_var005.bin";
//    const char *bsax_answer_filename = "../../dataset_answers/shift1e8_btree_bsax.bin";
//    const char *isax_answer_filename = "../../dataset_answers/shift1e8_btree_isax.bin";

    ESBTree esb_tree(filename);
    esb_tree.BuildFromSAX(sax_filename);

    cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;

    TS* search_ts = new TS();
    FILE * data_file = fopen(query_filename,"r");
    FILE* bsax_answer_file = fopen(bsax_answer_filename,"w");
    FILE* isax_answer_file = fopen(isax_answer_filename,"w");
    if (!data_file) {
        cout << "can't find" << query_filename << endl;
        exit(-1);
    }

    for(int ii=0;ii<100;ii++) {

        fread(search_ts->ts, sizeof(TS) , 1, data_file);
        ts_type search_paa[SEGMENTS];
        paa_from_ts(search_ts->ts, search_paa);


        LeafNode* leaf_nodes = esb_tree.leaf_nodes;
        int num_leafs = esb_tree.num_leafs;
        float isax_ratios[num_leafs], bsax_ratios[num_leafs];

        for (int i = 0; i < num_leafs; i ++ ) {
            float min_key_dis = MAXFLOAT;

            CARD key_card;  key_card.set_max_card();

            SAX node_sax_lb, node_sax_ub;
            node_sax_lb.set_max_value();
            node_sax_ub.set_min_value();

            SAX node_isax;
            CARD node_card;

            for (int j = 0; j < leaf_nodes[i].len; j ++ ) {
                SAX key_sax = leaf_nodes[i].leaf_keys[j].sax_;

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
                bsax_ratios[i] = 1;
                isax_ratios[i] = 1;
            } else {
                float bsax_rations_2 = node_bsax_dis / min_key_dis;
                float isax_ratios_2 = node_isax_dis / min_key_dis;
                if (bsax_rations_2 > 1e-6) bsax_ratios[i] = pow(bsax_rations_2, 0.5);
                else bsax_ratios[i] = 0;
                if (isax_ratios_2 > 1e-6) isax_ratios[i] = pow(isax_ratios_2, 0.5);
                else isax_ratios[i] = 0;
            }
        }

        cout<<ii<<" "<<num_leafs<<endl;
        for(float & bsax_ratio : bsax_ratios) {
            fwrite(&bsax_ratio, 4, 1, bsax_answer_file);
        }
        for(float & isax_ratio : isax_ratios) {
            fwrite(&isax_ratio, 4, 1, isax_answer_file);
        }
    }

    fclose(data_file);
    return 0;
}
