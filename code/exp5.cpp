#include <vector>
#include <algorithm>
#include "sax.h"

#if exp5_use_tree == 0
#include "b_binary/bSAX2Tree.h"
#elif exp5_use_tree == 1
#include "i_binary/iSAX2Tree.h"
#elif exp5_use_tree == 2
#include "btree/ESBTree.h"
#elif exp5_use_tree == 3
#include "sigtree/SigTree.h"
#endif

#include "util/TimeRecord.h"
#include "CntRecord.h"

using namespace std;

int main() {


#if dataset_type == 0
    const char *filename = "../../dataset/ts.bin";
    const char *sax_filename = "../../dataset/ts_sax.bin";
    const char *query_filename = "../../dataset/ts_query.bin";
//    const char *query_filename = "../../dataset/ts.bin";
    const char *m_filename = "../../dataset_m/ts.bin";

#elif dataset_type == 1
    const char *filename = "../../dataset/shift1e8.bin";
    const char *sax_filename = "../../dataset/shift1e8_sax.bin";
    const char *query_filename = "../../dataset/shift1b_query_var005.bin";
//    const char *query_filename = "../../dataset/shift1e8.bin";
    const char *m_filename = "../../dataset_m/shift1e8.bin";

#elif dataset_type == 2
    const char *filename = "../../dataset/deep1b_99974688.bin";
    const char *sax_filename = "../../dataset/deep1b_sax_99974688.bin";
    const char *query_filename = "../../dataset/deep1b_query_var005.bin";
    const char *m_filename = "../../dataset_m/deep1b_99974688.bin";
#elif dataset_type == 3
    const char *filename = "../../dataset/glove.100d.bin";
    const char *sax_filename = "../../dataset/glove_sax.100d.bin";
    const char *query_filename = "../../dataset/glove_query.100d.bin";
    const char *m_filename = "../../dataset_m/glove.100d.bin";
#endif

#if exp5_use_tree == 0
#if is_query_m
    bSAX2Tree tree(m_filename);
#else
    bSAX2Tree tree(filename);
#endif
#elif exp5_use_tree == 1
#if is_query_m
    iSAX2Tree tree(m_filename);
#else
    iSAX2Tree tree(filename);
#endif
#elif exp5_use_tree == 2
#if is_query_m
    ESBTree tree(m_filename);
#else
    ESBTree tree(filename);
#endif
#elif exp5_use_tree == 3
#if is_query_m
    SigTree tree(m_filename);
#else
    SigTree tree(filename);
#endif
#endif

//    tree.Build();
    tree.BuildFromSAX(sax_filename);

#if is_build_m
    BUILD_M_START
    tree.materialized(m_filename);
    BUILD_M_END
#endif

    PRINT_TIME
    cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;
//    exit(0);

#if is_query_m
    tree.materialized();
#endif


    FILE * query_data_file = fopen(query_filename, "r");

    vector<TS> query_ts(NUM_SEARCH);
    for(int i = 0; i < NUM_SEARCH; i++) {
        fread(&query_ts[i], sizeof(TS) , 1, query_data_file);
    }
    fclose(query_data_file);

    vector<TS*> ts_ans;
    vector<float> top_dis;


//    for(int i = 0; i < NUM_SEARCH; i++ ) {
//        tree.ApproximateSearch(K, &query_ts[i], ts_ans, top_dis);
//        for (float dis : top_dis) {
//            cout << dis << " ";
//        }
//        cout << endl;
//        cout << "--------------------------------------------" << endl;
//        ts_ans.clear();
//        top_dis.clear();
//    }

    for(int i = 0; i < NUM_SEARCH; i++ ) {
        cout << ":" << i + 1 << endl;
#if is_use_m
        tree.ExactSearchM(K, &query_ts[i], ts_ans, top_dis);
#else
        tree.ExactSearch(K, &query_ts[i], ts_ans, top_dis);
#endif

        for (float dis : top_dis) {
            cout << dis << " ";
        }
        cout << endl;
        cout << "--------------------------------------------" << endl;
        ts_ans.clear();
        top_dis.clear();
    }
    PRINT_AVG_COUNT
    PRINT_TIME
    cout << "================================" << endl;
    return 0;
}

