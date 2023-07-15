


#include <vector>
#include <algorithm>
#include "sax.h"

#if exp3_use_tree == 0
#include "b_binary/bSAX2Tree.h"
#elif exp3_use_tree == 1
#include "i_binary/iSAX2Tree.h"
#elif exp3_use_tree == 2
#include "btree/ESBTree.h"
#endif

#include "util/TimeRecord.h"
#include "CntRecord.h"

using namespace std;

int main() {


#if dataset_type == 0

    const char *filename = "../../dataset/ts.bin";
    const char *sax_filename = "../../dataset/ts_sax.bin";
    const char *query_filename = "../../dataset/ts_query.bin";



#elif dataset_type == 1
    const char *filename = "../../dataset/shift1e8.bin";
    const char *sax_filename = "../../dataset/shift1e8_sax.bin";
    const char *query_filename = "../../dataset/shift1b_query_var005.bin";


#elif dataset_type == 2
    const char *filename = "../../dataset/deep1b_99974688.bin";
    const char *sax_filename = "../../dataset/deep1b_sax_99974688.bin";
    const char *query_filename = "../../dataset/deep1b_query_var005.bin";
#elif dataset_type == 3
    const char *filename = "../../dataset/glove.100d.bin";
    const char *sax_filename = "../../dataset/glove_sax.100d.bin";
    const char *query_filename = "../../dataset/glove_query.100d.bin";
#endif

#if exp3_use_tree == 0
    bSAX2Tree tree(filename);
#elif exp3_use_tree == 1
    iSAX2Tree tree(filename);
#elif exp3_use_tree == 2
    ESBTree tree(filename);
#endif

    tree.BuildFromSAX(sax_filename);

    cout << "~~~~~~~~~~~~~~~~~~build finish~~~~~~~~~~~~~~~~~~~~" << endl;

    FILE * query_data_file = fopen(query_filename, "r");

    vector<TS> query_ts(NUM_SEARCH);
    for(int i = 0; i < NUM_SEARCH; i++) {
        fread(&query_ts[i], sizeof(TS) , 1, query_data_file);

    }
    fclose(query_data_file);

    vector<TS*> ts_ans;
    vector<float> top_dis;

    vector<float> bsf;
    float bsf_sum = 0;

#if exp3_write_index_ans
#if dataset_type == 0
    string bsf_answer_filename = "../../dataset_bsf_answers/ts";
#elif dataset_type == 1
    string bsf_answer_filename = "../../dataset_bsf_answers/shift";
#elif dataset_type == 2
        string bsf_answer_filename = "../../dataset_bsf_answers/deep";
#elif dataset_type == 3
        string bsf_answer_filename = "../../dataset_bsf_answers/glove";
#endif

    FILE * bsf_answer_file = fopen(bsf_answer_filename.c_str(), "w");

#endif

    for(int i = 0; i < NUM_SEARCH; i++) {
        tree.ApproximateSearch(K, &query_ts[i], ts_ans, top_dis);
        cout << "appro bsf dis:" << endl;
        cout << i << "  "<< top_dis[0] << endl;
        bsf.emplace_back(top_dis[0]);
        bsf_sum += pow(top_dis[0], 0.5);



#if exp3_write_index_ans
        for (int j = 0; j < K; j ++) {
            fwrite(&top_dis[j], 4, 1, bsf_answer_file);
            fwrite(&ts_ans[j], 8, 1, bsf_answer_file);
        }
#endif

        ts_ans.clear();
        top_dis.clear();
    }
    PRINT_COUNT
    PRINT_AVG_COUNT
    cout << "===============approximate search finish=================" << endl;

    exit(0);


    for(int i = 0; i < NUM_SEARCH; i++ ) {
        float top_k_dis = bsf[i];
        vector<pair<SAX, u_int64_t>> sax_and_p;
        tree.ExactSearchExp3(&query_ts[i], top_k_dis, sax_and_p, top_dis);

#if exp3_write_index_ans
#if dataset_type == 0
        string sax_answer_filename = "../../dataset_sax_answers/ts" + to_string(i);
#elif dataset_type == 1
        string sax_answer_filename = "../../dataset_sax_answers/shift" + to_string(i);
#elif dataset_type == 2
        string sax_answer_filename = "../../dataset_sax_answers/deep" + to_string(i);
#elif dataset_type == 3
        string sax_answer_filename = "../../dataset_sax_answers/glove" + to_string(i);
#endif
        FILE * sax_data_file = fopen(sax_answer_filename.c_str(), "w");
        if (!sax_data_file) {
            cout << "can't find" << sax_answer_filename << endl;
            exit(-1);
        }
        u_int64_t len = sax_and_p.size();
        fwrite(&len, sizeof(len), 1, sax_data_file);
        for (auto& sax_p: sax_and_p) {
            fwrite(sax_p.first.sax, sizeof(SAX), 1, sax_data_file);
            fwrite((void*) (&sax_p.second), sizeof(u_int64_t), 1, sax_data_file);
        }
        fclose (sax_data_file);
#endif
    }


    cout << "===============exact search finish=================" << endl;

    cout << "number of approximate:" << num_approximate_search_nodes * NUM_SEARCH << endl;
    cout << "avg dis:" << bsf_sum / bsf.size() << endl;

    PRINT_AVG_COUNT
    PRINT_TIME
    return 0;
}

