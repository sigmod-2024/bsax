#include <vector>
#include <algorithm>
#include "util/TopKHeap.h"
#include "sax.h"
#include "CntRecord.h"
#include "TimeRecord.h"

using namespace std;

struct LeafKey {
    SAX sax_;
    u_int64_t p;
};

//
//void GetTopKAns(const char* filename, int k, vector<LeafKey>& leaf_ans, TS* search_ts, ts_type* search_paa,
//                vector<TS*>& ts_ans, vector<float>& top_dis) {
//#if sort_strategy
//    vector<pair<float, long>> leaf_keys_pair(leaf_ans.size());  // min_dis, p
//#else
//    vector<pair<long, float>> leaf_keys_pair(leaf_ans.size()); // p, min_dis
//#endif
//    COMPUTE_MIN_DIS_START
//    for (int i = 0; i < leaf_ans.size(); i ++ ) {
//        LeafKey* leaf_key = &leaf_ans[i];
//        float dis = min_dist_paa_to_sax(search_paa, leaf_key->sax_);
//#if sort_strategy
//        leaf_keys_pair[i] = {dis, leaf_key->p};
//#else
//        leaf_keys_pair[i] = {leaf_key->p,  dis};
//#endif
//    }
//    COMPUTE_MIN_DIS_END
//    COUNT_EXACT_ANS(leaf_ans.size())
//
//    sort(leaf_keys_pair.begin(), leaf_keys_pair.end());
//
//    TopKHeap* heap = new TopKHeap(k);
//    TS* ts_arr = new TS[K+1];
//    TS* read_ts = ts_arr;
//    FILE * file;
//    file = fopen(filename, "r");
//    if (!file) {
//        cout << "can't find" << filename << endl;
//        exit(-1);
//    }
//#if sort_strategy == 2
//    int cnt = 0;
//    vector<pair<long, float>> leaf_keys_batch(sort_batch_size);
//#endif
//    for (int i = 0; i < leaf_keys_pair.size(); i ++ ) {
//        auto& leaf_key_pair = leaf_keys_pair[i];
//
//#if sort_strategy
//        long p = leaf_key_pair.second;
//        float dis = leaf_key_pair.first;
//#else
//        long p = leaf_key_pair.first;
//        float dis = leaf_key_pair.second;
//#endif
//
//        if (!heap->check_approximate(dis)) continue;
//#if sort_strategy <= 1
//        long offset = p * TS_LENGTH * sizeof(ts_type);
//        fseek(file, offset, SEEK_SET);
//        size_t bytes_read = fread(read_ts->ts, sizeof(TS), 1, file);
//        COUNT_EXACT_READ_TS(1)
//        float true_dis = ts_euclidean_distance(search_ts->ts, read_ts->ts);
//        heap->push_ans_approximate(true_dis, &read_ts);
//#else
//        leaf_keys_batch[cnt ++] = {p, dis};
//        if (cnt >= sort_batch_size || i == leaf_keys_pair.size() - 1) {
//            sort(leaf_keys_batch.begin(), leaf_keys_batch.begin() + cnt);
//
//            for (auto& one_leaf_key: leaf_keys_batch) {
//                long offset = one_leaf_key.first * TS_LENGTH * sizeof(ts_type);
//                fseek(file, offset, SEEK_SET);
//                size_t bytes_read = fread(read_ts->ts, sizeof(TS), 1, file);
//                COUNT_EXACT_READ_TS(1)
//                float true_dis = ts_euclidean_distance(search_ts->ts, read_ts->ts);
//                heap->push_ans_approximate(true_dis, &read_ts);
//            }
//            cnt = 0;
//        }
//#endif
//
//    }
//
//    fclose (file);
//
//    while(!heap->pq.empty()) {
//        top_dis.emplace_back(heap->pq.top().first);
//        ts_ans.emplace_back(heap->pq.top().second);
//        heap->pq.pop();
//    }
//    delete heap;
//}

void GetTopKAns(const char* filename, int k, vector<LeafKey>& leaf_ans, TS* search_ts, ts_type* search_paa,
                vector<TS*>& ts_ans, vector<float>& top_dis, float bsf) {
#if sort_strategy
    vector<pair<float, long>> leaf_keys_pair;  // min_dis, p
#else
    vector<pair<long, float>> leaf_keys_pair; // p, min_dis
#endif
    leaf_keys_pair.reserve(leaf_ans.size());
    COMPUTE_MIN_DIS_START
    for (int i = 0; i < leaf_ans.size(); i ++ ) {
        LeafKey* leaf_key = &leaf_ans[i];
        float dis = min_dist_paa_to_sax(search_paa, leaf_key->sax_);
        if (dis > bsf) continue;
#if sort_strategy
        leaf_keys_pair.emplace_back(dis, leaf_key->p);
#else
        leaf_keys_pair.emplace_back(leaf_key->p,  dis);
#endif
    }
    COMPUTE_MIN_DIS_END
    COUNT_EXACT_ANS(leaf_ans.size())

    sort(leaf_keys_pair.begin(), leaf_keys_pair.end());

    TopKHeap* heap = new TopKHeap(k);
    TS* ts_arr = new TS[K+1];
    TS* read_ts = ts_arr;
    FILE * file;
    file = fopen(filename, "r");
    if (!file) {
        cout << "can't find " << filename << endl;
        exit(-1);
    }
#if sort_strategy == 2
    int cnt = 0;
    int isbreak = false;
    vector<pair<long, float>> leaf_keys_batch(sort_batch_size);
#endif
    for (int i = 0; i < leaf_keys_pair.size(); i ++ ) {
        auto& leaf_key_pair = leaf_keys_pair[i];

#if sort_strategy
        long p = leaf_key_pair.second;
        float dis = leaf_key_pair.first;
#else
        long p = leaf_key_pair.first;
        float dis = leaf_key_pair.second;
#endif

#if sort_strategy == 1
        if (!heap->check_exact(dis, bsf)) break;
#elif sort_strategy == 2
    if (!heap->check_exact(dis, bsf)) {
        isbreak = true;
    }
    else leaf_keys_batch[cnt ++] = {p, dis};
#else
        if (!heap->check_exact(dis, bsf)) continue;
#endif


#if sort_strategy <= 1
        long offset = p * TS_LENGTH * sizeof(ts_type);
        fseek(file, offset, SEEK_SET);
        size_t bytes_read = fread(read_ts->ts, sizeof(TS), 1, file);
        COUNT_EXACT_READ_TS(1)
        float true_dis = ts_euclidean_distance(search_ts->ts, read_ts->ts);
        heap->push_ans_exact(true_dis, &read_ts, bsf);
#else

        if (cnt >= sort_batch_size || i == leaf_keys_pair.size() - 1 || isbreak) {

            sort(leaf_keys_batch.begin(), leaf_keys_batch.begin() + cnt);
            cout<<"once"<<endl;
            for (int j = 0; j < cnt; j ++) {
                if (!heap->check_exact(dis, bsf)) {
                    isbreak = true;
                    continue;
                }
                long offset = leaf_keys_batch[j].first * TS_LENGTH * sizeof(ts_type);
                fseek(file, offset, SEEK_SET);
                size_t bytes_read = fread(read_ts->ts, sizeof(TS), 1, file);

                COUNT_EXACT_READ_TS(1)
                float true_dis = ts_euclidean_distance(search_ts->ts, read_ts->ts);
                heap->push_ans_exact(true_dis, &read_ts, bsf);
            }

            cnt = 0;
        }
        if (isbreak) break;
#endif

    }

    fclose (file);



    while(!heap->pq.empty()) {
        top_dis.emplace_back(heap->pq.top().first);
        ts_ans.emplace_back(heap->pq.top().second);
        heap->pq.pop();
    }
    delete heap;
}

int main() {

#if dataset_type == 0
    const char *filename = "../../dataset/ts.bin";
    const char *query_filename = "../../dataset/ts_query.bin";
#elif dataset_type == 1
    const char *filename = "../../dataset/shift1e8.bin";
    const char *query_filename = "../../dataset/shift1b_query_var005.bin";

#elif dataset_type == 2
    const char *filename = "../../dataset/deep1b_99974688.bin";
    const char *query_filename = "../../dataset/deep1b_query_var005.bin";
#elif dataset_type == 3
    const char *filename = "../../dataset/glove.100d.bin";
    const char *query_filename = "../../dataset/glove_query.100d.bin";
#endif


    TS* search_ts = new TS();
    FILE * query_data_file = fopen(query_filename,"r");
    if (!query_data_file) {
        cout << "can't find " << query_filename << endl;
        exit(-1);
    }

    for (int i = 0; i < NUM_SEARCH; i ++ ) {
        fread(search_ts->ts, sizeof(TS) , 1, query_data_file);



#if dataset_type == 0
        string bsf_answer_filename = "../../dataset_bsf_answers/ts";
#elif dataset_type == 1
        string bsf_answer_filename = "../../dataset_bsf_answers/shift";
#elif dataset_type == 2
        string bsf_answer_filename = "../../dataset_bsf_answers/deep";
#elif dataset_type == 3
        string bsf_answer_filename = "../../dataset_bsf_answers/glove";
#endif
    vector<float> bsf(NUM_SEARCH);
    FILE * bsf_answer_file = fopen(bsf_answer_filename.c_str(), "r");

    for (int j = 0; j < NUM_SEARCH; j ++) {
        fread(&bsf[j], 4, 1, bsf_answer_file);
    }

    fclose(bsf_answer_file);


#if dataset_type == 0
        string sax_answer_filename = "../../dataset_sax_answers/ts" + to_string(i);
#elif dataset_type == 1
        string sax_answer_filename = "../../dataset_sax_answers/shift" + to_string(i);
#elif dataset_type == 2
        string sax_answer_filename = "../../dataset_sax_answers/deep" + to_string(i);
#elif dataset_type == 3
        string sax_answer_filename = "../../dataset_sax_answers/glove" + to_string(i);
#endif
        FILE * sax_data_file = fopen(sax_answer_filename.c_str(), "r");
        if (!sax_data_file) {
            cout << "can't find:" << sax_answer_filename << endl;
            exit(-1);
        }
        u_int64_t len;
        fread(&len, sizeof(len), 1, sax_data_file);
        vector<LeafKey> leaf_ans(len);
        for (int j = 0; j < len; j ++ ) {
            fread(&leaf_ans[j].sax_, sizeof(SAX), 1, sax_data_file);
            fread(&leaf_ans[j].p, sizeof(u_int64_t), 1, sax_data_file);
        }
        fclose (sax_data_file);


        ts_type search_paa[SEGMENTS];
        paa_from_ts(search_ts->ts, search_paa);
        vector<TS*> ts_ans;
        vector<float> top_dis;
        cout << "BSF ans: " << bsf[i] << endl;
        EXACT_GET_ANS_START
        GetTopKAns(filename, K, leaf_ans, search_ts, search_paa, ts_ans, top_dis, bsf[i]);
        EXACT_GET_ANS_END


        cout << i << endl;

        PRINT_COUNT
        COUNT_CLEAR
    }
    cout << "get ans total time:" << total_exact_get_ans_time / NUM_SEARCH << endl << " compute dis time:" << total_compute_min_dis / NUM_SEARCH << endl;
    cout << "access time" << total_exact_get_ans_time / NUM_SEARCH - total_compute_min_dis / NUM_SEARCH <<endl;
    cout << "number of access raw data:" << (double) total_exact_count_read_ts / NUM_SEARCH / TOTAL_TS << endl;
    fclose(query_data_file);
}