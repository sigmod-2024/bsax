



#include <vector>
#include <algorithm>
#include "util/TopKHeap.h"
#include "sax.h"
#include "CntRecord.h"
#include "TimeRecord.h"
#include <chrono>
#include <unordered_set>
#include "set"
using namespace std;

struct LeafKey {
    SAX sax_;
    u_int64_t p;
};

bool cmpp(LeafKey& a, LeafKey& b) {
    return  a.p < b.p;
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
        cout << "can't find:" << filename << endl;
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
        COUNT_EXACT_READ_TS(1)
        float true_dis = ts_euclidean_distance(search_ts->ts, read_ts->ts);
        heap->push_ans_approximate(true_dis, &read_ts);
    }
#endif


#if sort_strategy == 1
    sort(dis_p.begin(), dis_p.end(), DisCmp);
    for(int i = 0; i < dis_p.size(); i ++) {

        if (!heap->check_approximate(dis_p[i].dis)) break;


        long offset = dis_p[i].p * TS_LENGTH * sizeof(ts_type);
        fseek(file, offset, SEEK_SET);
        size_t bytes_read = fread(read_ts->ts, sizeof(TS), 1, file);
        COUNT_EXACT_READ_TS(1)
        float true_dis = ts_euclidean_distance(search_ts->ts, read_ts->ts);
        heap->push_ans_approximate(true_dis, &read_ts);
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
            COUNT_EXACT_READ_TS(1)
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
    std::chrono::steady_clock::time_point t1;
    std::chrono::steady_clock::time_point t2;
    t1 = std::chrono::steady_clock::now();
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


    delete heap;
//    delete[] ts_arr;
    delete read_ts;
    t2 = std::chrono:: steady_clock::now();
    std::cout <<"total"<< std::chrono::duration_cast<std::chrono::milliseconds>( t2-t1 ).count() <<"ms==========================================="<< std::endl;
}





int main() {

#if dataset_type == 0
    const char *filename = "../../dataset/ts.bin";
    const char *query_filename = "../../dataset/ts_query.bin";
//    const char *query_filename = "../../dataset/ts_query_new.bin";
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
        cout << "can't find:" << query_filename << endl;
        exit(-1);
    }


#if dataset_type == 0
    string bsf_answer_filename = "../../dataset_bsf_answers/ts";
#elif dataset_type == 1
    string bsf_answer_filename = "../../dataset_bsf_answers/shift";
#elif dataset_type == 2
        string bsf_answer_filename = "../../dataset_bsf_answers/deep";
#elif dataset_type == 3
        string bsf_answer_filename = "../../dataset_bsf_answers/glove";
#endif

    FILE * bsf_answer_file = fopen(bsf_answer_filename.c_str(), "r");




//    for (int i = 0; i < NUM_SEARCH - 31; i ++ ) {
//        fread(search_ts->ts, sizeof(TS), 1, query_data_file);
//    }

//    for (int i = NUM_SEARCH - 31; i < NUM_SEARCH; i ++ ) {
    for (int i = 0; i < NUM_SEARCH; i ++ ) {
        fread(search_ts->ts, sizeof(TS) , 1, query_data_file);

        vector<float> bsf(K);
        vector<TS*> bsf_p(K);
        for (int j = 0; j < K; j ++) {
            fread(&bsf[j], 4, 1, bsf_answer_file);
            fread(&bsf_p[j], 8, 1, bsf_answer_file);
//            cout<<bsf[j]<<" "<<(u_int64_t)bsf_p[j]<<endl;
        }



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
            cout << "" << sax_answer_filename << endl;
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

    sort(leaf_ans.begin(), leaf_ans.end(), cmpp);
//    for(int i=0;i<leaf_ans.size();i++) {
//        if(i%10000 == 0)
//        cout<<(u_int64_t)leaf_ans[i].p<<endl;
//    }
        ts_type search_paa[SEGMENTS];
        paa_from_ts(search_ts->ts, search_paa);
        vector<TS*> ts_ans;
        vector<float> top_dis;

        EXACT_GET_ANS_START
        GetTopKAns(filename, K, leaf_ans, search_ts, search_paa, ts_ans, top_dis, bsf, bsf_p);
        EXACT_GET_ANS_END


        cout << i << endl;
        cout << ":" << total_exact_get_ans_time / (i + 1) << endl << " :" << total_compute_min_dis / (i + 1) << endl;

        PRINT_COUNT
        COUNT_CLEAR
    }
    cout << ":" << total_exact_get_ans_time / NUM_SEARCH << endl << " :" << total_compute_min_dis / NUM_SEARCH << endl;
    cout << "" << total_exact_get_ans_time / NUM_SEARCH - total_compute_min_dis / NUM_SEARCH <<endl;
    cout << ":" << (double) total_exact_count_read_ts / NUM_SEARCH / TOTAL_TS << endl;
    fclose(query_data_file);
}



