#include <iostream>
#include <vector>
#include <queue>
#include <algorithm>
#include "globals.h"
#include "sax.h"
#include "util/TopKHeap.h"

using  namespace  std;
int main() {
    const char *filename = "../../dataset/ts.bin";
    const char *sax_filename = "../../dataset/ts_sax.bin";
    const char *filename1 = "../../dataset/ts_query_new.bin";

//    FILE* data_file = fopen(filename,"r");
    FILE* data_file_new = fopen(filename1,"r");


    vector<TS> ts(NUM_SEARCH);
    TS t;
    SAX s;
    ts_type search_paa[SEGMENTS];
    SAX news;

    float bsf = 16.4898;

    for(int i=0;i<NUM_SEARCH;i++) {
        fread(&ts[i], sizeof(TS), 1, data_file_new);

        paa_from_ts(ts[i].ts, search_paa);

        FILE* data_file = fopen(filename,"r");
        FILE* sax_data_file = fopen(sax_filename,"r");
        TopKHeap* heap = new TopKHeap(K);
        for(int j=0;j<TOTAL_TS;j++) {
            fread(&t, sizeof(TS), 1, data_file);
            fread(&s, sizeof(SAX), 1, sax_data_file);
//            sax_from_ts(t.ts, news.sax);

            float low_dis = min_dist_paa_to_sax(search_paa, s);
//            float low_dis1 = min_dist_paa_to_sax(search_paa, news);

//            if(low_dis != low_dis1) exit(2);

            float true_dis = ts_euclidean_distance(ts[i].ts, t.ts, TS_LENGTH);
//            if(low_dis > 16.4898) continue;
//            if(low_dis > true_dis) exit(1);
//            if(!heap->check_approximate(low_dis)) continue;
            if(bsf < low_dis) continue;
//            heap->push_ans_approximate(true_dis, j);
            heap->push_ans_exact(true_dis, j, bsf);
        }

        cout<<"ppp:"<<endl;
        vector<u_int64_t> ppp;
        for(int k=0;k<K;k++) {
            ppp.push_back((u_int64_t)heap->pq.top().second);
            heap->pq.pop();
        }

        sort(ppp.begin(), ppp.end());
        for(auto item: ppp) {
            cout<<item<<endl;
        }
    }


}
