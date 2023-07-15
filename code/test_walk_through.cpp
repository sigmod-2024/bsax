#include <iostream>
#include <vector>
#include <algorithm>
#include <queue>
#include "globals.h"
#include "sax.h"
#include "../util/TopKHeap.h"

using namespace std;
int main() {

    const char *filename = "../ts.bin";
    FILE * data_file;
    data_file = fopen(filename,"r");

    int n = 1000000;
    vector<TS> ts_vec(n);
    for(int i = 0; i < n; i ++) {
        fread(&ts_vec[i], sizeof(TS) , 1, data_file);

//        SAX sax1;
//        SAX sax2;
//        SAXT saxt;
//        saxt_from_ts(ts_vec[i].ts, saxt.saxt);
//
//        sax_from_saxt(saxt.saxt, sax1.sax);
//        sax_from_saxt_simd(saxt.saxt, sax2.sax);
//        for (int j = 0; j < SEGMENTS; j ++ ) {
////            if (sax1.sax[j] != sax2.sax[j]) {
//                saxt_print_bit(saxt.saxt);
//                sax_print_bit(sax1);
//                sax_print_bit(sax2);
//                cout << i << endl;
//                exit(0);
////            }
//        }
    }
    cout << " finish " << endl;
    fclose (data_file);

    int k = 100;
    TopKHeap* heap = new TopKHeap(k);
    TS* search_ts = &ts_vec[0];
    for (int i = 0; i < n; i ++ ) {
        TS* ts = new TS;

        (*ts) = ts_vec[i];
        float dis = ts_euclidean_distance(ts->ts, search_ts->ts);
        if (!heap->check_approximate(dis)) continue;
        heap->push_ans_approximate(dis, ts);
    }

    while(!heap->pq.empty()) {
        cout << heap->pq.top().first << endl;

//        TS* ts = heap->pq.top().second;
//        SAX* sax = new SAX();
//        sax_from_ts(ts->ts, sax->sax);
//        sax_print_bit(sax->sax);

        heap->pq.pop();
    }

    return 0;
}
