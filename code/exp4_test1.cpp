

#include <vector>
#include <algorithm>
#include "util/TopKHeap.h"
#include "sax.h"
#include "CntRecord.h"
#include "TimeRecord.h"
#include <chrono>
using namespace std;

int main() {
    const char *filename1 = "../../dataset_p/ts.bin";
    FILE * file1;
    file1 = fopen(filename1, "r");
    u_int64_t *pp = (u_int64_t*)malloc(8*801088);
    for(int kkk=0;kkk<801088;kkk++) fread(&pp[kkk], 8, 1, file1);

    const char *filename = "../../dataset/ts.bin";
    FILE * file;
    file = fopen(filename, "r");
    float* qts = (float *)malloc(sizeof(float )* TS_LENGTH);
    float true_dis;
    std::chrono::steady_clock::time_point t1;
    std::chrono::steady_clock::time_point t2;
    t1 = std::chrono::steady_clock::now();
    for(int i=0;i<801088;i++) {

        fseek(file, pp[i] * 1024, SEEK_SET);
        fread(qts, sizeof(ts_type), TS_LENGTH, file);
//        COUNT_EXACT_READ_TS(1)
//        float true_dis = ts_euclidean_distance(search_ts->ts, read_ts->ts);
//        true_dis = my_ts_euclidean_distance(search_ts->ts, read_ts->ts, 1);
        true_dis = qts[0];
    }

    t2 = std::chrono:: steady_clock::now();
    std::cout <<""<< std::chrono::duration_cast<std::chrono::milliseconds>( t2-t1 ).count() <<"ms==========================================="<< std::endl;
}
