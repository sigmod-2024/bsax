#include <iostream>
#include <vector>
#include <queue>
#include <algorithm>
#include "globals.h"
#include "sax.h"


using namespace std;




int main() {
    const char *filename = "../../dataset/ts.bin";
    const char *sax_filename = "../../dataset/ts_sax.bin";
//    const char *filename = "../../dataset/glove.100d.bin";
//    const char *sax_filename = "../../dataset/glove_sax_16.100d.bin";
//    const char *filename = "../../dataset/deep1b.bin";
//    const char *sax_filename = "../../dataset/deep1b_5m_sax.bin";
//    const char *filename = "../../dataset/shift1e8.bin";
//    const char *sax_filename = "../../dataset/shift1e8_sax.bin";

    FILE* data_file = fopen(filename,"r");
    FILE* sax_data_file = fopen(sax_filename,"w");
    if (!data_file) {
        cout << "" << filename << endl;
        exit(-1);
    }
//    vector<SAX> sax_vec(READ_TS_BATCH);
//    vector<TS> ts_vec(READ_TS_BATCH);

    TS ts;
    SAX s;

    for(int i=0;i<TOTAL_TS;i++) {
        fread(&ts, sizeof(TS), 1, data_file);
        sax_from_ts(ts.ts, s.sax);
        fwrite(&s, sizeof(SAX), 1, sax_data_file);
    }


    return 0;
}
