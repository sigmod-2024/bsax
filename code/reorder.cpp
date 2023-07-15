#include <vector>
#include <algorithm>
#include "sax.h"

using namespace std;

int main() {

#if dataset_type == 0
    const char *filename = "../../dataset/ts.bin";
    const char *sax_filename = "../../dataset/ts_sax.bin";
//    const char *query_filename = "../../dataset/ts_query.bin";
//    const char *query_filename = "../../dataset/ts.bin";

    const char *filename_re = "../../dataset/ts.bin1";
    const char *sax_filename_re = "../../dataset/ts_sax.bin1";

#elif dataset_type == 1
    const char *filename = "../../dataset/shift1e8.bin";
    const char *sax_filename = "../../dataset/shift1e8_sax.bin";
    const char *query_filename = "../../dataset/shift1b_query_var005.bin";
//    const char *query_filename = "../../dataset/shift1e8.bin";

    const char *filename_re = "../../dataset/shift1e8.bin1";
    const char *sax_filename_re = "../../dataset/shift1e8_sax.bin1";

#elif dataset_type == 2
    const char *filename = "../../dataset/new_deep1b.bin";
    const char *sax_filename = "../../dataset/new_deep1b_sax.bin";
    const char *query_filename = "../../dataset/deep1b_query_var005.bin";

    const char *filename_re = "../../dataset/deep1b_99974688.bin";
    const char *sax_filename_re = "../../dataset/deep1b_sax_99974688.bin";
#elif dataset_type == 3
    const char *filename = "../../dataset/glove.100d.bin";
    const char *sax_filename = "../../dataset/glove_sax.100d.bin";
    const char *query_filename = "../../dataset/glove_query.100d.bin";

    const char *filename_re = "../../dataset/glove.100d.bin1";
    const char *sax_filename_re = "../../dataset/glove_sax.100d.bin1";
#endif

    long long num = 20000000;

    vector<TS> re_ts(num);
    FILE * data_file = fopen(filename,"r");
    FILE * data_file_re = fopen(filename_re,"w");


    for(int i=0;i<num;i++) {
        fread(&re_ts[i], sizeof(TS), 1, data_file);
    }



    TS tmp_ts;
    for(int i=0;i<num*4;i++) {
        fread(&tmp_ts, sizeof(TS), 1, data_file);
        fwrite(&tmp_ts, sizeof(TS), 1, data_file_re);
    }

    for(int i=0;i<num;i++) {
        fwrite(&re_ts[i], sizeof(TS), 1, data_file_re);
    }
//
//    for(int i=0;i<half_num_2;i++) {
//        fread(&tmp_ts, sizeof(TS), 1, data_file);
//        fwrite(&tmp_ts, sizeof(TS), 1, data_file_re);
//    }
//
//
//    vector<SAX> re_sax(20000000);
//    FILE * sax_file = fopen(sax_filename,"r");
//    FILE * sax_file_re = fopen(sax_filename_re,"w");
//
//    for(int i=0;i<20000000;i++) {
//        fread(&re_sax[i], sizeof(SAX), 1, sax_file);
//    }
//
//    num = TOTAL_TS - 20000000;
//    half_num = num / 2;
//    half_num_2 = num - half_num;
//
//
//    SAX tmp_sax;
//    for(int i=0;i<half_num;i++) {
//        fread(&tmp_sax, sizeof(SAX), 1, sax_file);
//        fwrite(&tmp_sax, sizeof(SAX), 1, sax_file_re);
//    }
//
//    for(int i=0;i<20000000;i++) {
//        fwrite(&re_sax[i], sizeof(SAX), 1, sax_file_re);
//    }
//
//    for(int i=0;i<half_num_2;i++) {
//        fread(&tmp_sax, sizeof(SAX), 1, sax_file);
//        fwrite(&tmp_sax, sizeof(SAX), 1, sax_file_re);
//    }


}
