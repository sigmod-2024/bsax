
#include <iostream>
#include <vector>
#include <queue>
#include <algorithm>
#include "globals.h"
#include "sax.h"

int main() {
    const char *filename = "../../dataset/ts.bin";
    const char *filename1 = "../../dataset/ts_query_new.bin";

    FILE* data_file = fopen(filename,"r");
    FILE* data_file_new = fopen(filename1,"w");

    TS ts;
    for(int i=0;i<TOTAL_TS;i++) {
        fread(&ts, sizeof(TS), 1, data_file);
        if(i>90000000 && i % 1000000 == 0) {
            fwrite(&ts, sizeof(TS), 1, data_file_new);
        }
    }





}
