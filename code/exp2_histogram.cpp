#include <iostream>
#include <vector>
#include <queue>
#include <algorithm>
#include <cstring>
using namespace std;
int main() {
//    const char *filename = "../../dataset_answers/ts_btree_isax.bin";
//    const char *filename = "../../dataset_answers/glove_bbinary_bsax.bin";
//    const char *filename = "../../dataset_answers/glove_btree_bsax.bin";
//    const char *filename = "../../dataset_answers/deep1b_ibinary_isax.bin";
    const char *filename = "../../dataset_answers/deep1b_btree_isax.bin";
//    const char *filename = "../../dataset_answers/shift1e8_ibinary_bsax.bin";
//    const char *filename = "../../dataset_answers/shift1e8_btree_isax.bin";
    const int num = 5;
    float b[] = {0.2, 0.4, 0.6, 0.8, 1};
//    float b[] = {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1};
    FILE *data_file = fopen(filename, "r");
    if (!data_file) {
        cout << "can't find" << filename << endl;
        exit(-1);
    }
    long long all = 99974700;
    vector<float> f(all);
    fread(f.data(), 4, f.size(), data_file);
    long long res[num];
    memset(res, 0, sizeof res);
    for (auto item: f) {
        for (int i = 0; i < num; i++) {
            if (item <= b[i]) {
                res[i]++;
                break;
            }
        }
    }
    double sum = 0;
    for (long long re: res) {
        cout << (double) re / all << endl;
        sum += (double) re / all;
    }
    cout << sum;
    return 0;
}
