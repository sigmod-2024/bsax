



#ifndef BSAX_TOPKHEAP_H
#define BSAX_TOPKHEAP_H

#include <vector>
#include <queue>
#include "globals.h"
using namespace std;

class TopKHeap {
public:
    TopKHeap(int _k): k(_k) {}

    bool check_approximate(float dis) const;
    bool check_exact(float dis, float bsf) const;
    void push_ans_approximate(float dis, TS** ts);
    void push_ans_approximate(float dis, u_int64_t p);
    void push_ans_exact(float dis, TS** ts, float &bsf);
    void push_ans_exact(float dis, u_int64_t p, float& bsf);

    int k;
    priority_queue<pair<float, TS*>, vector<pair<float, TS*>>> pq;

};


#endif 
