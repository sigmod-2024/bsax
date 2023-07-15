



#include "TopKHeap.h"


void TopKHeap::push_ans_approximate(float dis, TS** ts) {
    if (pq.size() < k) {
        pq.push({dis, *ts});
        *ts = *ts + 1;
    }
    else {
        if (pq.top().first > dis) {
            TS * old_ts = pq.top().second;
            pq.pop();
            pq.push({dis, *ts});
            *ts = old_ts;
        }
    }
}

void TopKHeap::push_ans_approximate(float dis, u_int64_t p) {
    if (pq.size() < k) {
        pq.push({dis, (TS*)p});
    }
    else {
        if (pq.top().first > dis) {

            pq.pop();
            pq.push({dis, (TS*)p});
        }
    }
}

void TopKHeap::push_ans_exact(float dis, TS** ts, float& bsf) {
    if (bsf >= dis) {
        if (pq.size() < k) {
            pq.push({dis, *ts});

            if (pq.size() == k && bsf > dis) bsf = pq.top().first;
        }
        else {
            TS * old_ts = pq.top().second;
            pq.pop();
            pq.push({dis, *ts});

            bsf = pq.top().first;
        }
    }
}

void TopKHeap::push_ans_exact(float dis, u_int64_t p, float& bsf) {

    if (pq.top().first > dis) {

        pq.pop();
        pq.push({dis, (TS*)p});
        bsf = pq.top().first;
    }
}


bool TopKHeap::check_approximate(float dis) const {
    if (pq.size() < k) {
        return true;
    }
    else {
        return pq.top().first > dis;
    }
}

bool TopKHeap::check_exact(float dis, float bsf) const {
    if (dis > bsf) return false;
    return true;
}
