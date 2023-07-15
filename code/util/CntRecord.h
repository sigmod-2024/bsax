



#ifndef BSAX_CNTRECORD_H
#define BSAX_CNTRECORD_H


unsigned long long int count_approximate_ans;
unsigned long long int total_count_approximate_ans;
unsigned long long int count_exact_ans;
unsigned long long int total_count_exact_ans;


unsigned long long int count_node_compute_min_dis;
unsigned long long int total_count_node_compute_min_dis;

unsigned long long int count_approximate_read_ts;
unsigned long long int total_count_approximate_read_ts;
unsigned long long int count_exact_read_ts;
unsigned long long int total_exact_count_read_ts;

#define COUNT_APPROXIMATE_ANS(cnt)  count_approximate_ans += cnt, total_count_approximate_ans += cnt;
#define COUNT_EXACT_ANS(cnt)  count_exact_ans += cnt, total_count_exact_ans += cnt;

#define COUNT_NODE_COMPUTE_MIN_DIS(cnt)  count_node_compute_min_dis += cnt, total_count_node_compute_min_dis += cnt;

#define COUNT_APPROXIMATE_READ_TS(cnt)  count_approximate_read_ts += cnt, total_count_approximate_read_ts += cnt;
#define COUNT_EXACT_READ_TS(cnt)  count_exact_read_ts += cnt, total_exact_count_read_ts += cnt;


#define COUNT_CLEAR count_approximate_ans = 0, count_exact_ans = 0, count_node_compute_min_dis = 0, count_approximate_read_ts = 0, count_exact_read_ts = 0;

#define PRINT_COUNT printf("number of approximate index ans: %llu \t number of exact index ans: %llu\n\
number of exact search compute min dis: %llu\n\
number of approximate search access raw data: %llu \t number of exact search access raw data: %llu\n\n",             \
count_approximate_ans, count_exact_ans,  \
count_node_compute_min_dis, \
count_approximate_read_ts, count_exact_read_ts);


#define PRINT_AVG_COUNT printf("avg number of approximate index ans:%lf \t avg number of exact index ans:%lf\n\
avg number of exact search compute min dis:%lf\n\
avg number of approximate search access raw data:%lf \t avg number of exact search access raw data:%lf\n\n",             \
(double)total_count_approximate_ans/ ((u_int64_t)NUM_SEARCH * TOTAL_TS), (double)total_count_exact_ans/ ((u_int64_t)NUM_SEARCH * TOTAL_TS), \
(double)total_count_node_compute_min_dis / ((u_int64_t)NUM_SEARCH * TOTAL_TS), \
(double)total_count_approximate_read_ts / ((u_int64_t)NUM_SEARCH * TOTAL_TS), (double)total_exact_count_read_ts / ((u_int64_t)NUM_SEARCH * TOTAL_TS));

#endif 
