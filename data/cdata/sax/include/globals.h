









#include "bitset"
#include "config.h"
#include "iostream"
#include <cassert>
#include "cstring"
#include "sax_bsearch.h"
#include "immintrin.h"
#ifndef isax_globals_h
#define isax_globals_h


#define daxiao 0 

#define istime 0

#define iscount_saxt_for_exact 0

#define isgreed 1

#define isprint 0

#define lookupi 2

#define init_st 1

#define cha 0

#define isap 1

#define ishash 0


#define Get_div 20

#define info_p_max_size 10000


#if isprint
#define out(a) std::cout<<a<<std::endl
#define out1(a,b) std::cout<<a<<" "<<to_string(b)<<std::endl
#define out2(a) std::cout<<a<<std::endl
#else
#define out(a) 
#define out1(a,b) 
#define out2(a) 
#endif


typedef unsigned char sax_type;
typedef sax_type* sax;


#if daxiao
typedef unsigned char saxt_type;
#else
typedef unsigned short saxt_type;
#endif
typedef saxt_type* saxt;
typedef saxt_type* saxt_prefix;
typedef float ts_type;
typedef ts_type *ts;
typedef time_t ts_time;

typedef unsigned char cod;

#define Cardinality 256
#define Bit_cardinality 8
#if daxiao
#define Segments 8
#define nchuw 32 
#define Ts_values_per_segment 32
#else
#define Segments 16
#define nchuw 16 
#define Ts_values_per_segment 16
#endif


#define Ts_length 256
#define Leaf_maxnum 2048
#define Leaf_minnum Leaf_maxnum/2

#define Leaf_maxnum_rebalance 10
#define Leaf_minnum_rebalance 5


#define init_num 4000000


#define Table_maxnum 2000000


#define pool_size init_num/Table_maxnum


#define get_exact_multiThread_file_size 1000*1024*1024

#define pool_get_size 32  
#define pool_compaction_size 2  


#define compaction_leaf_size 500

#define qiehuan 0

#define input_buffer_size 2048  



static const int Leaf_rebuildnum = Leaf_maxnum * 2;
static const int compaction_buffer_size = Leaf_rebuildnum;




static const int sax_offset = ((Cardinality - 1) * (Cardinality - 2)) / 2;
static int sax_offset_i[Bit_cardinality+1] = {0,0,3,21,105,465,1953,8001,32385};
static int cardinality_1_i[Bit_cardinality+1] = {0,1,3,7,15,31,63,127,255};

typedef struct {
  ts_type ts[Ts_length];
} ts_only;

typedef struct saxt_only_rep{
  saxt_type asaxt[Bit_cardinality];

  saxt_only_rep() {}

  saxt_only_rep(const void* saxt_) {
    memcpy(asaxt, saxt_, sizeof(saxt_only_rep));
  }

  bool operator< (const saxt_only_rep& a) const {
#if daxiao
    return *(uint64_t*)asaxt < *(uint64_t*)a.asaxt;
#else
    return *(((uint64_t*)asaxt)+1) == *(((uint64_t*)a.asaxt)+1) ?
           *(uint64_t*)asaxt < *(uint64_t*)a.asaxt : *(((uint64_t*)asaxt)+1) < *(((uint64_t*)a.asaxt)+1);
#endif
  }
  bool operator> (const saxt_only_rep& a) const {
#if daxiao
    return *(uint64_t*)asaxt > *(uint64_t*)a.asaxt;
#else
    return *(((uint64_t*)asaxt)+1) == *(((uint64_t*)a.asaxt)+1) ?
           *(uint64_t*)asaxt > *(uint64_t*)a.asaxt : *(((uint64_t*)asaxt)+1) > *(((uint64_t*)a.asaxt)+1);
#endif
  }

  bool operator<= (const saxt_only_rep& a) const {
#if daxiao
    return *(uint64_t*)asaxt <= *(uint64_t*)a.asaxt;
#else
    return *(((uint64_t*)asaxt)+1) == *(((uint64_t*)a.asaxt)+1) ?
           *(uint64_t*)asaxt <= *(uint64_t*)a.asaxt : *(((uint64_t*)asaxt)+1) < *(((uint64_t*)a.asaxt)+1);
#endif
  }
  bool operator>= (const saxt_only_rep& a) const {
#if daxiao
    return *(uint64_t*)asaxt >= *(uint64_t*)a.asaxt;
#else
    return *(((uint64_t*)asaxt)+1) == *(((uint64_t*)a.asaxt)+1) ?
           *(uint64_t*)asaxt >= *(uint64_t*)a.asaxt : *(((uint64_t*)asaxt)+1) > *(((uint64_t*)a.asaxt)+1);
#endif
  }

  bool operator== (const saxt_only_rep& a) const {
#if daxiao
    return *(uint64_t*)asaxt == *(uint64_t*)a.asaxt;
#else
    return *(((uint64_t*)asaxt)+1) == *(((uint64_t*)a.asaxt)+1) && *(uint64_t*)asaxt == *(uint64_t*)a.asaxt;
#endif
  }

} saxt_only;

typedef struct {
  ts_type apaa[Segments];
} paa_only;

typedef struct {
  ts_type ts[Ts_length];
#if istime
  ts_time tsTime;
#endif
} tsKey;



typedef struct {
  ts_type ts[Ts_length];
#if istime
  ts_time startTime;
  ts_time endTime;
#endif
} aquery_rep;

typedef struct {
  aquery_rep rep;
  int k;
  ts_type paa[Segments];
  saxt_only asaxt;
} aquery;

typedef struct ares_exact_rep{
  tsKey atskey;
  float dist;

  bool operator< (const ares_exact_rep& a) const {
    return dist < a.dist;
  }
  bool operator> (const ares_exact_rep& a) const {
    return dist > a.dist;
  }
} ares_exact;

typedef struct ares{
  ares_exact rep;
  void* p;

  bool operator< (const ares& a) const {
    return rep < a.rep;
  }
  bool operator> (const ares& a) const {
    return rep > a.rep;
  }
} ares;



typedef std::pair<float, void*> dist_p;

static const size_t send_size1 = 1+sizeof(int)*2+sizeof(uint64_t)+sizeof(saxt_only)*2+sizeof(ts_time)*2;
static const size_t send_size2 = 1+sizeof(int)*3;
static const size_t send_size2_add = sizeof(uint64_t) + sizeof(saxt_only)*2 + sizeof(ts_time)*2;

static const size_t sizeinfo_pos = sizeof(aquery_rep) + sizeof(int)*2 + sizeof(float);

#if isap
static const size_t to_find_size_leafkey = sizeof(aquery_rep) + sizeof(int)*3 + sizeof(float) + sizeof(void*);
#else
static const size_t to_find_size_leafkey = sizeof(aquery_rep) + sizeof(int)*3 + sizeof(float);
#endif


static inline int compare_saxt(const void* a, const void* b) {
  if (*(saxt_only*)a < *(saxt_only*)b) return -1;
  if (*(saxt_only*)a > *(saxt_only*)b) return 1;
  return 0;
}

typedef struct to_bsear_rep {

  to_bsear_rep() {
    a1 = _mm256_loadu_ps(sax_a1);
    for(int i=0;i<8;i++) a2[i] = _mm256_loadu_ps(sax_a2[i]);
    for(int i=0;i<8;i++)
      for(int j=0;j<8;j++)
        a3[i][j] = _mm_loadu_ps(sax_a3[i][j]);
  }
  __m256 a1;
  __m256 a2[8];
  __m128 a3[8][8];
} to_bsear;

static to_bsear BM;



typedef unsigned long long file_position_type;
typedef unsigned long long root_mask_type;

enum response {OUT_OF_MEMORY_FAILURE, FAILURE, SUCCESS};
enum insertion_mode {PARTIAL = 1,
                     TMP = 2,
                     FULL = 4,
                     NO_TMP = 8};

enum buffer_cleaning_mode {FULL_CLEAN, TMP_ONLY_CLEAN, TMP_AND_TS_CLEAN};
enum node_cleaning_mode {DO_NOT_INCLUDE_CHILDREN = 0,
                         INCLUDE_CHILDREN = 1};






























































































































































































































































#endif
