#include "bitset"
#include "config.h"
#include "iostream"
#include <cassert>
#include "cstring"
#include "sax_bsearch.h"
#include "immintrin.h"
#ifndef isax_globals_h
#define isax_globals_h



#define istime 0

#define iscount_saxt_for_exact 0

#define isgreed 1

#define isprint 1

#define lookupi 2

#define init_st 1

#define cha 0

#define Get_div 20

///// TYPES /////
#if isprint
#define out(a) std::cout<<a<<std::endl
#define out1(a,b) std::cout<<a<<" "<<to_string(b)<<std::endl
#else
#define out(a) //std::cout<<a<<std::endl
#define out1(a,b) //std::cout<<a<<" "<<to_string(b)<<std::endl
#endif


#define CARDINALITY 256
#define BIT_CARDINALITY 8

#if BIT_CARDINALITY == 8
typedef unsigned char sax_type;
#endif


#define build_batch_size 20
#define is_build_m 0
#define is_query_m 1
#define is_reorder_m 1


#define dataset_type 3
#define exp5_use_tree 1
#define btree_use_bsax 0
#if dataset_type
#define daxiao 1
#else
#define daxiao 0
#endif



#if daxiao
typedef unsigned int saxt_type;
#else
typedef unsigned short saxt_type;
#endif
typedef float ts_type;
typedef time_t ts_time;

typedef unsigned char cod;



#if  dataset_type == 3
#define TS_LENGTH 96
#define TOTAL_TS 1000000
#elif dataset_type == 2
#define TS_LENGTH 96
#define TOTAL_TS 99974603  // 
#elif dataset_type == 1
#define TS_LENGTH 128
#define TOTAL_TS 100000000  // 
//#define TOTAL_TS 1000000  // 
#else
#define TS_LENGTH 256
#define TOTAL_TS 100000000  // 
//#define TOTAL_TS 1000000  // 
#endif

#define READ_TS_BATCH 1000000


#define LEAF_MAX_NUM 100

#if daxiao
#define SEGMENTS 32
#else
#define SEGMENTS 16

#endif
#define NUM_PER_SEGMENT (TS_LENGTH / SEGMENTS)

#if dataset_type
#define K 1
#else
#define K 100  // k
#endif

#if dataset_type == 3
#define num_approximate_search_nodes 10
#else
#define num_approximate_search_nodes 500 // b+
#endif


#define NUM_SEARCH 100



#define binary_tree_root_full 0 // 2^segment
#define num_approximate_search_key (num_approximate_search_nodes * LEAF_MAX_NUM) // key
#define b_binary_use_breakpoint_to_split 0

//#define sort_strategy 2
//#define sort_batch_num 10

#define sigtree_leaf_keys_use_vector 1
#define sigtree_delay_create_leafnode 1


#define exp2_binary_use_isax 0

#define exp3_use_tree 0
#define exp3_b_binary_use_isax_to_prune 0
#define exp3_i_binary_use_bsax_to_prune 0
#define exp3_write_index_ans 1


#if is_query_m

#if dataset_type != 0 && (exp5_use_tree == 0 || (exp5_use_tree == 2 && btree_use_bsax == 1))
#define sort_strategy 1 // 0p,1,2(SBB),batchp(SBS)
#define sort_batch_num 5 // 2,batch
#define is_use_m 1
#else
#define sort_strategy 0
#define sort_batch_num 2
#define is_use_m 0
#endif

#else


#if exp5_use_tree == 0 && dataset_type == 0
#define sort_strategy 2 // 0p,1,2(SBB),batchp(SBS)
#define sort_batch_num 10 // 2,batch
#elif exp5_use_tree == 2 && dataset_type == 0 && btree_use_bsax == 1
#define sort_strategy 2 // 0p,1,2(SBB),batchp(SBS)
#define sort_batch_num 10 // 2,batch
#elif exp5_use_tree == 0
#define sort_strategy 0 // 0p,1,2(SBB),batchp(SBS)
#elif exp5_use_tree == 2 && btree_use_bsax == 1
#define sort_strategy 1 // 0p,1,2(SBB),batchp(SBS)
#else
#define sort_strategy 0 // 0p,1,2(SBB),batchp(SBS)
#define sort_batch_num 2 // 2,batch
#endif

#endif


static const int sax_offset = ((CARDINALITY - 1) * (CARDINALITY - 2)) / 2;
static int sax_offset_i[BIT_CARDINALITY + 1] = {0, 0, 3, 21, 105, 465, 1953, 8001, 32385};
static int cardinality_1_i[BIT_CARDINALITY + 1] = {0, 1, 3, 7, 15, 31, 63, 127, 255};

struct TS{
  ts_type ts[TS_LENGTH];
};

/**
 * sax，，>256，，()
 */
struct SAX {
    void set_min_value() {
        memset(sax, 0, sizeof(sax));
    }
    void set_max_value() {
        memset(sax, 0xff, sizeof(sax));
    }

    bool operator== (const SAX& a) const {
        for (int i = SEGMENTS - 1; i >= 0; i--) {
            if (sax[i] == a.sax[i]) continue;
            return false;
        }
        return true;
    }

    sax_type sax[SEGMENTS];
};

struct CARD {
    u_int8_t card[SEGMENTS];
    void set_min_card() {
        memset(card, 0, sizeof(card));
    }
    void set_max_card() {
        for (int i = 0; i < SEGMENTS; i ++ ) {
            card[i] = BIT_CARDINALITY;
        }
    }
};

/**
 * saxt，()，8
**/
struct SAXT {
  saxt_type saxt[BIT_CARDINALITY];

//    SAXT() {}
//    SAXT(const void* saxt_) {
//    memcpy(saxt, saxt_, sizeof(SAXT));
//  }

  bool operator< (const SAXT& a) const {
#if daxiao
    for(int i=BIT_CARDINALITY-1;i>=0;i--) {
        if(saxt[i] == a.saxt[i]) continue;
        return saxt[i] < a.saxt[i];
    }
    return false;
#else
    return *(((uint64_t*)saxt) + 1) == *(((uint64_t*)a.saxt) + 1) ?
           *(uint64_t*)saxt < *(uint64_t*)a.saxt : *(((uint64_t*)saxt) + 1) < *(((uint64_t*)a.saxt) + 1);
#endif
  }
  bool operator> (const SAXT& a) const {
#if daxiao
      for(int i=BIT_CARDINALITY-1;i>=0;i--) {
          if(saxt[i] == a.saxt[i]) continue;
          return saxt[i] > a.saxt[i];
      }
      return false;
#else
    return *(((uint64_t*)saxt) + 1) == *(((uint64_t*)a.saxt) + 1) ?
           *(uint64_t*)saxt > *(uint64_t*)a.saxt : *(((uint64_t*)saxt) + 1) > *(((uint64_t*)a.saxt) + 1);
#endif
  }

  bool operator<= (const SAXT& a) const {
#if daxiao
      for(int i=BIT_CARDINALITY-1;i>=0;i--) {
          if(saxt[i] == a.saxt[i]) continue;
          return saxt[i] < a.saxt[i];
      }
      return true;
#else
    return *(((uint64_t*)saxt) + 1) == *(((uint64_t*)a.saxt) + 1) ?
           *(uint64_t*)saxt <= *(uint64_t*)a.saxt : *(((uint64_t*)saxt) + 1) < *(((uint64_t*)a.saxt) + 1);
#endif
  }
  bool operator>= (const SAXT& a) const {
#if daxiao
      for(int i=BIT_CARDINALITY-1;i>=0;i--) {
          if(saxt[i] == a.saxt[i]) continue;
          return saxt[i] > a.saxt[i];
      }
      return true;
#else
    return *(((uint64_t*)saxt) + 1) == *(((uint64_t*)a.saxt) + 1) ?
           *(uint64_t*)saxt >= *(uint64_t*)a.saxt : *(((uint64_t*)saxt) + 1) > *(((uint64_t*)a.saxt) + 1);
#endif
  }

  bool operator== (const SAXT& a) const {
#if daxiao
      for(int i=BIT_CARDINALITY-1;i>=0;i--) {
          if(saxt[i] == a.saxt[i]) continue;
          return false;
      }
      return true;
#else
    return *(((uint64_t*)saxt) + 1) == *(((uint64_t*)a.saxt) + 1) && *(uint64_t*)saxt == *(uint64_t*)a.saxt;
#endif
  }

};

typedef struct {
  ts_type apaa[SEGMENTS];
} paa_only;

typedef struct {
  ts_type ts[TS_LENGTH];
#if istime
  ts_time tsTime;
#endif
} tsKey;



typedef struct {
  ts_type ts[TS_LENGTH];
#if istime
  ts_time startTime;
  ts_time endTime;
#endif
} aquery_rep;

typedef struct {
  aquery_rep rep;
  int k;
  ts_type paa[SEGMENTS];
  SAXT asaxt;
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

static const size_t send_size1 = 1 +sizeof(int)*2 + sizeof(uint64_t) + sizeof(SAXT) * 2 + sizeof(ts_time) * 2;
static const size_t send_size2 = 1+sizeof(int)*3;
static const size_t send_size2_add = sizeof(uint64_t) + sizeof(SAXT) * 2 + sizeof(ts_time) * 2;

static const size_t sizeinfo_pos = sizeof(aquery_rep) + sizeof(int)*2 + sizeof(float);

static const size_t to_find_size_leafkey = sizeof(aquery_rep) + sizeof(int)*3 + sizeof(float);

static inline int compare_saxt(const void* a, const void* b) {
  if (*(SAXT*)a < *(SAXT*)b) return -1;
  if (*(SAXT*)a > *(SAXT*)b) return 1;
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

struct DisP {

    DisP(float dis, u_int64_t p): dis(dis), p(p) {}
    DisP() {}

    float dis;
    u_int64_t p;
};

struct DisNode {
    DisNode(float dis, u_int64_t begin, u_int64_t size): dis(dis), begin(begin), size(size) {}
    DisNode() {}

    float dis;
    u_int64_t begin;
    u_int64_t size;
};


static bool DisCmp(DisP& a, DisP& b) {
    return a.dis < b.dis;
}

static bool PCmp(DisP& a, DisP& b) {
    return a.p < b.p;
}

static bool DisNodeCmp(DisNode& a, DisNode& b) {
    return a.dis < b.dis;
}




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
