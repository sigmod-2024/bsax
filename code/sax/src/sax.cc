#include "../include/config.h"
#include "../include/globals.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#ifdef VALUES
	#include <values.h>
#include <iostream>

#endif

#include "../include/sax.h"
#include "../include/ts.h"
#include "../include/sax_breakpoints.h"
#include "algorithm"
#include "immintrin.h"
#include "sax_bsearch.h"

///////////////////////////////////////////////////////////////////////////////////////
/** 
 Calculate paa.
 */
void paa_from_ts_simd (ts_type *ts_in, ts_type *paa_out) {

    __m256 xfsSum;
    __m256 a;
    for (int i=0; i < SEGMENTS; i++) {
        xfsSum = _mm256_setzero_ps();
        int off = i * NUM_PER_SEGMENT;
        for (int j=0; j < NUM_PER_SEGMENT; j+=8) {
            a = _mm256_loadu_ps(ts_in + off + j);
            xfsSum = _mm256_add_ps(xfsSum, a);
        }
        const auto* q = (const float*)&xfsSum;
        paa_out[i] = (q[0] + q[1] + q[2] + q[3] + q[4] + q[5] + q[6] + q[7]) / NUM_PER_SEGMENT;
    }
}


void sax_from_paa_simd(ts_type *paa, sax_type *sax) {

    //printf("FROM %lf TO %lf\n", sax_breakpoints[offset], sax_breakpoints[offset + cardinality - 2]);

  __m256 vmask;
  __m256 vt;
  __m128 vmask1;

  for(int i=0; i < SEGMENTS; i++){

    vt = _mm256_set1_ps(paa[i]);
    vmask = _mm256_cmp_ps(vt, BM.a1, _CMP_GT_OQ);
    int bit_mask = _mm256_movemask_ps(vmask);

    int d = sax_c1[bit_mask];
    vmask = _mm256_cmp_ps(vt, BM.a2[d], _CMP_GT_OQ);
    bit_mask = _mm256_movemask_ps(vmask);

    int d1 = sax_c1[bit_mask];
    vmask1 = _mm_cmp_ps(*(__m128*)&vt, BM.a3[d][d1], _CMP_GT_OQ);
    bit_mask = _mm_movemask_ps(vmask1);

    sax[i] = (d << 5) + (d1 << 2) + sax_c1[bit_mask];
  }

  
}

/**
 This function converts a ts record into its sax representation.
 */
void sax_from_ts_simd(ts_type *ts_in, sax_type *sax_out)
{
    // Create PAA representation
    float paa[SEGMENTS];


  __m256 xfsSum;
  __m256 a;
  for (int i=0; i < SEGMENTS; i++) {
    xfsSum = _mm256_setzero_ps();
    int off = i * NUM_PER_SEGMENT;
    for (int j=0; j < NUM_PER_SEGMENT; j+=8) {
      a = _mm256_loadu_ps(ts_in + off + j);
      xfsSum = _mm256_add_ps(xfsSum, a);
    }
    const auto* q = (const float*)&xfsSum;
    paa[i] = (q[0] + q[1] + q[2] + q[3] + q[4] + q[5] + q[6] + q[7]) / NUM_PER_SEGMENT;
  }
    
    // Convert PAA to SAX
    // Note: Each cardinality has cardinality - 1 break points if c is cardinality
    //       the breakpoints can be found in the following array positions: 
    //       FROM (c - 1) * (c - 2) / 2 
    //       TO   (c - 1) * (c - 2) / 2 + c - 1

    //printf("FROM %lf TO %lf\n", sax_breakpoints[offset], sax_breakpoints[offset + cardinality - 2]);

    __m256 vmask;
    __m256 vt;
    __m128 vmask1;

    for(int i=0; i < SEGMENTS; i++){

      vt = _mm256_set1_ps(paa[i]);
      vmask = _mm256_cmp_ps(vt, BM.a1, _CMP_GT_OQ);
      int bit_mask = _mm256_movemask_ps(vmask);

      int d = sax_c1[bit_mask];
      vmask = _mm256_cmp_ps(vt, BM.a2[d], _CMP_GT_OQ);
      bit_mask = _mm256_movemask_ps(vmask);

      int d1 = sax_c1[bit_mask];
      vmask1 = _mm_cmp_ps(*(__m128*)&vt, BM.a3[d][d1], _CMP_GT_OQ);
      bit_mask = _mm_movemask_ps(vmask1);

      sax_out[i] = (d << 5) + (d1 << 2) + sax_c1[bit_mask];
    }
    //sax_print_bit(sax_out, segments, cardinality);
}



void saxt_from_ts_simd(ts_type *ts_in, saxt_type *saxt_out) {
    // Create PAA representation
    float paa[SEGMENTS];
    sax_type sax_out[SEGMENTS];

    __m256 xfsSum;
    __m256 a;
    for (int i=0; i < SEGMENTS; i++) {
      xfsSum = _mm256_setzero_ps();
      int off = i * NUM_PER_SEGMENT;
      for (int j=0; j < NUM_PER_SEGMENT; j+=8) {
        a = _mm256_loadu_ps(ts_in + off + j);
        xfsSum = _mm256_add_ps(xfsSum, a);
      }
      const auto* q = (const float*)&xfsSum;
      paa[i] = (q[0] + q[1] + q[2] + q[3] + q[4] + q[5] + q[6] + q[7]) / NUM_PER_SEGMENT;
    }

    // Convert PAA to SAX
    // Note: Each cardinality has cardinality - 1 break points if c is cardinality
    //       the breakpoints can be found in the following array positions:
    //       FROM (c - 1) * (c - 2) / 2
    //       TO   (c - 1) * (c - 2) / 2 + c - 1
    //printf("FROM %lf TO %lf\n", sax_breakpoints[offset], sax_breakpoints[offset + cardinality - 2]);

    __m256 vmask;
    __m256 vt;
    __m128 vmask1;

    for(int i=0; i < SEGMENTS; i++){

      vt = _mm256_set1_ps(paa[i]);
      vmask = _mm256_cmp_ps(vt, BM.a1, _CMP_GT_OQ);
      int bit_mask = _mm256_movemask_ps(vmask);

      int d = sax_c1[bit_mask];
      vmask = _mm256_cmp_ps(vt, BM.a2[d], _CMP_GT_OQ);
      bit_mask = _mm256_movemask_ps(vmask);

      int d1 = sax_c1[bit_mask];
      vmask1 = _mm_cmp_ps(*(__m128*)&vt, BM.a3[d][d1], _CMP_GT_OQ);
      bit_mask = _mm_movemask_ps(vmask1);

      sax_out[i] = (d << 5) + (d1 << 2) + sax_c1[bit_mask];
    }

//    uint64_t t = (*(uint64_t*)sax_out);
//    (*(uint64_t*)saxt_out) = ((t & 0x8080808080808080) * 0x2040810204081) & 0xff00000000000000 |
//                              (((((t) & 0x4040404040404040) * 0x2040810204081) & 0x7f80000000000000) >> 7) |
//                              (((((t) & 0x2020202020202020) * 0x2040810204081) & 0x3fc0000000000000) >> 14) |
//                              (((((t) & 0x1010101010101010) * 0x2040810204081) & 0x1fe0000000000000) >> 21) |
//                              (((((t) & 0x0808080808080808) * 0x2040810204081) & 0x0ff0000000000000) >> 28) |
//                              (((((t) & 0x0404040404040404) * 0x2040810204081) & 0x07f8000000000000) >> 35) |
//                              (((((t) & 0x0202020202020202) * 0x2040810204081) & 0x03fc000000000000) >> 42) |
//                              (((((t) & 0x0101010101010101) * 0x2040810204081) & 0x01fe000000000000) >> 49);

#if daxiao
  uint64_t t = (*(uint64_t*)sax_out);
  (*(uint64_t*)saxt_out) = ((t & 0x8080808080808080) * 0x2040810204081) & 0xff00000000000000 |
                           (((((t<<1) & 0x8080808080808080) * 0x2040810204081) & 0xff00000000000000) >> 8) |
                           (((((t<<2) & 0x8080808080808080) * 0x2040810204081) & 0xff00000000000000) >> 16) |
                           (((((t<<3) & 0x8080808080808080) * 0x2040810204081) & 0xff00000000000000) >> 24) |
                           (((((t<<4) & 0x8080808080808080) * 0x2040810204081) & 0xff00000000000000) >> 32) |
                           (((((t<<5) & 0x8080808080808080) * 0x2040810204081) & 0xff00000000000000) >> 40) |
                           (((((t<<6) & 0x8080808080808080) * 0x2040810204081) & 0xff00000000000000) >> 48) |
                           (((((t<<7) & 0x8080808080808080) * 0x2040810204081) & 0xff00000000000000) >> 56);
#else
  for(int i=0; i < BIT_CARDINALITY; i++) {
    saxt_out[i] = _mm_movemask_epi8(_mm_slli_epi64(_mm_loadu_si128(reinterpret_cast<const __m128i_u *>(sax_out)), BIT_CARDINALITY - i - 1));
  }
#endif

    
}

void paa_saxt_from_ts_simd(ts_type *ts_in, saxt_type *saxt_out, ts_type *paa) {
  // Create PAA representation
  sax_type sax_out[SEGMENTS];

  __m256 xfsSum;
  __m256 a;
  for (int i=0; i < SEGMENTS; i++) {
    xfsSum = _mm256_setzero_ps();
    int off = i * NUM_PER_SEGMENT;
    for (int j=0; j < NUM_PER_SEGMENT; j+=8) {
      a = _mm256_loadu_ps(ts_in + off + j);
      xfsSum = _mm256_add_ps(xfsSum, a);
    }
    const auto* q = (const float*)&xfsSum;
    paa[i] = (q[0] + q[1] + q[2] + q[3] + q[4] + q[5] + q[6] + q[7]) / NUM_PER_SEGMENT;
  }

  // Convert PAA to SAX
  // Note: Each cardinality has cardinality - 1 break points if c is cardinality
  //       the breakpoints can be found in the following array positions:
  //       FROM (c - 1) * (c - 2) / 2
  //       TO   (c - 1) * (c - 2) / 2 + c - 1
  //printf("FROM %lf TO %lf\n", sax_breakpoints[offset], sax_breakpoints[offset + cardinality - 2]);

  __m256 vmask;
  __m256 vt;
  __m128 vmask1;

  for(int i=0; i < SEGMENTS; i++){

    vt = _mm256_set1_ps(paa[i]);
    vmask = _mm256_cmp_ps(vt, BM.a1, _CMP_GT_OQ);
    int bit_mask = _mm256_movemask_ps(vmask);

    int d = sax_c1[bit_mask];
    vmask = _mm256_cmp_ps(vt, BM.a2[d], _CMP_GT_OQ);
    bit_mask = _mm256_movemask_ps(vmask);

    int d1 = sax_c1[bit_mask];
    vmask1 = _mm_cmp_ps(*(__m128*)&vt, BM.a3[d][d1], _CMP_GT_OQ);
    bit_mask = _mm_movemask_ps(vmask1);

    sax_out[i] = (d << 5) + (d1 << 2) + sax_c1[bit_mask];
  }

#if daxiao
  uint64_t t = (*(uint64_t*)sax_out);
  (*(uint64_t*)saxt_out) = ((t & 0x8080808080808080) * 0x2040810204081) & 0xff00000000000000 |
                           (((((t<<1) & 0x8080808080808080) * 0x2040810204081) & 0xff00000000000000) >> 8) |
                           (((((t<<2) & 0x8080808080808080) * 0x2040810204081) & 0xff00000000000000) >> 16) |
                           (((((t<<3) & 0x8080808080808080) * 0x2040810204081) & 0xff00000000000000) >> 24) |
                           (((((t<<4) & 0x8080808080808080) * 0x2040810204081) & 0xff00000000000000) >> 32) |
                           (((((t<<5) & 0x8080808080808080) * 0x2040810204081) & 0xff00000000000000) >> 40) |
                           (((((t<<6) & 0x8080808080808080) * 0x2040810204081) & 0xff00000000000000) >> 48) |
                           (((((t<<7) & 0x8080808080808080) * 0x2040810204081) & 0xff00000000000000) >> 56);
#else
  for(int i=0; i < BIT_CARDINALITY; i++) {
    saxt_out[i] = _mm_movemask_epi8(_mm_slli_epi64(_mm_loadu_si128(reinterpret_cast<const __m128i_u *>(sax_out)), BIT_CARDINALITY - i - 1));
  }
#endif

  
}

void saxt_from_sax_simd(sax_type *sax_in, saxt_type *saxt_out) {
#if daxiao
  uint64_t t = (*(uint64_t*)sax_in);
  (*(uint64_t*)saxt_out) = ((t & 0x8080808080808080) * 0x2040810204081) & 0xff00000000000000 |
                           (((((t<<1) & 0x8080808080808080) * 0x2040810204081) & 0xff00000000000000) >> 8) |
                           (((((t<<2) & 0x8080808080808080) * 0x2040810204081) & 0xff00000000000000) >> 16) |
                           (((((t<<3) & 0x8080808080808080) * 0x2040810204081) & 0xff00000000000000) >> 24) |
                           (((((t<<4) & 0x8080808080808080) * 0x2040810204081) & 0xff00000000000000) >> 32) |
                           (((((t<<5) & 0x8080808080808080) * 0x2040810204081) & 0xff00000000000000) >> 40) |
                           (((((t<<6) & 0x8080808080808080) * 0x2040810204081) & 0xff00000000000000) >> 48) |
                           (((((t<<7) & 0x8080808080808080) * 0x2040810204081) & 0xff00000000000000) >> 56);
#else
  for(int i=0; i < BIT_CARDINALITY; i++) {
    saxt_out[i] = _mm_movemask_epi8(_mm_slli_epi64(_mm_loadu_si128(reinterpret_cast<const __m128i_u *>(sax_in)), BIT_CARDINALITY - i - 1));
  }
#endif
    
}

void sax_from_saxt_simd(saxt_type *saxt_in, sax_type *sax_out) {
//  uint64_t t = 0;
//  memcpy(&t, saxt_in, sizeof(uint64_t));
#if daxiao
  uint64_t t = (*(uint64_t*)saxt_in);
  (*(uint64_t*)sax_out) = ((t & 0x8080808080808080) * 0x2040810204081) & 0xff00000000000000 |
                          (((((t<<1) & 0x8080808080808080) * 0x2040810204081) & 0xff00000000000000) >> 8) |
                          (((((t<<2) & 0x8080808080808080) * 0x2040810204081) & 0xff00000000000000) >> 16) |
                          (((((t<<3) & 0x8080808080808080) * 0x2040810204081) & 0xff00000000000000) >> 24) |
                          (((((t<<4) & 0x8080808080808080) * 0x2040810204081) & 0xff00000000000000) >> 32) |
                          (((((t<<5) & 0x8080808080808080) * 0x2040810204081) & 0xff00000000000000) >> 40) |
                          (((((t<<6) & 0x8080808080808080) * 0x2040810204081) & 0xff00000000000000) >> 48) |
                          (((((t<<7) & 0x8080808080808080) * 0x2040810204081) & 0xff00000000000000) >> 56);
#else
  uint64_t t[2];
  memset(t, 0, sizeof t);
  memcpy(t, saxt_in, sizeof(SAXT));
  for(int i=0; i < SEGMENTS; i++) {
    int daoi = SEGMENTS - i - 1;
    sax_out[i] = ((((t[0] << daoi) & 0x8000800080008000) * 0x200040008001) >> 60) |
             ((((t[1] << daoi) & 0x8000800080008000) * 0x200040008001) >> 56) ;
  }
#endif

    
}



float minidist_paa_to_saxt_simd(const float* paa, saxt_type* saxt_, cod co_d) {
    ts_type distance = 0;
    // TODO: Store offset in index settings. and pass index settings as parameter.

    int offset = sax_offset_i[co_d];
    sax_type sax[SEGMENTS];
    // For each sax record find the break point

#if daxiao
  uint64_t t = 0;
    memcpy(&t, saxt_, sizeof(saxt_type) * co_d);
    (*(uint64_t*)sax) = ((t & 0x8080808080808080) * 0x2040810204081) & 0xff00000000000000 |
                        (((((t<<1) & 0x8080808080808080) * 0x2040810204081) & 0xff00000000000000) >> 8) |
                        (((((t<<2) & 0x8080808080808080) * 0x2040810204081) & 0xff00000000000000) >> 16) |
                        (((((t<<3) & 0x8080808080808080) * 0x2040810204081) & 0xff00000000000000) >> 24) |
                        (((((t<<4) & 0x8080808080808080) * 0x2040810204081) & 0xff00000000000000) >> 32) |
                        (((((t<<5) & 0x8080808080808080) * 0x2040810204081) & 0xff00000000000000) >> 40) |
                        (((((t<<6) & 0x8080808080808080) * 0x2040810204081) & 0xff00000000000000) >> 48) |
                        (((((t<<7) & 0x8080808080808080) * 0x2040810204081) & 0xff00000000000000) >> 56);
#else
  uint64_t t[2];
  memset(t, 0, sizeof t);
  memcpy(t, saxt_, sizeof(saxt_type) * co_d);
  for(int i=0; i < SEGMENTS; i++) {
    int daoi = SEGMENTS - i - 1;
    sax[i] = ((((t[0] << daoi) & 0x8000800080008000) * 0x200040008001) >> 60) |
             ((((t[1] << daoi) & 0x8000800080008000) * 0x200040008001) >> 56) ;
  }
#endif

    ts_type breakpoint_lower[SEGMENTS];
    memset(breakpoint_lower, 0, sizeof breakpoint_lower);
    ts_type breakpoint_upper[SEGMENTS];
    memset(breakpoint_upper, 0, sizeof breakpoint_upper);

    for (int i=0; i < SEGMENTS; i++) {
//        for(int j=co_d-1; j>=0; j--) {
//          sax[i] |= (saxt_[j]>>(i) & 1) << (j);
//        }
        sax_type region = sax[i];

        /*
            int region_lower = v << (c_m - c_c);
            int region_upper = (~((int)MAXFLOAT << (c_m - c_c)) | region_lower);
        */

        breakpoint_lower[i] = region == 0 ? -MAXFLOAT : sax_breakpoints[offset + region - 1];
        breakpoint_upper[i] = region == cardinality_1_i[co_d] ? MAXFLOAT : sax_breakpoints[offset + region];

        //printf("FROM: \n");
        //sax_print_bit(&region_lower, 1, c_m);
        //printf("TO: \n");
        //sax_print_bit(&region_upper, 1, c_m);
        //printf("\n%d.%d is from %d to %d, %lf - %lf\n", v, c_c, region_lower, region_upper,
        //       breakpoint_lower, breakpoint_upper);
//        if (breakpoint_lower[i] > paa[i]) {
//          distance += pow(breakpoint_lower[i] - paa[i], 2);
//        }
//        else if(breakpoint_upper[i] < paa[i]) {
//          distance += pow(breakpoint_upper[i] - paa[i], 2);
//        }

//        else {
//            printf("%lf is between: %lf and %lf\n", paa[i], breakpoint_lower, breakpoint_upper);
//        }
    }
//
//    float distance1 = 0;
#if daxiao
  __m256 dis8 = _mm256_setzero_ps();
  __m256 paa8 = _mm256_loadu_ps(paa);

  __m256 lower8 = _mm256_loadu_ps(breakpoint_lower);
  __m256 tosub = _mm256_sub_ps(lower8, paa8);
  __m256 f = _mm256_mul_ps(tosub, tosub);
  __m256 mask = _mm256_cmp_ps(lower8, paa8, _CMP_GT_OQ);
  dis8 = _mm256_add_ps(dis8, _mm256_and_ps(f, mask));

  __m256 upper8 = _mm256_loadu_ps(breakpoint_upper);
  tosub = _mm256_sub_ps(upper8, paa8);
  f = _mm256_mul_ps(tosub, tosub);
  mask = _mm256_cmp_ps(upper8, paa8, _CMP_LT_OQ);
  dis8 = _mm256_add_ps(dis8, _mm256_and_ps(f, mask));

  const auto* q = (const float*)&dis8;
  distance = q[0] + q[1] + q[2] + q[3] + q[4] + q[5] + q[6] + q[7];

#else
  __m256 dis8 = _mm256_setzero_ps();

  for(int i=0; i < SEGMENTS; i+=8) {
    __m256 paa8 = _mm256_loadu_ps(paa + i);

    __m256 lower8 = _mm256_loadu_ps(breakpoint_lower + i);
    __m256 tosub = _mm256_sub_ps(lower8, paa8);
    __m256 f = _mm256_mul_ps(tosub, tosub);
    __m256 mask = _mm256_cmp_ps(lower8, paa8, _CMP_GT_OQ);
    dis8 = _mm256_add_ps(dis8, _mm256_and_ps(f, mask));

    __m256 upper8 = _mm256_loadu_ps(breakpoint_upper + i);
    tosub = _mm256_sub_ps(upper8, paa8);
    f = _mm256_mul_ps(tosub, tosub);
    mask = _mm256_cmp_ps(upper8, paa8, _CMP_LT_OQ);
    dis8 = _mm256_add_ps(dis8, _mm256_and_ps(f, mask));

  }

  const auto* q = (const float*)&dis8;
  distance = q[0] + q[1] + q[2] + q[3] + q[4] + q[5] + q[6] + q[7];
#endif

//  if (abs(distance - distance1) > 1e6) {
//    out("");
//    out(distance);
//    out(distance1);
//    ts_print((float*)paa, 16);
//    sax_print_bit(sax);
//    exit(1);
//  }

    //distance = ratio_sqrt * sqrtf(distance);
    return NUM_PER_SEGMENT * distance;
}


/**
 * simdtspaa
 */
void paa_from_ts(const ts_type *ts_in, ts_type *paa)
{
    assert(BIT_CARDINALITY == 8);
    for (int i = 0, s = 0; s < SEGMENTS; i += NUM_PER_SEGMENT, s ++) {
        paa[s] = 0;
        for (int j = i; j < i + NUM_PER_SEGMENT; j ++ ) {
            paa[s] += ts_in[j];
        }
        paa[s] /= NUM_PER_SEGMENT;
    }
}

/**
 * simdtssax
 */
void sax_from_ts(const ts_type *ts_in, sax_type *sax_out)
{
    assert(BIT_CARDINALITY == 8);
    float paa[SEGMENTS];
    for (int i = 0, s = 0; s < SEGMENTS; i += NUM_PER_SEGMENT, s ++) {
        paa[s] = 0;
        for (int j = i; j < i + NUM_PER_SEGMENT; j ++ ) {
            paa[s] += ts_in[j];
        }
        paa[s] /= NUM_PER_SEGMENT;
    }

    for (int s = 0; s < SEGMENTS; s ++ ) {
        int l = 0, r = 256;
        while(l < r) {
            int mid = (l + r + 1) >> 1;
            if (sax_a[mid] < paa[s]) l = mid;
            else r = mid - 1;
        }
        sax_out[s] = l;
    }
}


/**
 * paris
 * @param ts_in
 * @param sax_out
 */
int compare(const void *a, const void *b)
{
    float * c = (float *) b - 1;
    if (*(float*)a>*(float*)c && *(float*)a<=*(float*)b) {
        //printf("Found %lf between %lf and %lf\n",*(float*)a,*(float*)c,*(float*)b);
        return 0;
    }
    else if (*(float*)a<=*(float*)c) {
        return -1;
    }
    else
    {
        return 1;
    }
}
void sax_from_ts_paris(ts_type *ts_in, sax_type *sax_out)
{
    // Create PAA representation
    float *paa = static_cast<float *>(malloc(sizeof(float) * SEGMENTS));
    if(paa == NULL) {
        fprintf(stderr,"error: could not allocate memory for PAA representation.\n");
    }

    int s, i;
    for (s=0; s<SEGMENTS; s++) {
        paa[s] = 0;
        for (i=0; i<NUM_PER_SEGMENT; i++) {
            paa[s] += ts_in[(s * NUM_PER_SEGMENT)+i];
        }
        paa[s] /= NUM_PER_SEGMENT;
        //#ifdef DEBUG
        //printf("%d: %lf\n", s, paa[s]);
        //#endif
    }

    // Convert PAA to SAX
    // Note: Each cardinality has cardinality - 1 break points if c is cardinality
    //       the breakpoints can be found in the following array positions:
    //       FROM (c - 1) * (c - 2) / 2
    //       TO   (c - 1) * (c - 2) / 2 + c - 1
    int cardinality = 1 << BIT_CARDINALITY;
    int offset = ((cardinality - 1) * (cardinality - 2)) / 2;
    //printf("FROM %lf TO %lf\n", sax_breakpoints[offset], sax_breakpoints[offset + cardinality - 2]);

    int si;
    for (si=0; si<SEGMENTS; si++) {
        sax_out[si] = 0;

        // First object = sax_breakpoints[offset]
        // Last object = sax_breakpoints[offset + cardinality - 2]
        // Size of sub-array = cardinality - 1

        float *res = (float *) bsearch(&paa[si], &sax_breakpoints[offset], cardinality - 1,
                                       sizeof(ts_type), compare);
        if(res != NULL)
            sax_out[si] = (int) (res -  &sax_breakpoints[offset]);
        else if (paa[si] > 0)
            sax_out[si] = cardinality-1;
    }

    //sax_print_bit(sax_out, segments, cardinality);
    free(paa);
}

/**
 * simdsaxtsax
 */
void sax_from_saxt(const saxt_type *saxt_in, sax_type *sax_out)
{
    for (int i = 0; i < SEGMENTS; i ++ ) {
        sax_type a_segment_sax = 0;
        saxt_type bit = 1 << (SEGMENTS - i - 1);
        for (int j = 0; j < BIT_CARDINALITY; j ++ ) {
            a_segment_sax |= ( ((saxt_in[j] & bit) >> (SEGMENTS - i - 1)) << j );
        }
        sax_out[SEGMENTS - i - 1] = a_segment_sax;
    }
}


/**
 * simdsaxsaxt
 */
void saxt_from_sax(const sax_type *sax_in, saxt_type *saxt_out)
{
    for (int i = 0; i < BIT_CARDINALITY; i ++ ) {
        saxt_type a_word_saxt = 0;
        sax_type bit = 1 << (BIT_CARDINALITY - i - 1);
        for (int j = 0; j < SEGMENTS; j ++ ) {
            a_word_saxt |= ( ((sax_in[j] & bit) >> (BIT_CARDINALITY - i - 1)) << j);
        }
        saxt_out[BIT_CARDINALITY - i - 1] = a_word_saxt;
    }
}

/**
 * simdtssaxt
 */
void saxt_from_ts(const ts_type *ts_in, saxt_type *saxt_out)
{
    assert(BIT_CARDINALITY == 8);
    float paa[SEGMENTS];
    for (int i = 0, s = 0; s < SEGMENTS; i += NUM_PER_SEGMENT, s ++) {
        paa[s] = 0;
        for (int j = i; j < i + NUM_PER_SEGMENT; j ++ ) {
            paa[s] += ts_in[j];
        }
        paa[s] /= NUM_PER_SEGMENT;
    }

    SAX sax;
    for (int s = 0; s < SEGMENTS; s ++ ) {
        int l = 0, r = 256;
        while(l < r) {
            int mid = (l + r + 1) >> 1;
            if (sax_a[mid] < paa[s]) l = mid;
            else r = mid - 1;
        }
        sax.sax[s] = l;
    }


    for (int i = 0; i < BIT_CARDINALITY; i ++ ) {
        saxt_type a_word_saxt = 0;
        sax_type bit = 1 << (BIT_CARDINALITY - i - 1);
        for (int j = 0; j < SEGMENTS; j ++ ) {
            a_word_saxt |= ( ((sax.sax[j] & bit) >> (BIT_CARDINALITY - i - 1)) << j);
        }
        saxt_out[BIT_CARDINALITY - i - 1] = a_word_saxt;
    }
}

///////////////////////////////////////////////////////////////////////////////////////




/**
 * saxbreakpoint
 */
float dist_breakpoint_to_sax(float breakpoint, saxt_type sax_s) {
    assert(BIT_CARDINALITY == 8);
    float dis = 0;
    ts_type breakpoint_lower = sax_a[sax_s];
    ts_type breakpoint_upper = sax_a[sax_s + 1];

    if (breakpoint_lower > breakpoint) {
        dis += breakpoint_lower - breakpoint;
    }
    else if(breakpoint_upper < breakpoint) {
        dis += breakpoint - breakpoint_upper;
    }
    return dis;
}

/**
 * sax,isax,saxBIT_CARDINALITY
 * sax_aBIT_CARDINALITY = 8
 */
float min_dist_paa_to_sax(const ts_type *paa, SAX sax_) {
    assert(BIT_CARDINALITY == 8);
    float dis = 0;
    for (int i = 0; i < SEGMENTS; i ++ ) {
        sax_type region = sax_.sax[i];
        ts_type breakpoint_lower = sax_a[region];
        ts_type breakpoint_upper = sax_a[region + 1];

        if (breakpoint_lower > paa[i]) {
            dis += (breakpoint_lower - paa[i]) * (breakpoint_lower - paa[i]);
        }
        else if(breakpoint_upper < paa[i]) {
            dis += (breakpoint_upper - paa[i]) * (breakpoint_upper - paa[i]);
        }
    }
    return dis * NUM_PER_SEGMENT;
}

/**
 * isaxsax_card，0，1
 *  sax_aBIT_CARDINALITY = 8
 */
float min_dist_paa_to_isax(const ts_type *paa, SAX isax, CARD sax_card) {
    assert(BIT_CARDINALITY == 8);
    float dis = 0;

    for (int i = 0; i < SEGMENTS; i ++ ) {
        sax_type region_lb = isax.sax[i] & (0xff << (BIT_CARDINALITY - sax_card.card[i]));
        sax_type region_ub = isax.sax[i] | (0xff >> sax_card.card[i]);
        ts_type breakpoint_lower = sax_a[region_lb];
        ts_type breakpoint_upper = sax_a[region_ub + 1];

        if (breakpoint_lower > paa[i]) {
            dis += (breakpoint_lower - paa[i]) * (breakpoint_lower - paa[i]);
        }
        else if(breakpoint_upper < paa[i]) {
            dis += (breakpoint_upper - paa[i]) * (breakpoint_upper - paa[i]);
        }
    }
    return dis * NUM_PER_SEGMENT;
}

/**
 * sax，
 * sax_aBIT_CARDINALITY = 8
 */
float min_dist_paa_to_bsax(const ts_type *paa, SAX sax_lb, SAX sax_ub) {
    assert(BIT_CARDINALITY == 8);
    float dis = 0;
    for (int i = 0; i < SEGMENTS; i ++ ) {
        sax_type region_lb = sax_lb.sax[i];
        sax_type region_ub = sax_ub.sax[i];
        ts_type breakpoint_lower = sax_a[region_lb];
        ts_type breakpoint_upper = sax_a[region_ub + 1];

        if (breakpoint_lower > paa[i]) {
            dis += (breakpoint_lower - paa[i]) * (breakpoint_lower - paa[i]);
        }
        else if(breakpoint_upper < paa[i]) {
            dis += (breakpoint_upper - paa[i]) * (breakpoint_upper - paa[i]);
        }
    }
    return dis * NUM_PER_SEGMENT;
}

/**
 * saxtcard，，，01
 * sax_abit_cardinality = 8
 */
float min_dist_paa_to_saxt(const ts_type *paa, SAXT saxt, u_int8_t card) {
    assert(BIT_CARDINALITY == 8);
    SAX sax;
    sax_from_saxt(saxt.saxt, sax.sax);
    CARD card_;
    for (int i = 0; i < SEGMENTS; i ++ ) {
        card_.card[i] = card;
    }
    return min_dist_paa_to_isax(paa, sax, card_);
}

/**
 * saxt，，，01
 * sax_abit_cardinality = 8
 * ，word
 */
float min_dist_paa_to_saxt(const ts_type *paa, SAXT saxt_lb, SAXT saxt_ub) {
    exit(123);
//    assert(BIT_CARDINALITY == 8);
//    int cod_word = 0;   // word
//    int cod_bit = 0;  // word
//    for (int i = BIT_CARDINALITY - 1; i >= 0; i -- ) {
//        if (saxt_lb.saxt[i] == saxt_ub.saxt[i]) cod_word ++;
//        else {
//            for (int j = SEGMENTS - 1; j >= 0; j -- ) {
//                if (((saxt_lb.saxt[i] >> j) & 1) == ((saxt_ub.saxt[i] >> j) & 1)) {
//                    cod_bit ++ ;
//                }
//                else {
//                    break;
//                }
//            }
//            break;
//        }
//    }
//
//    for (int i = 0; i < BIT_CARDINALITY - cod_word - 1; i ++ ) {
//        saxt_lb.saxt[i] = 0;
//        saxt_ub.saxt[i] = 0xff;
//    }
//    if (cod_bit > 0) {
//        saxt_lb.saxt[BIT_CARDINALITY - cod_word - 1] &= (0xff << (BIT_CARDINALITY - cod_bit));
//        saxt_ub.saxt[BIT_CARDINALITY - cod_word - 1] |= (0xff >> cod_bit);
//    }
//
//    sax_type sax_lb[SEGMENTS];
//    sax_type sax_ub[SEGMENTS];
//
//#if daxiao
//    uint64_t sax_lb_tmp = 0;
//    memcpy(&sax_lb_tmp, saxt_lb.asaxt, sizeof(saxt_only));
//    (*(uint64_t*)sax_lb) = ((sax_lb_tmp & 0x8080808080808080) * 0x2040810204081) & 0xff00000000000000 |
//                        (((((sax_lb_tmp<<1) & 0x8080808080808080) * 0x2040810204081) & 0xff00000000000000) >> 8) |
//                        (((((sax_lb_tmp<<2) & 0x8080808080808080) * 0x2040810204081) & 0xff00000000000000) >> 16) |
//                        (((((sax_lb_tmp<<3) & 0x8080808080808080) * 0x2040810204081) & 0xff00000000000000) >> 24) |
//                        (((((sax_lb_tmp<<4) & 0x8080808080808080) * 0x2040810204081) & 0xff00000000000000) >> 32) |
//                        (((((sax_lb_tmp<<5) & 0x8080808080808080) * 0x2040810204081) & 0xff00000000000000) >> 40) |
//                        (((((sax_lb_tmp<<6) & 0x8080808080808080) * 0x2040810204081) & 0xff00000000000000) >> 48) |
//                        (((((sax_lb_tmp<<7) & 0x8080808080808080) * 0x2040810204081) & 0xff00000000000000) >> 56);
//
//    uint64_t sax_ub_tmp = 0;
//    memcpy(&sax_ub_tmp, saxt_ub.asaxt, sizeof(saxt_only));
//    (*(uint64_t*)sax_ub) = ((sax_ub_tmp & 0x8080808080808080) * 0x2040810204081) & 0xff00000000000000 |
//                        (((((sax_ub_tmp<<1) & 0x8080808080808080) * 0x2040810204081) & 0xff00000000000000) >> 8) |
//                        (((((sax_ub_tmp<<2) & 0x8080808080808080) * 0x2040810204081) & 0xff00000000000000) >> 16) |
//                        (((((sax_ub_tmp<<3) & 0x8080808080808080) * 0x2040810204081) & 0xff00000000000000) >> 24) |
//                        (((((sax_ub_tmp<<4) & 0x8080808080808080) * 0x2040810204081) & 0xff00000000000000) >> 32) |
//                        (((((sax_ub_tmp<<5) & 0x8080808080808080) * 0x2040810204081) & 0xff00000000000000) >> 40) |
//                        (((((sax_ub_tmp<<6) & 0x8080808080808080) * 0x2040810204081) & 0xff00000000000000) >> 48) |
//                        (((((sax_ub_tmp<<7) & 0x8080808080808080) * 0x2040810204081) & 0xff00000000000000) >> 56);
//#else
//
//    uint64_t sax_lb_tmp[2];
//    memcpy(sax_lb_tmp, saxt_lb.saxt, sizeof(SAXT));
//    for(int i = 0; i < SEGMENTS; i ++ ) {
//        int trans = SEGMENTS - i - 1;
//        sax_lb[i] = ((((sax_lb_tmp[0] << trans) & 0x8000800080008000) * 0x200040008001) >> 60) |
//                    ((((sax_lb_tmp[1] << trans) & 0x8000800080008000) * 0x200040008001) >> 56) ;
//    }
//
//    uint64_t sax_ub_tmp[2];
//    memcpy(sax_ub_tmp, saxt_ub.saxt, sizeof(SAXT));
//    for(int i = 0; i < SEGMENTS; i ++ ) {
//        int trans = SEGMENTS - i - 1;
//        sax_ub[i] = ((((sax_ub_tmp[0] << trans) & 0x8000800080008000) * 0x200040008001) >> 60) |
//                    ((((sax_ub_tmp[1] << trans) & 0x8000800080008000) * 0x200040008001) >> 56) ;
//    }
//#endif
//    float dis = 0;
//
//    for (int i = 0; i < SEGMENTS; i ++ ) {
//        sax_type region_lb = sax_lb[i];
//        sax_type region_ub = sax_ub[i];
//        ts_type breakpoint_lower = sax_a[region_lb];
//        ts_type breakpoint_upper = sax_a[region_ub + 1];
//
//        if (breakpoint_lower > paa[i]) {
//            dis += (breakpoint_lower - paa[i]) * (breakpoint_lower - paa[i]);
//        }
//        else if(breakpoint_upper < paa[i]) {
//            dis += (breakpoint_upper - paa[i]) * (breakpoint_upper - paa[i]);
//        }
//    }
//    return dis * NUM_PER_SEGMENT;
}



///////////////////////////////////////////////////////////////////////////////////////



void printbin(long long unsigned int n, int size) {
#if isprint
    char *b = static_cast<char *>(malloc(sizeof(char) * (size + 1)));
    int i;

    for (i = 0; i < size; i ++ ) {
        b[i] = '0';
    }

    for (i = 0; i < size; i ++, n /= 2)
        if (n & 1) b[size - 1 - i] = '1';

    b[size] = '\0';
    printf("%s", b);
    free(b);
#endif
}

void serial_printbin (unsigned long long n, int size) {
    char *b = static_cast<char *>(malloc(sizeof(char) * (size + 1)));
    int i;

    for (i=0; i<size; i++) {
        b[i] = '0';
    }

    for (i=0; i<size; i++, n=n/2)
        if (n%2) b[size-1-i] = '1';

    b[size] = '\0';
    printf(" %s ", b);
}


void sax_print(sax_type *sax) {
    int i;
    for (i=0; i < SEGMENTS; i++) {
        std::cout<<(int)sax[i]<<" ";
    }
    printf("\n");
}

void sax_print(SAX sax) {
    sax_print(sax.sax);
}

void sax_print_bit(SAX sax)
{
    int i;
    for (i = 0; i < SEGMENTS; i++) {
        printf("%d:\t", i);
        printbin(sax.sax[i], BIT_CARDINALITY);
        printf("\n");
    }
    printf("\n");
}

void sax_print(sax_type *sax, int segments, int bit_cardinality, CARD* card)
{
    int i;
    for (i = 0; i < segments; i ++) {
        printf("%d:\t", i);
        printbin(sax[i], bit_cardinality);
        printf(" %d\n", card->card[i]);
    }
    printf("\n");
}

void saxt_print_bit(SAXT saxt) {
#if isprint
    int i;
    for (i = 0; i < BIT_CARDINALITY; i ++) {
        printbin(saxt.saxt[i], SEGMENTS);
        printf("\n");
    }
    printf("\n");
#endif
}

void saxt_print(saxt_type *saxt) {
#if isprint
    int i;
    for (i=0; i < BIT_CARDINALITY; i++) {
//        printf("%d:\t", i);
        std::cout<<(int)saxt[i]<<" ";
//        printbin(saxt[i], Segments);
    }
    printf("\n");
#endif
}


void saxt_print(SAXT saxt) {
    saxt_print(saxt.saxt);
}

void leafkey_print(void* saxt) {
    int i;
    for (i=0; i < BIT_CARDINALITY * 2; i++) {
//        printf("%d:\t", i);
        std::cout<<(int)((unsigned char*)saxt)[i]<<" ";
//        printbin(saxt[i], Segments);
    }
    printf("\n");
}

void saxt_print(saxt_type *saxt, saxt_type *prefix, cod co_d) {
#if isprint
    int i;
    for (i=0; i < co_d; i++) {
//        printf("%d:\t", i);
//        printbin(prefix[i], Segments);
        std::cout<<(int)prefix[i]<<" ";
    }
    for (i=co_d; i < BIT_CARDINALITY; i++) {
//        printf("%d:\t", i);
//        printbin(saxt[i-co_d], Segments);
        std::cout<<(int)saxt[i-co_d]<<" ";
    }
    printf("\n");
#endif
}


