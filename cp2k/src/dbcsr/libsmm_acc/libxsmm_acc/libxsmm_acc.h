/*****************************************************************************
 *  CP2K: A general program to perform molecular dynamics simulations        *
 *  Copyright (C) 2000 - 2016  CP2K developers group                         *
 *****************************************************************************/

//! **************************************************************************
//!> \author Hans Pabst (Intel Corp.)
//! **************************************************************************

#ifndef LIBXSMM_ACC_HPP
#define LIBXSMM_ACC_HPP

#if defined(__LIBXSMM) || (defined(__ACC) && defined(__ACC_MIC) && defined(__DBCSR_ACC) && defined(__LIBXSTREAM))
#include "../include/libsmm_acc.h"

#if defined(__ACC) && defined(__ACC_MIC) && defined(__DBCSR_ACC) && defined(__LIBXSTREAM)
# include "../../../acc/mic/libmicacc.h"
#endif

#if defined(__LIBXSMM)
# include <libxsmm.h>
# if (0 != LIBXSMM_ROW_MAJOR)
#   error Please compile LIBXSMM using "make ROW_MAJOR=0 ..."!
# endif
# define LIBXSMM_ACC_ASSERT(A) assert(A)
# define LIBXSMM_ACC_NOT_SUPPORTED 1
# define LIBXSMM_ACC_ERROR_NONE 0
# define LIBXSMM_ACC_ERROR_CONDITION -2
# define LIBXSMM_ACC_CHECK_CONDITION(CONDITION) if (!(CONDITION)) return LIBXSMM_ACC_ERROR_CONDITION
# define LIBXSMM_ACC_CHECK_CALL_ASSERT(EXPRESSION) EXPRESSION
# define LIBXSMM_ACC_PRAGMA_LOOP_COUNT LIBXSMM_PRAGMA_LOOP_COUNT
# define LIBXSMM_ACC_ALIGN_VALUE LIBXSMM_ALIGN_VALUE
# define LIBXSMM_ACC_HASH2 LIBXSMM_HASH2
# define LIBXSMM_ACC_ASSUME_ALIGNED LIBXSMM_ASSUME_ALIGNED
# define LIBXSMM_ACC_ALIGN_LDST LIBXSMM_ALIGN_LDST
# define LIBXSMM_ACC_ALIGNED_STORES (0 != (LIBXSMM_GEMM_FLAG_ALIGN_C & LIBXSMM_FLAGS) ? LIBXSMM_ALIGNMENT : 1)
# define LIBXSMM_ACC_ALIGNMENT LIBXSMM_ALIGNMENT
# define LIBXSMM_ACC_LOOP_MAX_M LIBXSMM_MAX_M
# define LIBXSMM_ACC_LOOP_AVG_M LIBXSMM_AVG_M
# define LIBXSMM_ACC_LOOP_MAX_N LIBXSMM_MAX_N
# define LIBXSMM_ACC_LOOP_AVG_N LIBXSMM_AVG_N
# define LIBXSMM_ACC_RESTRICT LIBXSMM_RESTRICT
# if defined(LIBXSMM_EXTERN)
#   define LIBXSMM_ACC_EXTERN LIBXSMM_EXTERN
# else /*legacy*/
#   define LIBXSMM_ACC_EXTERN LIBXSMM_EXTERN_C
# endif
# define LIBXSMM_ACC_FSYMBOL LIBXSMM_FSYMBOL
# define LIBXSMM_ACC_UNUSED LIBXSMM_UNUSED
# define LIBXSMM_ACC_MOD2 LIBXSMM_MOD2
# define LIBXSMM_ACC_MAX LIBXSMM_MAX
# define LIBXSMM_ACC_ATTRIBUTE LIBXSMM_ATTRIBUTE
# define LIBXSMM_ACC_RETARGETABLE LIBXSMM_RETARGETABLE
# if defined(LIBXSMM_OFFLOAD_BUILD)
#   define LIBXSMM_ACC_OFFLOAD_BUILD LIBXSMM_OFFLOAD_BUILD
# endif
/** Prefetch strategy used when processing an MM stack. */
LIBXSMM_ACC_EXTERN int libxsmm_acc_prefetch;
#else // defined(__LIBXSTREAM)
# define LIBXSMM_ACC_ASSERT LIBXSTREAM_ASSERT
# define LIBXSMM_ACC_NOT_SUPPORTED LIBXSTREAM_NOT_SUPPORTED
# define LIBXSMM_ACC_ERROR_NONE LIBXSTREAM_ERROR_NONE
# define LIBXSMM_ACC_ERROR_CONDITION LIBXSTREAM_ERROR_CONDITION
# define LIBXSMM_ACC_CHECK_CONDITION LIBXSTREAM_CHECK_CONDITION
# define LIBXSMM_ACC_CHECK_CALL_ASSERT LIBXSTREAM_CHECK_CALL_ASSERT
# define LIBXSMM_ACC_PRAGMA_LOOP_COUNT LIBXSTREAM_PRAGMA_LOOP_COUNT
# define LIBXSMM_ACC_ALIGN_VALUE LIBXSTREAM_ALIGN_VALUE
# define LIBXSMM_ACC_HASH2 LIBXSTREAM_HASH2
# define LIBXSMM_ACC_ASSUME_ALIGNED LIBXSTREAM_ASSUME_ALIGNED
# define LIBXSMM_ACC_ALIGN_LDST(POINTER) (POINTER)
# define LIBXSMM_ACC_ALIGNED_STORES 1
# define LIBXSMM_ACC_ALIGNMENT LIBXSTREAM_MAX_SIMD
# define LIBXSMM_ACC_LOOP_MAX_M LIBXSMM_ACC_MAX_M
# define LIBXSMM_ACC_LOOP_AVG_M 23
# define LIBXSMM_ACC_LOOP_MAX_N LIBXSMM_ACC_MAX_N
# define LIBXSMM_ACC_LOOP_AVG_N 23
# define LIBXSMM_ACC_RESTRICT LIBXSTREAM_RESTRICT
# define LIBXSMM_ACC_EXTERN LIBXSTREAM_EXTERN_C
# define LIBXSMM_ACC_FSYMBOL LIBXSTREAM_FSYMBOL
# define LIBXSMM_ACC_UNUSED LIBXSTREAM_UNUSED
# define LIBXSMM_ACC_MOD2 LIBXSTREAM_MOD2
# define LIBXSMM_ACC_MAX LIBXSTREAM_MAX
# define LIBXSMM_ACC_ATTRIBUTE LIBXSTREAM_ATTRIBUTE
# define LIBXSMM_ACC_RETARGETABLE LIBXSTREAM_RETARGETABLE
# if defined(LIBXSTREAM_OFFLOAD_BUILD)
#   define LIBXSMM_ACC_OFFLOAD_BUILD LIBXSTREAM_OFFLOAD_BUILD
# endif
#endif

#define LIBXSMM_ACC_ABORT(MESSAGE) xsmm_acc_abort(__FILE__, __LINE__, MESSAGE)


/** Upper limits for the supported matrix sizes. */
#define LIBXSMM_ACC_MAX_M 368
#define LIBXSMM_ACC_MAX_N 368
#define LIBXSMM_ACC_MAX_K 368

/** Number of parameters per stack entry. */
#define LIBXSMM_ACC_NPARAMS 7

/** Nested parallelism. */
#if !defined(LIBXSMM_ACC_OPENMP) && ((defined(_OPENMP) && defined(__ACC) && defined(__ACC_MIC) && defined(__DBCSR_ACC) && defined(__LIBXSTREAM)) || defined(__NESTED_OPENMP))
# define LIBXSMM_ACC_OPENMP
#endif

#if !defined(LIBXSMM_ACC_SORT)
/** Disable sort at compile time */
# define LIBXSMM_ACC_SORT
#endif

/**
 * Ensures an amortized synchronization overhead.
 * >=2: maximum number of locally processed MM
 * >=1: assuming unsorted C indexes
 */
#if defined(LIBXSMM_ACC_OPENMP)
# if defined(__MIC__)
#   define LIBXSMM_ACC_NLOCAL 128
# else
#   define LIBXSMM_ACC_NLOCAL 16
# endif
#else
# define LIBXSMM_ACC_NLOCAL 1
#endif

/** OpenMP scheduling policy (and chunk size) */
#if defined(__MIC__)
# define LIBXSMM_ACC_SCHEDULE schedule(dynamic)
#else
# define LIBXSMM_ACC_SCHEDULE
#endif

/**
 * Synchronization mechanism.
 * >1: number of locks (POT)
 * =1: omp critical
 * =0: atomic
 */
#if defined(__MIC__)
# define LIBXSMM_ACC_SYNCHRONIZATION 1
#else
# define LIBXSMM_ACC_SYNCHRONIZATION 16
#endif

#if defined(__RECONFIGURE)
# define LIBXSMM_ACC_MM_DRIVER 4 /*mm_driver_xsmm*/
# define LIBXSMM_ACC_STACKSIZE 1000000
# define LIBXSMM_ACC_MULTREC_LIMIT 64
# define LIBXSMM_ACC_COMM_THREAD_LOAD 99
# if defined(__ACC) && defined(__ACC_MIC) && defined(__DBCSR_ACC) && defined(__LIBXSTREAM)
#   define LIBXSMM_ACC_ACCDRV_POSTERIOR_STREAMS 1
#   define LIBXSMM_ACC_ACCDRV_POSTERIOR_BUFFERS 1
#   define LIBXSMM_ACC_ACCDRV_PRIORITY_STREAMS 1
#   define LIBXSMM_ACC_ACCDRV_PRIORITY_BUFFERS 1
#   define LIBXSMM_ACC_ACCDRV_MIN_MFLOPS_PERSTACK 120
#   define LIBXSMM_ACC_ACCDRV_MIN_NFLOPS_PERMM 0
# endif
#endif // defined(__RECONFIGURE)

/*#define LIBXSMM_ACC_PRETRANSPOSE*/
/*#define LIBXSMM_ACC_MKLTRANS*/

/** Must match stack_descriptor_type. */
struct libxsmm_acc_stackdesc_type {
  int m, n, k, max_m, max_n, max_k;
  int defined_mnk/*boolean*/;
};

union LIBXSMM_ACC_RETARGETABLE libxsmm_acc_param_type {
  struct component_type { int m, n, k, ia, ib, ic, s; } component;
  int param[LIBXSMM_ACC_NPARAMS];
};

enum libxsmm_acc_params {
  LIBXSMM_ACC_PARAM_M = 0, LIBXSMM_ACC_PARAM_N = 1, LIBXSMM_ACC_PARAM_K = 2,
  LIBXSMM_ACC_PARAM_A = 3, LIBXSMM_ACC_PARAM_B = 4, LIBXSMM_ACC_PARAM_C = 5,
  LIBXSMM_ACC_PARAM_COUNT
};

enum libxsmm_acc_elem_type {
  LIBXSMM_ACC_ELEM_UNKNOWN = 0,
  LIBXSMM_ACC_ELEM_F32 = 1, LIBXSMM_ACC_ELEM_F64 = 3,
  LIBXSMM_ACC_ELEM_C32 = 5, LIBXSMM_ACC_ELEM_C64 = 7
};

template<typename T, bool Complex> struct libxsmm_acc_elem {
                                                    static const libxsmm_acc_elem_type type = LIBXSMM_ACC_ELEM_UNKNOWN;
                                                    static const char* name() { return "unknown"; } };
template<> struct libxsmm_acc_elem<float,false>   { static const libxsmm_acc_elem_type type = LIBXSMM_ACC_ELEM_F32;
                                                    static const char* name() { return "f32"; } };
template<> struct libxsmm_acc_elem<double,false>  { static const libxsmm_acc_elem_type type = LIBXSMM_ACC_ELEM_F64;
                                                    static const char* name() { return "f64"; } };
template<> struct libxsmm_acc_elem<float,true>    { static const libxsmm_acc_elem_type type = LIBXSMM_ACC_ELEM_C32;
                                                    static const char* name() { return "c32"; } };
template<> struct libxsmm_acc_elem<double,true>   { static const libxsmm_acc_elem_type type = LIBXSMM_ACC_ELEM_C64;
                                                    static const char* name() { return "c64"; } };

LIBXSMM_ACC_EXTERN void xsmm_acc_abort(const char* filename, int line_number, const char* message);

#endif // defined(__LIBXSMM) || (defined(__ACC) && defined(__ACC_MIC) && defined(__DBCSR_ACC) && defined(__LIBXSTREAM))
#endif // LIBXSMM_ACC_HPP
