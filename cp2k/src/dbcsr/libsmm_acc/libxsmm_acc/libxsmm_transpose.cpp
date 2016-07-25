/*****************************************************************************
 *  CP2K: A general program to perform molecular dynamics simulations        *
 *  Copyright (C) 2000 - 2016  CP2K developers group                         *
 *****************************************************************************/

//! **************************************************************************
//!> \author Hans Pabst (Intel Corp.)
//! **************************************************************************

#if defined(__LIBXSMM) || (defined(__ACC) && defined(__ACC_MIC) && defined(__DBCSR_ACC) && defined(__LIBXSTREAM))
#include "libxsmm_acc.h"

#if defined(__ACC) && defined(__ACC_MIC) && defined(__DBCSR_ACC) && defined(__LIBXSTREAM)
# include <libxstream_begin.h>
#endif
#include <cstddef>
#include <cassert>
#if defined(LIBXSMM_ACC_OPENMP)
# include <omp.h>
#endif
#if defined(__ACC) && defined(__ACC_MIC) && defined(__DBCSR_ACC) && defined(__LIBXSTREAM)
# include <libxstream_end.h>
#endif

#define LIBXSMM_ACC_MAX_MATRIX_SIZE LIBXSMM_ACC_MAX(LIBXSMM_ACC_MAX( \
  LIBXSMM_ACC_MAX_M * LIBXSMM_ACC_MAX_K,  \
  LIBXSMM_ACC_MAX_K * LIBXSMM_ACC_MAX_N), \
  LIBXSMM_ACC_MAX_M * LIBXSMM_ACC_MAX_N)

#define LIBXSMM_ACC_MAX_MNK LIBXSMM_ACC_MAX(LIBXSMM_ACC_MAX( \
  LIBXSMM_ACC_MAX_M, LIBXSMM_ACC_MAX_N), \
  LIBXSMM_ACC_MAX_K)

#if defined(LIBXSMM_ACC_MKLTRANS) && defined(__MKL)
LIBXSMM_ACC_EXTERN void MKL_Simatcopy(const char, const char, size_t, size_t, float, float*, size_t, size_t);
LIBXSMM_ACC_EXTERN void MKL_Dimatcopy(const char, const char, size_t, size_t, double, double*, size_t, size_t);
#endif


namespace libxsmm_transpose_private {

LIBXSMM_ACC_RETARGETABLE inline/*IPO*/ void mkl_imatcopy(size_t m, size_t n, float* matrix)
{
#if defined(LIBXSMM_ACC_MKLTRANS) && defined(__MKL)
  MKL_Simatcopy('R', 'T', m, n, 1.f, matrix, n, m);
#endif
}


LIBXSMM_ACC_RETARGETABLE inline/*IPO*/ void mkl_imatcopy(size_t m, size_t n, double* matrix)
{
#if defined(LIBXSMM_ACC_MKLTRANS) && defined(__MKL)
  MKL_Dimatcopy('R', 'T', m, n, 1.0, matrix, n, m);
#endif
}


template<typename T, typename U>
LIBXSMM_ACC_RETARGETABLE void kernel(const U* stack, const U* pstacksize, const U* pm, const U* pn, T* matrix)
{
  const U stacksize = *pstacksize, m = *pm, n = *pn;

#if defined(LIBXSMM_ACC_OPENMP)
# pragma omp parallel for LIBXSMM_ACC_SCHEDULE
#endif
  for (U s = 0; s < stacksize; ++s) {
    T *const mat = matrix + stack[s];

#if defined(LIBXSMM_ACC_MKLTRANS) && defined(__MKL)
    mkl_imatcopy(static_cast<size_t>(m), static_cast<size_t>(n), mat);
#else
    T buffer[LIBXSMM_ACC_MAX_MATRIX_SIZE];
    T *const tmp = LIBXSMM_ACC_ALIGN_LDST(buffer);

    for (U i = 0; i < m; ++i) {
      LIBXSMM_ACC_PRAGMA_LOOP_COUNT(1, LIBXSMM_ACC_MAX_N, 23)
      for (U j = 0; j < n; ++j) {
        tmp[j*m+i] = mat[i*n+j];
      }
    }

    for (U i = 0; i < m; ++i) {
      LIBXSMM_ACC_PRAGMA_LOOP_COUNT(1, LIBXSMM_ACC_MAX_N, 23)
      for (U j = 0; j < n; ++j) {
        mat[i*n+j] = tmp[i*n+j];
      }
    }
#endif
  }
}


template<typename T, bool Complex, typename U>
int transpose(const U* stack, U offset, U nblocks, U m, U n, void* data, void* stream)
{
  LIBXSMM_ACC_CHECK_CONDITION(
    0 != stack && 0 != data && 0 <= offset && 0 <= nblocks &&
    LIBXSMM_ACC_MAX_MATRIX_SIZE >= (m * n) && 0 <= m && 0 <= n);

  if (1 < m || 1 < n) {
#if defined(__ACC) && defined(__ACC_MIC) && defined(__DBCSR_ACC) && defined(__LIBXSTREAM)
    if (0 != stream) {
      const size_t stacksize = nblocks;
      libxstream_argument* signature = 0;
      LIBXSMM_ACC_CHECK_CALL_ASSERT(libxstream_fn_signature(&signature));
      LIBXSMM_ACC_CHECK_CALL_ASSERT(libxstream_fn_input(signature, 0, stack + offset, libxstream_map_to<U>::type(), 1, &stacksize));
      LIBXSMM_ACC_CHECK_CALL_ASSERT(libxstream_fn_input(signature, 1, &nblocks, libxstream_map_to<U>::type(), 0, 0));
      LIBXSMM_ACC_CHECK_CALL_ASSERT(libxstream_fn_input(signature, 2,       &m, libxstream_map_to<U>::type(), 0, 0));
      LIBXSMM_ACC_CHECK_CALL_ASSERT(libxstream_fn_input(signature, 3,       &n, libxstream_map_to<U>::type(), 0, 0));
      LIBXSMM_ACC_CHECK_CALL_ASSERT(libxstream_fn_inout(signature, 4,     data, libxstream_map_to<T>::type(), 1, 0/*unknown*/));
      const libxstream_function libxsmm_acc_transpose_function = reinterpret_cast<libxstream_function>(kernel<T,U>);
      LIBXSMM_ACC_CHECK_CALL_ASSERT(libxstream_fn_call(libxsmm_acc_transpose_function, signature, static_cast<libxstream_stream*>(stream), LIBXSTREAM_CALL_DEFAULT));
    }
    else
#endif
    {
      kernel(stack + offset, &nblocks, &m, &n, static_cast<T*>(data))
    }
  }

  return LIBXSMM_ACC_ERROR_NONE;
}

} // namespace libxsmm_transpose_private

#if defined(__ACC) && defined(__ACC_MIC) && defined(__DBCSR_ACC) && defined(__LIBXSTREAM)
// workaround for issue "cannot find address of function" (use unoptimized build or apply mic attribute globally)
const libxstream_function libxsmm_acc_transpose_function = reinterpret_cast<libxstream_function>(libxsmm_transpose_private::kernel<double,int>);
#endif


LIBXSMM_ACC_EXTERN int libsmm_acc_transpose(void* trs_stack, int offset, int nblks, void* buffer, int datatype, int m, int n, void* stream)
{
  int result = LIBXSMM_ACC_ERROR_NONE;

#if defined(LIBXSMM_ACC_PRETRANSPOSE)
  const int *const stack = static_cast<const int*>(trs_stack);
  switch(static_cast<libxsmm_acc_elem_type>(datatype)) {
    case LIBXSMM_ACC_ELEM_F32: {
      result = libxsmm_transpose_private::transpose<float,false>(stack, offset, nblks, m, n, buffer, stream);
    } break;
    case LIBXSMM_ACC_ELEM_F64: {
      result = libxsmm_transpose_private::transpose<double,false>(stack, offset, nblks, m, n, buffer, stream);
    } break;
    case LIBXSMM_ACC_ELEM_C32: {
      LIBXSMM_ACC_ASSERT(false/*TODO: not implemented yet*/);
      result = LIBXSMM_ACC_ERROR_CONDITION;
    } break;
    case LIBXSMM_ACC_ELEM_C64: {
      LIBXSMM_ACC_ASSERT(false/*TODO: not implemented yet*/);
      result = LIBXSMM_ACC_ERROR_CONDITION;
    } break;
    default:
      result = LIBXSMM_ACC_ERROR_CONDITION;
  }
#endif

  return result;
}

#endif // defined(__LIBXSMM) || (defined(__ACC) && defined(__ACC_MIC) && defined(__DBCSR_ACC) && defined(__LIBXSTREAM))
