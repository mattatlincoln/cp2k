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
#elif defined(LIBXSMM_OFFLOAD_TARGET)
# pragma offload_attribute(push,target(LIBXSMM_OFFLOAD_TARGET))
#endif
#include <iostream>
#include <cstdlib>
#if defined(__MKL) || defined(MKL_DIRECT_CALL_SEQ) || defined(MKL_DIRECT_CALL)
# include <mkl_service.h>
#endif
#if defined(__TBBMALLOC)
# include <tbb/scalable_allocator.h>
#endif
#if defined(__ACC) && defined(__ACC_MIC) && defined(__DBCSR_ACC) && defined(__LIBXSTREAM)
# include <libxstream_end.h>
#elif defined(LIBXSMM_OFFLOAD_TARGET)
# pragma offload_attribute(pop)
#endif


#if defined(__LIBXSMM)
namespace libxsmm_acc_private {
  /** Internal type-agnostic call-forwarding to CP2K/intel stack processing; this is called by dbcsr_mm_hostdrv_mp_xsmm_process_mm_stack_[s|d]. */
  template<typename T>
  void process_mm_stack(const libxsmm_acc_stackdesc_type* descriptor, /*const*/ int* params, const int* myvalue, const T* a, const T* b, T* c, int* efficient/*Boolean*/)
  {
    int result = LIBXSMM_ACC_ERROR_CONDITION;

    if (0 != descriptor && 0 != params && 0 != myvalue && 0 != a && 0 != b && 0 != c) {
      result = libsmm_acc_process( // TODO: fix const-correctness in libsmm_acc.h
        params, *myvalue, LIBXSMM_ACC_NPARAMS, libxsmm_acc_elem<T,false>::type, const_cast<T*>(a), const_cast<T*>(b), c,
        descriptor->max_m, descriptor->max_n, descriptor->max_k, descriptor->defined_mnk, 0/*stream*/);
      if (efficient) *efficient = 1;
    }

    switch (result) {
      case LIBXSMM_ACC_ERROR_CONDITION: LIBXSMM_ACC_ABORT("incorrect argument(s)"); break;
      default: if (LIBXSMM_ACC_ERROR_NONE != result) LIBXSMM_ACC_ABORT("unknown error");
    }
  }

  const char *const prefetch_env = getenv("CP2K_PREFETCH");
  const char *const reconf_env = getenv("CP2K_RECONFIGURE");
  const bool explicit_configure = (reconf_env && *reconf_env);
  const bool reconfigure = explicit_configure
    ? 0 != atoi(reconf_env)
#if defined(LIBXSMM_ACC_OFFLOAD_BUILD)
    : true;
#else
    : false;
#endif

} // namespace libxsmm_acc_private


int libxsmm_acc_prefetch = (libxsmm_acc_private::prefetch_env && *libxsmm_acc_private::prefetch_env)
  ? atoi(libxsmm_acc_private::prefetch_env)
  /* Select automatic prefetch strategy if no default prefetch was selected at build time of LIBXSMM. */
  : (0 <= LIBXSMM_PREFETCH ? LIBXSMM_PREFETCH : -1);


LIBXSMM_ACC_EXTERN void xsmm_acc_abort(const char* filename, int line_number, const char* message)
{
  if (filename && *filename) {
    std::cerr << filename << ':' << line_number << " - " << ((message && *message) ? message : "unknown error") << std::endl/*includes flush*/;
  }
  exit(-1);
}


#if defined(__RECONFIGURE)

# if defined(__GNUC__)
LIBXSMM_ACC_EXTERN LIBXSMM_ACC_ATTRIBUTE(weak)
# else
LIBXSMM_ACC_EXTERN
# endif
void LIBXSMM_ACC_FSYMBOL(__real_dbcsr_config_mp_dbcsr_set_conf_mm_stacksize)(const int*);
LIBXSMM_ACC_EXTERN void LIBXSMM_ACC_FSYMBOL(__wrap_dbcsr_config_mp_dbcsr_set_conf_mm_stacksize)(const int* value)
{
  int myvalue = value ? *value : -1;
  if (!libxsmm_acc_private::explicit_configure || libxsmm_acc_private::reconfigure) {
    const char *const env = getenv("CP2K_STACKSIZE");
    if (env && *env) {
      myvalue = atoi(env);
#if defined(LIBXSMM_ACC_STACKSIZE)
      if (0 >= myvalue) myvalue = LIBXSMM_ACC_STACKSIZE;
#endif
    }
  }
  if (0 < myvalue) {
    LIBXSMM_ACC_FSYMBOL(__real_dbcsr_config_mp_dbcsr_set_conf_mm_stacksize)(&myvalue);
  }
}


# if defined(__GNUC__)
LIBXSMM_ACC_EXTERN LIBXSMM_ACC_ATTRIBUTE(weak)
# else
LIBXSMM_ACC_EXTERN
# endif
void LIBXSMM_ACC_FSYMBOL(__real_dbcsr_config_mp_dbcsr_set_conf_mm_driver)(const int*);
LIBXSMM_ACC_EXTERN void LIBXSMM_ACC_FSYMBOL(__wrap_dbcsr_config_mp_dbcsr_set_conf_mm_driver)(const int* driver)
{
#if defined(MKL_ENABLE_AVX512_MIC) // only necessary for the version carrying KNL support for the first time
  mkl_enable_instructions(MKL_ENABLE_AVX512_MIC);
#endif
#if defined(__LIBXSMM)
  // eventually pre-generate dispatch tables for any static code
  libxsmm_init();
#endif

#if defined(LIBXSMM_ACC_MM_DRIVER)
  int mm_driver = LIBXSMM_ACC_MM_DRIVER;
#else
  LIBXSMM_ACC_ASSERT(driver);
  int mm_driver = *driver;
#endif

  // better leave "CP2K_DRIVER" environment variable undocumented
  // variable takes the internal literal/number representing MM driver
  const char *const driver_env = getenv("CP2K_DRIVER");
  if (driver_env && *driver_env) { // environment variable is present
    mm_driver = atoi(driver_env);
  }

  // make sure to reconfigure *after* the original configuration procedure ran
  LIBXSMM_ACC_FSYMBOL(__real_dbcsr_config_mp_dbcsr_set_conf_mm_driver)(&mm_driver);

  // define stacksize as part of the MM driver config (otherwise it remains undefined)
  LIBXSMM_ACC_FSYMBOL(__wrap_dbcsr_config_mp_dbcsr_set_conf_mm_stacksize)(NULL);

  const char *const multrec_env = getenv("CP2K_MULTREC");
  if (multrec_env && *multrec_env) { // environment variable is present
    extern int LIBXSMM_ACC_FSYMBOL(dbcsr_config_mp_multrec_limit);
    int multrec_limit = atoi(multrec_env);
    if (0 < multrec_limit) {
      LIBXSMM_ACC_FSYMBOL(dbcsr_config_mp_multrec_limit) = multrec_limit;
    }
#if defined(LIBXSMM_ACC_MULTREC_LIMIT)
    else {
      LIBXSMM_ACC_FSYMBOL(dbcsr_config_mp_multrec_limit) = LIBXSMM_ACC_MULTREC_LIMIT;
    }
#endif
  }

  const char *const dense_mult_env = getenv("CP2K_DENSE");
  int dense_mult = (dense_mult_env && *dense_mult_env) ? atoi(dense_mult_env) : -1;

  if (libxsmm_acc_private::reconfigure) {
#if defined(__TBBMALLOC)
    scalable_allocation_mode(TBBMALLOC_USE_HUGE_PAGES, 1);
#endif
    if (0 > dense_mult) dense_mult = 0;

#if defined(__ACC) && defined(__ACC_MIC) && defined(__DBCSR_ACC) && defined(__LIBXSTREAM)
# if defined(LIBXSMM_ACC_ACCDRV_POSTERIOR_STREAMS)
    extern int LIBXSMM_ACC_FSYMBOL(dbcsr_config_mp_accdrv_posterior_streams);
    LIBXSMM_ACC_FSYMBOL(dbcsr_config_mp_accdrv_posterior_streams) = LIBXSMM_ACC_ACCDRV_POSTERIOR_STREAMS;
# endif
# if defined(LIBXSMM_ACC_ACCDRV_POSTERIOR_BUFFERS)
    extern int LIBXSMM_ACC_FSYMBOL(dbcsr_config_mp_accdrv_posterior_buffers);
    LIBXSMM_ACC_FSYMBOL(dbcsr_config_mp_accdrv_posterior_buffers) = LIBXSMM_ACC_ACCDRV_POSTERIOR_BUFFERS;
# endif
# if defined(LIBXSMM_ACC_ACCDRV_PRIORITY_STREAMS)
    extern int LIBXSMM_ACC_FSYMBOL(dbcsr_config_mp_accdrv_priority_streams);
    LIBXSMM_ACC_FSYMBOL(dbcsr_config_mp_accdrv_priority_streams) = LIBXSMM_ACC_ACCDRV_PRIORITY_STREAMS;
# endif
# if defined(LIBXSMM_ACC_ACCDRV_PRIORITY_BUFFERS)
    extern int LIBXSMM_ACC_FSYMBOL(dbcsr_config_mp_accdrv_priority_buffers);
    LIBXSMM_ACC_FSYMBOL(dbcsr_config_mp_accdrv_priority_buffers) = LIBXSMM_ACC_ACCDRV_PRIORITY_BUFFERS;
# endif
# if defined(LIBXSMM_ACC_ACCDRV_MIN_NFLOPS_PERMM)
    extern int LIBXSMM_ACC_FSYMBOL(dbcsr_config_mp_accdrv_min_flop_process);
    LIBXSMM_ACC_FSYMBOL(dbcsr_config_mp_accdrv_min_flop_process) = LIBXSMM_ACC_ACCDRV_MIN_NFLOPS_PERMM;
# endif
#endif
  }

  if (0 <= dense_mult) {
    extern int LIBXSMM_ACC_FSYMBOL(dbcsr_mm_cannon_mp_dense_mult_default);
    LIBXSMM_ACC_FSYMBOL(dbcsr_mm_cannon_mp_dense_mult_default) = dense_mult;
  }
}


# if defined(__GNUC__)
LIBXSMM_ACC_EXTERN LIBXSMM_ACC_ATTRIBUTE(weak)
# else
LIBXSMM_ACC_EXTERN
# endif
void LIBXSMM_ACC_FSYMBOL(__real_dbcsr_config_mp_dbcsr_set_conf_comm_thread_load)(const int*);
LIBXSMM_ACC_EXTERN void LIBXSMM_ACC_FSYMBOL(__wrap_dbcsr_config_mp_dbcsr_set_conf_comm_thread_load)(const int* value)
{
  LIBXSMM_ACC_ASSERT(value);
  int myvalue = *value;
#if defined(LIBXSMM_ACC_COMM_THREAD_LOAD)
  if (libxsmm_acc_private::reconfigure) {
    myvalue = LIBXSMM_ACC_COMM_THREAD_LOAD;
  }
#endif
  LIBXSMM_ACC_FSYMBOL(__real_dbcsr_config_mp_dbcsr_set_conf_comm_thread_load)(&myvalue);
}


# if defined(__GNUC__)
LIBXSMM_ACC_EXTERN LIBXSMM_ACC_ATTRIBUTE(weak)
# else
LIBXSMM_ACC_EXTERN
# endif
void LIBXSMM_ACC_FSYMBOL(__real_dbcsr_config_mp_dbcsr_set_conf_use_mpi_filtering)(const int*);
LIBXSMM_ACC_EXTERN void LIBXSMM_ACC_FSYMBOL(__wrap_dbcsr_config_mp_dbcsr_set_conf_use_mpi_filtering)(const int* value)
{
  LIBXSMM_ACC_ASSERT(value);
  int myvalue = *value;
  if (!libxsmm_acc_private::explicit_configure || libxsmm_acc_private::reconfigure) {
    const char *const env = getenv("CP2K_FILTERING");
    if (env && *env) {
      myvalue = atoi(env);
    }
    else if (libxsmm_acc_private::reconfigure) {
      myvalue = 1;
    }
  }
  LIBXSMM_ACC_FSYMBOL(__real_dbcsr_config_mp_dbcsr_set_conf_use_mpi_filtering)(&myvalue);
}


# if defined(__GNUC__)
LIBXSMM_ACC_EXTERN LIBXSMM_ACC_ATTRIBUTE(weak)
# else
LIBXSMM_ACC_EXTERN
# endif
void LIBXSMM_ACC_FSYMBOL(__real_dbcsr_config_mp_dbcsr_set_conf_use_mpi_rma)(const int*);
LIBXSMM_ACC_EXTERN void LIBXSMM_ACC_FSYMBOL(__wrap_dbcsr_config_mp_dbcsr_set_conf_use_mpi_rma)(const int* value)
{
  LIBXSMM_ACC_ASSERT(value);
  int myvalue = *value;
#if defined(__MPI_VERSION) && (3 <= __MPI_VERSION)
  if (!libxsmm_acc_private::explicit_configure || libxsmm_acc_private::reconfigure) {
    const char *const env = getenv("CP2K_RMA");
    if (env && *env) {
      myvalue = atoi(env);
    }
    else if (libxsmm_acc_private::reconfigure) {
      myvalue = 0;
    }
  }
#endif
  LIBXSMM_ACC_FSYMBOL(__real_dbcsr_config_mp_dbcsr_set_conf_use_mpi_rma)(&myvalue);
}


LIBXSMM_ACC_EXTERN void LIBXSMM_ACC_FSYMBOL(dbcsr_mm_hostdrv_mp_xsmm_process_mm_stack_s)(
  const libxsmm_acc_stackdesc_type* descriptor, /*const*/ int* params, const int* myvalue,
  const float* a, const float* b, float* c, int* efficient/*Boolean*/);
LIBXSMM_ACC_EXTERN void LIBXSMM_ACC_FSYMBOL(xsmm_acc_process_mm_stack_s)(
  const libxsmm_acc_stackdesc_type* descriptor, /*const*/ int* params, const int* myvalue,
  const float* a, const float* b, float* c, int* efficient/*Boolean*/)
{
  if (!libxsmm_acc_private::explicit_configure || libxsmm_acc_private::reconfigure) {
    libxsmm_acc_private::process_mm_stack(
      descriptor, params, myvalue, a, b, c, efficient);
  }
  else { /* CP2K/trunk/master code path */
    LIBXSMM_ACC_FSYMBOL(dbcsr_mm_hostdrv_mp_xsmm_process_mm_stack_s)(
      descriptor, params, myvalue, a, b, c, efficient);
  }
}


LIBXSMM_ACC_EXTERN void LIBXSMM_ACC_FSYMBOL(dbcsr_mm_hostdrv_mp_xsmm_process_mm_stack_d)(
  const libxsmm_acc_stackdesc_type* descriptor, /*const*/ int* params, const int* myvalue,
  const double* a, const double* b, double* c, int* efficient/*Boolean*/);
LIBXSMM_ACC_EXTERN void LIBXSMM_ACC_FSYMBOL(xsmm_acc_process_mm_stack_d)(
  const libxsmm_acc_stackdesc_type* descriptor, /*const*/ int* params, const int* myvalue,
  const double* a, const double* b, double* c, int* efficient/*Boolean*/)
{
  if (!libxsmm_acc_private::explicit_configure || libxsmm_acc_private::reconfigure) {
    libxsmm_acc_private::process_mm_stack(
      descriptor, params, myvalue, a, b, c, efficient);
  }
  else { /* CP2K/trunk/master code path */
    LIBXSMM_ACC_FSYMBOL(dbcsr_mm_hostdrv_mp_xsmm_process_mm_stack_d)(
      descriptor, params, myvalue, a, b, c, efficient);
  }
}

#endif // defined(__RECONFIGURE)

#endif // defined(__LIBXSMM)
#endif // defined(__LIBXSMM) || (defined(__ACC) && defined(__ACC_MIC) && defined(__DBCSR_ACC) && defined(__LIBXSTREAM))
