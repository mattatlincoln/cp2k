# Arch file targeting Linux 64-bit using the Intel tool chain
#
PERL     = perl
CPP      = #cpp
AR       = xiar -r

# LIBINT: configure, build, and install
# Version 1.1.5 and 1.1.6 (tested)
#
# $ env \
#   AR=xiar CC=icc CXX=icpc \
#   ./configure --prefix=$HOME/libint \
#     --with-cc-optflags="-O1 -xHost" \
#     --with-cxx-optflags="-O1 -xHost" \
#     --with-libint-max-am=5 \
#     --with-libderiv-max-am1=4
# $ make
# $ make install
# $ make realclean
#
#LIBINTROOT = $(HOME)/libint

# LIBXC: configure, build, and install
# Version 2.2.2 (tested)
#
# $ autoreconf --force --install
# $ env \
#   AR=xiar \
#   FC=ifort F77=ifort F90=ifort FCFLAGS="-O2 -xHost" \
#   CC=icc CFLAGS="-O2 -xHost" \
#   ./configure --prefix=$HOME/libxc
# $ make
# $ make install
# $ make clean
#
#LIBXCROOT = $(HOME)/libxc

# LIBXSMM (https://github.com/hfp/libxsmm)
#
#LIBXSMMROOT = /path/to/libxsmm
ifneq (0,$(LIBXSMM))
  ifeq (,$(strip $(LIBXSMMROOT)))
    ifneq (,$(wildcard $(CP2KHOME)/../libxsmm/Makefile))
      LIBXSMMROOT = $(CP2KHOME)/../libxsmm
    else ifneq (,$(wildcard $(HOME)/libxsmm/Makefile))
      LIBXSMMROOT = $(HOME)/libxsmm
    else ifneq (,$(wildcard $(TOOLSRC)/toolchain/build/libxsmm*/Makefile))
      LIBXSMMROOT = $(dir $(wildcard $(TOOLSRC)/toolchain/build/libxsmm*/Makefile))
    endif
  endif
endif

# LIBXSTREAM: cp2k/tools/mic/libxstream or https://github.com/hfp/libxstream
# Please note that CP2K redistributes a tested version of LIBXSTREAM
# which is built automatically (if LIBXSTREAMROOT is set).
#
LIBXSTREAMROOT = $(TOOLSRC)/mic/libxstream

# Diagnostic message to be turned off
DIAG_DISABLE = 8290,8291,10010,10212,10411,11060

# DEFAULTS
#
BEEP ?= 1
JIT ?= 1
SSE ?= 0
AVX ?= 0
OPT ?= 2
MPI ?= 1
OMP ?= 1
ACC ?= 0
OCL ?= 0
SYM ?= 0
DBG ?= 0
IPO ?= 0
MKL_DIRECT ?= 0
MKL_STATIC ?= 1
RECONFIGURE ?= 1
MEMKIND ?= 1
OFFLOAD ?= 0
NESTED ?= 0
ELPA ?= 3

ifneq (0,$(OMP))
  MKL ?= 2
else ifneq (0,$(MPI))
  MKL ?= 1
else
  MKL ?= 2
endif

# consider more accurate -fp-model (C/C++: precise, Fortran: source)
FPFLAGS ?= -fp-model fast=2
FPFLAGS += -complex-limited-range
#FPFLAGS += -fast-transcendentals
#FPFLAGS += -fimf-domain-exclusion=1

# TBB malloc proxy is enabled if TBBROOT is set
TBBMALLOC ?= 1
# TBB runtime compatible with oldest supported GCC
TBBGCC_OLD ?= 1

ifeq (0,$(OFFLOAD))
  MIC ?= 0
else
  ACC = 1
endif

OPT1 = $(shell echo $$((1<$(OPT)?1:$(OPT))))
OPT2 = $(shell echo $$((2<$(OPT)?2:$(OPT))))

ifeq (0,$(ATTRIBS))
  DFLAGS += -DCP_DISABLE_ATTRIBS
endif

ifeq (1,$(shell echo $$((2 > $(DBG)))))
  ifeq (1,$(AVX))
    TARGET = -xAVX
  else ifeq (2,$(AVX))
    TARGET = -xCORE-AVX2
  else ifeq (3,$(AVX))
    ifeq (0,$(MIC))
      ifneq (file,$(origin MIC))
        TARGET = -xCORE-AVX512
      else
        TARGET = -xCOMMON-AVX512
      endif
    else
      TARGET = -xMIC-AVX512
    endif
  else ifneq (0,$(SSE))
    TARGET = -xSSE3
  else
    TARGET = -xHost
  endif
endif

# initial build flags
CPPFLAGS  = #
CXXFLAGS  = -std=c++0x
CFLAGS    = #
FCFLAGS   = -free -fpp
# heap-array threshold is set to 50% of the default OMP stacksize
FCFLAGS  += -heap-arrays 4096
LDFLAGS  += #-static-intel -static-libgcc -static-libstdc++
OPTFLAGS  = $(TARGET)

# workaround for issue "cannot find address of function"
#ATTRIBUTE = mic
#DIAG_DISABLE := $(DIAG_DISABLE),2571,3218

ifeq (0,$(DBG))
  OPTFLAGS  += -O$(OPT)
  DFLAGS    += -DNDEBUG
  CXXFLAGS  += -fno-alias -ansi-alias $(FPFLAGS)
  CFLAGS    += -fno-alias -ansi-alias $(FPFLAGS)
  FCFLAGS   += -align array64byte     $(FPFLAGS)
  #LDFLAGS   += $(NULL)
else
  OPTFLAGS += -O0
  CXXFLAGS := -debug $(CXXFLAGS)
  CFLAGS := -debug $(CFLAGS)
  FCFLAGS := -debug $(FCFLAGS)
  ifeq (2,$(DBG))
    FCFLAGS += -fpe0 # debugging NaNs
  endif
  SYM = $(DBG)
endif

ifneq (0,$(SYM))
  DFLAGS += -D__USE_CP2K_TRACE
  OPTFLAGS  += -traceback
  ifneq (1,$(SYM))
    CXXFLAGS := -ggdb3 $(CXXFLAGS)
    CFLAGS := -ggdb3 $(CFLAGS)
    FCFLAGS := -ggdb3 $(FCFLAGS)
  else
    CXXFLAGS := -g $(CXXFLAGS)
    CFLAGS := -g $(CFLAGS)
    FCFLAGS := -g $(FCFLAGS)
  endif
endif

ifneq (0,$(IPO))
  OPTFLAGS += -ipo-separate
else ifeq (0,$(IPO))
  LDFLAGS += -no-ipo
endif

SCALAPACK ?= 1
ifneq (0,$(MPI))
  CXX = mpiicpc
  CC  = mpiicc
  FC  = mpiifort
  LD  = mpiifort
  DFLAGS += -D__parallel
  ifneq (0,$(SCALAPACK))
    DFLAGS += -D__BLACS -D__SCALAPACK
    ifeq (0,$(shell echo $(SCALAPACK) | grep -q "^-*[0-9]\+$$"; echo $$?)) # number
      ifneq (-1,$(SCALAPACK))
        ifneq (1,$(SCALAPACK))
          DFLAGS += -D__SCALAPACK$(SCALAPACK)
        endif
      else
        SCALAPACKDIR = $(HOME)/scalapack-2.0.2
      endif
      ifeq (,$(SCALAPACKDIR))
        SCALAPACKDIR = $(MKLROOT)/lib/intel64
        SCALAPACKLIB = mkl_scalapack_lp64
      else
        ifeq (,$(SCALAPACKLIB))
          SCALAPACKLIB = scalapack
        endif
        LIBS += -L$(SCALAPACKDIR)
      endif
    else # NaN
      SCALAPACKDIR = $(SCALAPACK)
      SCALAPACKLIB = scalapack
      LIBS += -L$(SCALAPACKDIR)
    endif
  endif
  ifneq (1,$(MPI))
    DFLAGS += -D__MPI_VERSION=$(MPI)
  else # default MPI std. version
    DFLAGS += -D__MPI_VERSION=3
  endif
else
  CXX = icpc
  CC  = icc
  FC  = ifort
  LD  = ifort
endif

ifneq (0,$(OMP))
  FCFLAGS   += -threads
  LDFLAGS   += -threads
  OPTFLAGS  += -fopenmp -qoverride_limits
  ifneq (0,$(NESTED))
    DFLAGS += -D__NESTED_OPENMP
  endif
endif

ifneq (,$(LIBINTROOT))
  DFLAGS  += -D__LIBINT -D__LIBINT_MAX_AM=6 -D__LIBDERIV_MAX_AM1=5
  IFLAGS  += -I$(LIBINTROOT)/include
  LIBS    += $(LIBINTROOT)/lib/libderiv.a $(LIBINTROOT)/lib/libint.a
endif

ifneq (,$(LIBXCROOT))
  DFLAGS  += -D__LIBXC2
  IFLAGS  += -I$(LIBXCROOT)/include
  LIBS    += $(LIBXCROOT)/lib/libxcf90.a $(LIBXCROOT)/lib/libxc.a
endif

ifneq (,$(ELPAROOT))
  ifneq (0,$(ELPA))
    ifeq (1,$(ELPA))
      DFLAGS  += -D__ELPA
    else
      DFLAGS  += -D__ELPA$(ELPA)
    endif
    IFLAGS  += -I$(ELPAROOT)/include/elpa/modules
    LIBS    += $(ELPAROOT)/lib/libelpa.a
    # in case ELPA is built with OpenMP
    LIBS    += -Wl,--as-needed -liomp5 -Wl,--no-as-needed
  endif
endif

ifneq (0,$(TBBMALLOC))
  ifneq (,$(TBBROOT))
    GCC = $(notdir $(shell which gcc 2> /dev/null))
    ifneq (,$(GCC))
      GCC_VERSION_STRING = $(shell $(GCC) --version 2> /dev/null | head -n1 | sed "s/..* \([0-9][0-9]*\.[0-9][0-9]*\.*[0-9]*\)[ \S]*.*/\1/")
      GCC_VERSION_MAJOR = $(shell echo "$(GCC_VERSION_STRING)" | cut -d"." -f1)
      GCC_VERSION_MINOR = $(shell echo "$(GCC_VERSION_STRING)" | cut -d"." -f2)
      GCC_VERSION_PATCH = $(shell echo "$(GCC_VERSION_STRING)" | cut -d"." -f3)
      TBBLIBDIR = $(TBBROOT)/lib/intel64/gcc$(GCC_VERSION_MAJOR).$(GCC_VERSION_MINOR)
      TBBMALLOCLIB = $(wildcard $(TBBLIBDIR)/libtbbmalloc_proxy.so)
    endif
    ifeq (,$(TBBMALLOCLIB))
      ifneq (0,$(TBBGCC_OLD))
        TBBGCCDIR = $(shell ls -1 "$(TBBROOT)/lib/intel64" | tr "\n" " " | cut -d" " -f1)
      else
        TBBGCCDIR = $(shell ls -1 "$(TBBROOT)/lib/intel64" | tr "\n" " " | rev | cut -d" " -f2 | rev)
      endif
      TBBLIBDIR = $(TBBROOT)/lib/intel64/$(TBBGCCDIR)
      TBBMALLOCLIB = $(wildcard $(TBBLIBDIR)/libtbbmalloc_proxy.so)
    endif
    ifneq (,$(TBBMALLOCLIB))
      IFLAGS += -I$(TBBROOT)/include
      DFLAGS += -D__TBBMALLOC
      LIBS += $(TBBMALLOCLIB) $(TBBLIBDIR)/libtbbmalloc.so
      ifneq (1,$(TBBMALLOC)) # TBBMALLOC=2
        FCFLAGS += -heap-arrays
      endif
    endif
  endif
else ifneq (,$(TCMALLOCROOT))
  # configured using ./configure --enable-minimal --prefix=<TCMALLOCROOT>
  LIBS += $(TCMALLOCROOT)/lib/libtcmalloc_minimal.a
endif

ifneq (0,$(MEMKIND))
  ifneq (,$(MEMKINDROOT))
    #LIBS += -L$(MEMKINDROOT)/lib -lmemkind
    LIBS += $(MEMKINDROOT)/lib/libmemkind.a
  endif
endif

# Allow for LIBSMM to ease performance comparison...
LIBSMM ?= 0
ifneq (0,$(LIBSMM))
  LIBSMM_INSTALL := $(shell cd $(TOOLSRC)/toolchain; ./scripts/install_libsmm.sh)
  LIBSMM_LIB = $(TOOLSRC)/toolchain/install/lib/libsmm_dnn.a
endif
ifneq (,$(wildcard $(LIBSMM_LIB))) # LIBSMM successfully downloaded
  DFLAGS  += -D__HAS_smm_dnn
  LIBS    += $(LIBSMM_LIB)
endif

ifneq (,$(LIBXSMMROOT))
  LIBXSMM ?= 1
  ifneq (0,$(LIBXSMM))
    LIBXSMM_LIB = $(MAINLIBDIR)/$(ARCH)/$(ONEVERSION)/libxsmm/lib/libxsmm.a

    # substitute "big" xGEMM calls with LIBXSMM
    ifeq (1,$(shell echo $$((1 < $(LIBXSMM)))))
      LIBS += $(MAINLIBDIR)/$(ARCH)/$(ONEVERSION)/libxsmm/lib/libxsmmext.a
      LDFLAGS += -Wl,--wrap=sgemm_,--wrap=dgemm_
      ifeq (0,$(OMP))
        ifeq (1,$(MKL))
          LIBS += -liomp5
        else ifeq (0,$(MKL))
          LIBS += -liomp5
        endif
      endif
    endif

    ifeq (1,$(shell echo $$((0==$(JIT) || 1!=($(SSE)+1)))))
      LIBXSMM_MNK := "23, 6, 14 16 29, 14 32 29, 5 32 13 24 26, 9 32 22, 64, 78, 16 29 55, 32 29 55, 12, 4 5 7 9 13 25 26 28 32 45"
    endif
    LIBXSMM_ALIGNED_STORES := 0
    LIBXSMM_PREFETCH := 0
    ifneq (0,$(OMP))
      ifneq (0,$(NESTED))
        LIBXSMM_ALIGNED_STORES := 1
      endif
    endif
    ifneq (0,$(ACC))
      ifneq (0,$(OFFLOAD))
        LIBXSMM_ALIGNED_STORES := 1
        LIBXSMM_PREFETCH := 1
      endif
    endif
    ifneq (,$(filter %MIC-AVX512,$(TARGET)))
      LIBXSMM_PREFETCH := 1
    endif
    LIBXSMM_MPSS := 0
    ifneq (0,$(MIC))
      ifneq (3,$(AVX))
        LIBXSMM_MPSS := 1
      endif
    endif

$(LIBXSMM_LIB): .state
	$(info ================================================================================)
	$(info Automatically enabled LIBXSMM $(shell $(LIBXSMMROOT)/scripts/libxsmm_utilities.py 2> /dev/null))
	$(info LIBXSMMROOT=$(LIBXSMMROOT))
	$(info ================================================================================)
	@$(MAKE) --no-print-directory -f $(LIBXSMMROOT)/Makefile \
		INCDIR=$(MAINOBJDIR)/$(ARCH)/$(ONEVERSION)/libxsmm/include \
		BLDDIR=$(MAINOBJDIR)/$(ARCH)/$(ONEVERSION)/libxsmm/build \
		BINDIR=$(MAINOBJDIR)/$(ARCH)/$(ONEVERSION)/libxsmm/bin \
		OUTDIR=$(MAINLIBDIR)/$(ARCH)/$(ONEVERSION)/libxsmm/lib \
		MNK=$(LIBXSMM_MNK) M=$(LIBXSMM_M) N=$(LIBXSMM_N) K=$(LIBXSMM_K) PRECISION=2 \
		ALIGNED_STORES=$(LIBXSMM_ALIGNED_STORES) PREFETCH=$(LIBXSMM_PREFETCH) JIT=$(JIT) \
		PTHREAD=$(OMP) OPT=$(OPT) IPO=$(IPO) TARGET=$(TARGET) SSE=$(SSE) AVX=$(AVX) \
		SYM=$(SYM) DBG=$(DBG) MPSS=$(LIBXSMM_MPSS) OFFLOAD=$(OFFLOAD) MIC=$(MIC)
LIBXSMM_UPTODATE_CHECK := $(shell touch .state)
# translation unit (dummy) which allows to consider LIBXSMM_LIB as dep. in general
# below hack is not triggering minimal rebuild but "good enough" (relink)
$(SRCDIR)/base/base_uses.f90: $(LIBXSMM_LIB)
ifneq (,$(wildcard $(LIBXSMM_LIB)))
	@touch $(OBJDIR)/*.o
endif

    DFLAGS += -D__LIBXSMM=$(LIBXSMM)
    IFLAGS += -I$(MAINOBJDIR)/$(ARCH)/$(ONEVERSION)/libxsmm/include
    LIBS += $(LIBXSMM_LIB) $(MAINLIBDIR)/$(ARCH)/$(ONEVERSION)/libxsmm/lib/libxsmmf.a
  endif
endif

ifneq (0,$(ACC))
  DFLAGS += -D__ACC -D__DBCSR_ACC

  ifeq (0,$(OCL))
    DFLAGS += -D__ACC_MIC
    ifeq (0,$(OFFLOAD))
      OPTFLAGS += -qno-offload
      LDFLAGS += -qoffload-option,mic,ld,"--unresolved-symbols=ignore-all"
    else # also true if OFFLOAD is undefined
      #OPTFLAGS += -qoffload=mandatory
      # enable OpenMP for OFFLOAD regardless of wether OMP is enabled or not
      MIC_CXFLAGS += -qopenmp -qno-openmp -qoffload-option,mic,compiler,"-qopenmp"
      MIC_CCFLAGS += -qopenmp -qno-openmp -qoffload-option,mic,compiler,"-qopenmp"
      MIC_FCFLAGS += -qopenmp -qno-openmp -qoffload-option,mic,compiler,"-qopenmp"
      MIC_LDFLAGS += -qoffload-option,mic,ld,"--no-undefined"
      DIAG_DISABLE := $(DIAG_DISABLE),10121
      ifneq (,$(ATTRIBUTE))
        MIC_CXFLAGS += -qoffload-attribute-target=$(ATTRIBUTE)
        MIC_CCFLAGS += -qoffload-attribute-target=$(ATTRIBUTE)
        #MIC_FCFLAGS += -qoffload-attribute-target=$(ATTRIBUTE)
      endif
    endif
  else
    DFLAGS  += -D__OPENCL -D__USE_INTEL_CL
    LIBS    += -L/usr/lib64 -lOpenCL -lrt
  endif
else
  ifeq (0,$(OFFLOAD))
    OPTFLAGS += -qno-offload
    LDFLAGS += -qoffload-option,mic,ld,"--unresolved-symbols=ignore-all"
  endif

  # save some build time
  LIBXSTREAMROOT = $(NULL)
endif

ifneq (,$(LIBXSTREAMROOT))
  LIBXSTREAM_BUILD := $(shell $(MAKE) -f $(LIBXSTREAMROOT)/Makefile \
    BLDDIR=$(MAINOBJDIR)/$(ARCH)/$(ONEVERSION)/libxstream \
    OUTDIR=$(MAINLIBDIR)/$(ARCH)/$(ONEVERSION)/libxstream \
    SYM=$(SYM) DBG=$(DBG) IPO=$(IPO) \
    OFFLOAD=$(OFFLOAD) MIC=$(MIC) \
  >&2)

  DFLAGS  += -D__LIBXSTREAM
  IFLAGS  += -I$(LIBXSTREAMROOT)/include
  LIBS    += $(MAINLIBDIR)/$(ARCH)/$(ONEVERSION)/libxstream/libxstream.a
endif

ifeq (1,$(MKL_DIRECT))
  MKL_STATIC = 1
  DFLAGS += -DMKL_DIRECT_CALL_SEQ
endif

ifneq (1,$(MKL))
  ifneq (0,$(MKL)) # smp
    DFLAGS  += -D__MKL -D__FFTW3
    IFLAGS  +=-I$(MKLROOT)/include -I$(MKLROOT)/include/fftw
    ifeq (0,$(MKL_STATIC))
      LIBS += -L$(MKLROOT)/lib/intel64
      ifneq (0,$(MPI))
        MIC_LDFLAGS += -qoffload-option,mic,ld,"-L$(MKLROOT)/lib/mic -l$(SCALAPACKLIB) -lmkl_intel_lp64 -lmkl_core -lmkl_intel_thread -lmkl_blacs_intelmpi_lp64"
        LIBS += -l$(SCALAPACKLIB) -lmkl_intel_lp64 -lmkl_core -lmkl_intel_thread -lmkl_blacs_intelmpi_lp64
      else
        MIC_LDFLAGS += -qoffload-option,mic,ld,"-L$(MKLROOT)/lib/mic -lmkl_intel_lp64 -lmkl_core -lmkl_intel_thread"
        LIBS += -lmkl_intel_lp64 -lmkl_core -lmkl_intel_thread
      endif
    else # static
      ifneq (0,$(MPI))
        MIC_LDFLAGS += -qoffload-option,mic,ld," \
          --start-group \
            $(MKLROOT)/lib/mic/libmkl_scalapack_lp64.a \
            $(MKLROOT)/lib/mic/libmkl_intel_lp64.a \
            $(MKLROOT)/lib/mic/libmkl_core.a \
            $(MKLROOT)/lib/mic/libmkl_intel_thread.a \
            $(MKLROOT)/lib/mic/libmkl_blacs_intelmpi_lp64.a \
          --end-group"
        LIBS += \
          -Wl,--start-group \
            $(SCALAPACKDIR)/lib$(SCALAPACKLIB).a \
            $(MKLROOT)/lib/intel64/libmkl_intel_lp64.a \
            $(MKLROOT)/lib/intel64/libmkl_core.a \
            $(MKLROOT)/lib/intel64/libmkl_intel_thread.a \
            $(MKLROOT)/lib/intel64/libmkl_blacs_intelmpi_lp64.a \
          -Wl,--end-group
      else
        MIC_LDFLAGS += -qoffload-option,mic,ld," \
          --start-group \
            $(MKLROOT)/lib/mic/libmkl_intel_lp64.a \
            $(MKLROOT)/lib/mic/libmkl_core.a \
            $(MKLROOT)/lib/mic/libmkl_intel_thread.a \
          --end-group"
        LIBS += \
          -Wl,--start-group \
            $(MKLROOT)/lib/intel64/libmkl_intel_lp64.a \
            $(MKLROOT)/lib/intel64/libmkl_core.a \
            $(MKLROOT)/lib/intel64/libmkl_intel_thread.a \
          -Wl,--end-group
      endif
    endif
    ifeq (0,$(OMP))
      MIC_LDFLAGS += -qoffload-option,mic,ld,"-liomp5"
      LIBS += -liomp5
    endif
    MIC_LDFLAGS += -qoffload-option,mic,ld,"-lpthread -lm"
    LIBS += -Wl,--as-needed -lpthread -lm -Wl,--no-as-needed
  endif
else # sequential
  DFLAGS  += -D__MKL -D__FFTW3
  IFLAGS  +=-I$(MKLROOT)/include -I$(MKLROOT)/include/fftw
  ifeq (0,$(MKL_STATIC))
    LIBS += -L$(MKLROOT)/lib/intel64
    ifneq (0,$(MPI))
      MIC_LDFLAGS += -qoffload-option,mic,ld,"-L$(MKLROOT)/lib/mic -lmkl_scalapack_lp64 -lmkl_intel_lp64 -lmkl_core -lmkl_sequential -lmkl_blacs_intelmpi_lp64"
      LIBS += -l$(SCALAPACKLIB) -lmkl_intel_lp64 -lmkl_core -lmkl_sequential -lmkl_blacs_intelmpi_lp64
    else
      MIC_LDFLAGS += -qoffload-option,mic,ld,"-L$(MKLROOT)/lib/mic -lmkl_intel_lp64 -lmkl_core -lmkl_sequential"
      LIBS += -lmkl_intel_lp64 -lmkl_core -lmkl_sequential
    endif
  else # static
    ifneq (0,$(MPI))
      MIC_LDFLAGS += -qoffload-option,mic,ld," \
        --start-group \
          $(MKLROOT)/lib/mic/libmkl_scalapack_lp64.a \
          $(MKLROOT)/lib/mic/libmkl_intel_lp64.a \
          $(MKLROOT)/lib/mic/libmkl_core.a \
          $(MKLROOT)/lib/mic/libmkl_sequential.a \
          $(MKLROOT)/lib/mic/libmkl_blacs_intelmpi_lp64.a \
        --end-group"
      LIBS += \
        -Wl,--start-group \
          $(SCALAPACKDIR)/lib$(SCALAPACKLIB).a \
          $(MKLROOT)/lib/intel64/libmkl_intel_lp64.a \
          $(MKLROOT)/lib/intel64/libmkl_core.a \
          $(MKLROOT)/lib/intel64/libmkl_sequential.a \
          $(MKLROOT)/lib/intel64/libmkl_blacs_intelmpi_lp64.a \
        -Wl,--end-group
    else
      MIC_LDFLAGS += -qoffload-option,mic,ld," \
        --start-group \
          $(MKLROOT)/lib/mic/libmkl_intel_lp64.a \
          $(MKLROOT)/lib/mic/libmkl_core.a \
          $(MKLROOT)/lib/mic/libmkl_sequential.a \
        --end-group"
      MKL_LIBS = \
        -Wl,--start-group \
          $(MKLROOT)/lib/intel64/libmkl_intel_lp64.a \
          $(MKLROOT)/lib/intel64/libmkl_core.a \
          $(MKLROOT)/lib/intel64/libmkl_sequential.a \
        -Wl,--end-group
    endif
  endif
  MIC_LDFLAGS += -qoffload-option,mic,ld,"-lpthread -lm"
  LIBS += -Wl,--as-needed -lpthread -lm -Wl,--no-as-needed
endif

ifneq (,$(wildcard $(SRCDIR)/dbcsr/libsmm_acc/libxsmm_acc))
  ifeq (,$(LIBXSMMROOT))
    ifeq (,$(LIBXSTREAMROOT))
      RECONFIGURE = 0
    else ifeq (0,$(ACC))
      RECONFIGURE = 0
    endif
  endif
else # this is not the CP2K/intel branch 
  RECONFIGURE = 0
endif
ifneq (0,$(RECONFIGURE))
  DFLAGS  += -D__RECONFIGURE
  LDFLAGS += -Wl,--wrap=dbcsr_config_mp_dbcsr_set_conf_mm_driver_
  LDFLAGS += -Wl,--wrap=dbcsr_config_mp_dbcsr_set_conf_comm_thread_load_
  LDFLAGS += -Wl,--wrap=dbcsr_config_mp_dbcsr_set_conf_mm_stacksize_
  LDFLAGS += -Wl,--wrap=dbcsr_config_mp_dbcsr_set_conf_use_mpi_filtering_
  LDFLAGS += -Wl,--wrap=dbcsr_config_mp_dbcsr_set_conf_use_mpi_rma_
  DIAG_DISABLE := $(DIAG_DISABLE),11021
endif

DFLAGS  += -D__INTEL -D__HAS_ISO_C_BINDING
IFLAGS  += # general include paths

# Define __INTEL_COMPILER in case of external preprocessing because some source (pw/fft/fftw3_lib.F)
# toggles code using this symbol, but of course the cpp preprocessor is not defining this symbol.
CPPFLAGS  += #-C $(IFLAGS) $(DFLAGS) -D__INTEL_COMPILER -P -traditional

CXXFLAGS  += $(OPTFLAGS) -diag-disable $(DIAG_DISABLE) $(DFLAGS) $(IFLAGS)
CFLAGS    += $(OPTFLAGS) -diag-disable $(DIAG_DISABLE) $(DFLAGS) $(IFLAGS)
FCFLAGS   += $(OPTFLAGS) -diag-disable $(DIAG_DISABLE) $(DFLAGS) $(IFLAGS)
LDFLAGS   += $(OPTFLAGS) -diag-disable $(DIAG_DISABLE)
LDFLAGS_C += $(OPTFLAGS) -diag-disable $(DIAG_DISABLE) -nofor_main

LIBS += -Wl,--as-needed -lstdc++ -Wl,--no-as-needed
ifneq (0,$(ACC))
  ifneq (0,$(OFFLOAD))
    LIBS      += $(MIC_LDFLAGS)
    CXXFLAGS  += $(MIC_CXFLAGS)
    CFLAGS    += $(MIC_CCFLAGS)
    FCFLAGS   += $(MIC_FCFLAGS)
    #LDFLAGS   += $(MIC_LDFLAGS)
  endif
endif

ifeq (1,$(shell echo $$((1 <= $(BEEP)))))
mp2_optimize_ri_basis.o: mp2_optimize_ri_basis.F
	$(FC) -c $(FCFLAGS) -O0 $<
qs_dispersion_nonloc.o: qs_dispersion_nonloc.F
	$(FC) -c $(FCFLAGS) -O$(OPT1) $<
endif

# likely outdated or resolved
ifeq (1,$(shell echo $$((2 <= $(BEEP)))))
qs_vxc_atom.o: qs_vxc_atom.F
	$(FC) -c $(FCFLAGS) -O$(OPT1) $<
cp_fm_types.o: cp_fm_types.F
	$(FC) -c $(FCFLAGS) -O$(OPT1) $<
cube_utils.o: cube_utils.F
	$(FC) -c $(FCFLAGS) -O$(OPT1) $<
bibliography.o: bibliography.F
	$(FC) -c $(FCFLAGS) -O$(OPT2) $<
ifneq (0,$(OMP))
xc_tpss.o: xc_tpss.F
	$(FC) -c $(filter-out -openmp -qopenmp,$(FCFLAGS)) $<
realspace_grid_types.o: realspace_grid_types.F
	$(FC) -c $(filter-out -heap-arrays,$(FCFLAGS)) $<
matrix_exp.o: matrix_exp.F
	$(FC) -c $(filter-out -heap-arrays,$(FCFLAGS)) $<
cp_dbcsr_operations.o: cp_dbcsr_operations.F
	$(FC) -c $(filter-out -heap-arrays,$(FCFLAGS)) $<
dbcsr_util.o: dbcsr_util.F
	$(FC) -c $(filter-out -heap-arrays,$(FCFLAGS)) $<
qs_integrate_potential_product.o: qs_integrate_potential_product.F
	$(FC) -c $(FCFLAGS) -qno-openmp $<
dbcsr_work_operations.o: dbcsr_work_operations.F
	$(FC) -c $(FCFLAGS) -qno-openmp $<
endif
endif

