AC_PREREQ(2.62)
AC_INIT(DIDA, 1.0, hmohamadi@bcgsc.ca, dida,
		http://www.bcgsc.ca/platform/bioinfo/software/dida)
AM_INIT_AUTOMAKE([subdir-objects])
AC_CONFIG_SRCDIR([dida.cpp])
AC_CONFIG_HEADER([config.h])

# Checks for programs.
AC_PROG_AWK
AC_PROG_CC
AC_PROG_CPP
AC_PROG_CXX
AC_PROG_INSTALL

# Checks for header files.
AC_CHECK_HEADERS([dlfcn.h fcntl.h float.h limits.h \
	stddef.h stdint.h stdlib.h sys/param.h])
AC_HEADER_STDBOOL
AC_HEADER_STDC

# Checks for library functions.
AC_CHECK_FUNCS([dup2 gethostname getopt_long getpagesize \
				memset strdup strerror strtoul])
AC_FUNC_FORK
AC_FUNC_MALLOC
AC_FUNC_MEMCMP
AC_FUNC_REALLOC
AC_FUNC_SETVBUF_REVERSED
AC_FUNC_VPRINTF

# Checks for typedefs, structures, and compiler characteristics.
AC_C_BIGENDIAN
AC_C_CONST
AC_C_INLINE
AC_CHECK_TYPES([ptrdiff_t])
AC_TYPE_MODE_T
AC_TYPE_PID_T
AC_TYPE_SIZE_T
AC_TYPE_SSIZE_T
AC_TYPE_INT64_T
AC_TYPE_UINT8_T
AC_TYPE_UINT16_T
AC_TYPE_UINT32_T
AC_TYPE_UINT64_T

# Check for the math library.
AC_CHECK_LIB([m], [sqrt])
AC_CHECK_FUNCS([pow sqrt])
AC_CHECK_FUNC(ceilf, [], AC_DEFINE(ceilf, [ceil],
			  [Define if the system does not provide ceilf]))

# Check for the dynamic linking library.
AC_CHECK_LIB([dl], [dlsym])

# Check for OpenMP.
AC_LANG_PUSH([C++])
AC_OPENMP
if test -z $OPENMP_CXXFLAGS; then
	AC_MSG_ERROR([DIDA must be compiled with a C++ compiler that supports OpenMP threading.])
fi
AC_LANG_POP([C++])

# Make sure we have pthreads
AC_CHECK_LIB([pthread], [pthread_create])

# MPI
AC_ARG_WITH(mpi, AS_HELP_STRING([--with-mpi=PATH],
	[specify prefix directory for the installed MPI parallel
	computing library]))
if test "$with_mpi" -a -d "$with_mpi"; then
	mpi_cppflags="-isystem$with_mpi/include"
	mpi_ldflags="-L$with_mpi/lib"
fi

AC_ARG_ENABLE(mpich, AS_HELP_STRING([--enable-mpich],
	[use MPICH (default is to use Open MPI)]))
AC_ARG_ENABLE(lammpi, AS_HELP_STRING([--enable-lammpi],
	[use LAM/MPI (default is to use Open MPI)]))

# Find the absolute path to the source.
my_abs_srcdir=$(cd $srcdir; pwd)

# Set compiler flags
AC_SUBST(CPPFLAGS,
		 "-I$my_abs_srcdir $mpi_cppflags $CPPFLAGS")
AC_SUBST(LDFLAGS, "$mpi_ldflags $OPENMP_CXXFLAGS $LDFLAGS")

# Check for the MPI parallel computing library.
libs="$LIBS"
AC_DEFINE(MPICH_SKIP_MPICXX, 1,
		  [Define to disable MPICH C++ bindings])
AC_DEFINE(OMPI_SKIP_MPICXX, 1,
		  [Define to disable OpenMPI C++ bindings])
AC_CHECK_HEADERS([mpi.h])
if test "$enable_mpich"; then
	AC_CHECK_LIB([mpl], [MPL_env2int])
	AC_CHECK_LIB([mpich], [MPI_Init])
	ac_cv_lib_mpi_MPI_Init=$ac_cv_lib_mpich_MPI_Init
elif test "$enable_lammpi"; then
	AC_CHECK_LIB([dl], [dlopen])
	AC_CHECK_LIB([lam], [lam_mutex_lock])
	AC_CHECK_LIB([mpi], [MPI_Init])
	AC_LANG_PUSH([C++])
	AC_CHECK_LIB([lammpi++], [main])
	AC_LANG_POP([C++])
else
	AC_CHECK_LIB([mpi], [MPI_Init])
fi

if test $ac_cv_header_mpi_h != yes -o $ac_cv_lib_mpi_MPI_Init != yes; then
	AC_MSG_ERROR([DIDA must be compiled with MPI support.])
fi

AC_SUBST(MPI_LIBS, "$LIBS")
LIBS="$libs"

AC_CONFIG_FILES([Makefile])
AC_OUTPUT
