#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.
# Autoconf requirements
# AC_INIT(package, version, bug-report-address)
# information on the package
# checks for programs
# checks for libraries
# checks for header files
# checks for types
# checks for structures
# checks for compiler characteristics
# checks for library functions
# checks for system services
# AC_CONFIG_FILES([file...])
# AC_OUTPUT

AC_PREREQ([2.69])
AC_INIT
AM_INIT_AUTOMAKE(dmrgcc, 1.0)
#AC_CONFIG_SRCDIR([src/dmrg/realdmrg/real_tensor_network.h])
#AC_CONFIG_HEADERS([config.h])

# Checks for programs.
AC_PROG_CXX
AC_PROG_CC
AC_PROG_MAKE_SET
AC_PROG_RANLIB

# Checks for libraries.
# AC_CHECK_LIB (library, function, [action-if-found], [action-if-not-found], [other-libraries])
AC_CHECK_LIB(m, sin, MATH_LIBS="-lm", AC_MSG_WARN(math lib not found))
AC_CHECK_LIB(pthread, pthread_create, THREAD_LIBS="-lpthread", AC_MSG_WARN(pthread lib not found))

# Checks for header files.
# stdlib.h, stdarg.h, string.h, folat.h
AC_HEADER_STDC
AC_CHECK_HEADERS(sys/time.h)

# Checks for typedefs, structures, and compiler characteristics.
AC_C_INLINE

# Check for LAPACK
#AC_ARG_ENABLE (feature, help-string, [action-if-given], [action-if-not-given])
AC_ARG_ENABLE(lapack, [AC_HELP_STRING([--enable-lapack], [enable LAPACK support])], 
              [HAVE_PKG_LAPACK=$enableval], [HAVE_PKG_LAPACK=no])
# AC_ARG_WITH (package, help-string, [action-if-given], [action-if-not-given])
USER_DEFINED_BLAS=no
AC_ARG_WITH(blas_libs, [AC_HELP_STRING([--with-blas-libs="LIBS"], [BLAS libraries @<:@"-lblas"@:>@])],
[BLAS_LIBS=${withval}
 AC_MSG_CHECKING(user-defined BLAS libraries)
 AC_MSG_RESULT([${BLAS_LIBS}])
 USER_DEFINED_BLAS=yes]
)

USER_DEFINED_LAPACK=no
AC_ARG_WITH(lapack_libs, [AC_HELP_STRING([--with-lapack-libs="LIBS"], [LAPACK libraries @<:@"-llapack"@:>@])],
[LAPACK_LIBS=${withval}
 AC_MSG_CHECKING(user-defined LAPACK libraries)
 AC_MSG_RESULT([${LAPACK_LIBS}])
 USER_DEFINED_LAPACK=yes]
)

AC_MSG_CHECKING(whether we are using LAPACK)
AC_MSG_RESULT([${HAVE_PKG_LAPACK}])

if test X${HAVE_PKG_LAPACK} = Xyes; then
    AC_DEFINE(HAVE_LAPACK,,[define if we are using LAPACK])
    AC_LANG_PUSH(Fortran 77)
    if test X${USER_DEFINED_BLAS} = Xno; then
        AC_CHECK_LIB(blas, dgemm, [HAVE_PKG_BLAS=$enableval;BLAS_LIBS="-lblas"], AC_MSG_ERROR(BLAS is not found))
    fi
    
    if test X${USER_DEFINED_LAPACK} = Xno; then
        AC_CHECK_LIB(lapack, dsyev, [HAVE_PKG_LAPACK=$enableval;LAPACK_LIBS="-llpack"], AC_MSG_ERROR(LAPACK is not found))
    fi
    AC_LANG_POP(Fortran 77)
fi

# Check for Intel MKL
AC_ARG_ENABLE(mkl, [AC_HELP_STRING([--enable-mkl], [enable MKL support])], 
              [HAVE_PKG_MKL=$enableval], [HAVE_PKG_MKL=no])

AC_MSG_CHECKING(whether we are using MKL)
AC_MSG_RESULT([${HAVE_PKG_MKL}])

AC_ARG_WITH(mkl_libdir, [AC_HELP_STRING([--with-mkl-libdir="DIR"], [MKL library directory @<:@MKLROOT/lib@:>@])],
[MKL_LIBDIR=${withval}
 AC_MSG_CHECKING(user-defined MKL libraries)
 AC_MSG_RESULT([${MKL_LIBDIR}])]
)

if test X${HAVE_PKG_MKL} = Xyes; then
    AC_CHECK_HEADERS([mkl.h], [], [AC_MSG_ERROR([mkl.h not found or not working])], [])
    AC_CHECK_LIB(mkl_intel_lp64 -lmkl_intel_thread -lmkl_rt -lmkl_core -liomp5 -lpthread -lm, cblas_dgemm, MKL_LIBS="-L${MKL_LIBDIR} -lmkl_intel_lp64 -lmkl_intel_thread \
                 -lmkl_rt -lmkl_core -liomp5 -lpthread", 
    AC_MSG_ERROR([Cannot find MKL]))
fi

# Checks for library functions.
CXXFLAGS="-g -O3 -W -Wall -msse2 -msse3 -msse4"

AC_ARG_ENABLE(gprof, [AC_HELP_STRING([--enable-gprof], [enable gprof mode])], [HAVE_GPROF=$enableval 
CXXFLAGS="${CXXFLAGS} -pg"], 
[HAVE_GPROF=no])
AC_MSG_CHECKING(whether we are using gprof)
AC_MSG_RESULT([${HAVE_GPROF}])

AC_ARG_ENABLE(gcov, [AC_HELP_STRING([--enable-gcov], [enable gprof mode])], [HAVE_GCOV=$enableval 
CXXFLAGS="${CXXFLAGS} -fprofile-arcs -ftest-coverage"], 
[HAVE_GCOV=no])
AC_MSG_CHECKING(whether we are using gcov)
AC_MSG_RESULT([${HAVE_GCOV}])

AC_MSG_CHECKING(environment for CXXFLAGS)
if test -n "${CXXFLAGS}"; then
    AC_MSG_RESULT([${CXXFLAGS}])
else
    CXXFLAGS=" -O3 -W -Wall -msse2 -msse3 -msse4"
    AC_MSG_RESULT([using defaults])
fi

INCLUDES="-I\$(top_srcdir)/src -I\$(top_srcdir)/include -I\$(top_srcdir)"
SRCDIR="\$(top_srcdir)/src"

AC_FUNC_ERROR_AT_LINE
AC_CHECK_FUNCS([gettimeofday pow sqrt])

# Substitude the variables into Makefile
AC_SUBST(CXXFLAGS)
AC_SUBST(INCLUDES)
AC_SUBST(BLAS_LIBS)
AC_SUBST(LAPACK_LIBS)
AC_SUBST(MKL_LIBS)
AC_SUBST(SRCDIR)



AC_CONFIG_FILES([Makefile
                 model/Makefile
                 model/hamiltonian/Makefile
                 model/space/Makefile
                 src/Makefile
                 src/dmrg/Makefile
                 src/dmrg/complexdmrg/Makefile
                 src/dmrg/realdmrg/Makefile
                 src/tensor/Makefile
                 src/util/Makefile
                 test/Makefile])
AC_OUTPUT
