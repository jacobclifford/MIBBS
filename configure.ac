
AC_INIT([MIBBS], [0.1])

AM_INIT_AUTOMAKE

AC_PROG_RANLIB
AC_CHECK_LIB(m,main)
AC_CHECK_LIB(gslcblas,main)
AC_CHECK_LIB(gsl,main)

AC_LANG(C++)

AC_PROG_CXX

AC_CONFIG_FILES([Makefile
                 cluster/Makefile
                 data_dir/Makefile
                 logodds/Makefile
                 miclasses/Makefile
                 rocenhancer/Makefile
                 rocproj/Makefile])


AC_OUTPUT
