#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_INIT(PHYPATTERNPROB, 0.0.1, mtholder@gmail.com)
AC_CONFIG_HEADERS([config.h])

AC_PREREQ(2.57)


################################################################################
#	According to (http://www.mail-archive.com/autoconf@gnu.org/msg14232.html)
#		this macro should be after AC_INIT but before AM_INIT_AUTOMAKE
################################################################################
AM_INIT_AUTOMAKE([1.6 foreign dist-zip tar-ustar filename-length-max=299])

AC_CONFIG_SRCDIR([src/phy-pat-class-prob.cpp])

# Checks for programs.
AC_PROG_LN_S
AC_PROG_CPP
AC_PROG_CC
AC_PROG_CXX

AC_PROG_INSTALL

NCL_INC_DIR="/usr/local/include"
NCL_LIB_DIR="/usr/local/lib"
NCL_BIN_DIR="/usr/local/bin"

AC_ARG_WITH(
	[ncl],
	AC_HELP_STRING(
		[--with-ncl=DIR],
		[Specify the root directory for the ncl version 2.1 or greater library (parent of the include/ncl and lib directories).  Omit this flag and configure will look for a system level NCL installation.]
		),
	[
	if ! test "$withval" = "yes" -o "$withval" = "no" ; then
		NCL_INC_DIR="$withval/include"
		NCL_LIB_DIR="$withval/lib/ncl"
		NCL_BIN_DIR="$withval/bin"
	fi
	])
if ! test -f "$NCL_INC_DIR/ncl/nxsdefs.h" ; then
    AC_MSG_ERROR([ncl version 2.1 or higher is a prerequisite for building. Use the --with-ncl argument to configure to specify NCL's prefix directory.])
fi

PYTBEAGLEHON_HOME="/usr/local"
AC_ARG_WITH(
	[pytbeaglehon],
	AC_HELP_STRING(
		[--with-pytbeaglehon=DIR],
		[Specify the root directory for the pytbeaglehon]
		),
	[
	if ! test "$withval" = "yes" -o "$withval" = "no" ; then
		PYTBEAGLEHON_HOME="$withval"
	fi
	])
if ! test -d "$PYTBEAGLEHON_HOME/include/pytbeaglehon" ; then
    AC_MSG_ERROR([pytbeaglehon is a prerequisite for building. Use the --with-pytbeaglehon argument to configure to specify pytbeaglehon's prefix directory.])
fi


BEAGLE_HOME="/usr/local"
AC_ARG_WITH(
	[beagle],
	AC_HELP_STRING(
		[--with-beagle=DIR],
		[Specify the root directory for the beagle]
		),
	[
	if ! test "$withval" = "yes" -o "$withval" = "no" ; then
		BEAGLE_HOME="$withval"
	fi
	])
if ! test -d "$BEAGLE_HOME/include/libhmsbeagle-1" ; then
    AC_MSG_ERROR([beagle is a prerequisite for building. Use the --with-beagle argument to configure to specify beagle's prefix directory.])
fi

LIBS="$LIBS -lncl -lfatbeagle -lhmsbeagle"
LDFLAGS="$LDFLAGS -L$NCL_LIB_DIR -L$BEAGLE_HOME/lib -L$PYTBEAGLEHON_HOME/lib"
CPPFLAGS="-I\$(top_srcdir) -I$NCL_INC_DIR -I$PYTBEAGLEHON_HOME/include -I$BEAGLE_HOME/include/libhmsbeagle-1"
AC_SUBST(CFLAGS)
AC_SUBST(CPPFLAGS)

AC_HEADER_STDC
AC_CHECK_HEADERS([float.h stdlib.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_CONST
AC_C_INLINE

# Checks for library functions.
AC_CHECK_FUNCS([memmove strchr strtol])
AC_CONFIG_FILES([				\
	Makefile 					\
	src/Makefile
	])
AC_OUTPUT
