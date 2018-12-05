#!/bin/bash
# This file is part of the Givaro library.
# It is distributed under the terms of the CeCILL-B licence 
# (see COPYING)
# Created by AB - 2014/12/03
# Modified by AC - 2016/04/08
# Modified by CP - 2016/06/22

# Some influential environment variables:
#	CXX			C++ compiler command
#	CXXFLAGS	C++ compiler flags

# Note: This script is intended to be launched
# by the Jenkins web interface whenever it needs
# to compile the project.
# But should be stored in /<slave_jenkins_path>/makers/
# It is launched from the svn:trunk root directory.

SOURCE_DIRECTORY=$( cd "$( dirname "$0" )" && pwd )


#=============================#
# Change only these variables #
#=============================#
ARCH=`pwd | awk -F/ '{print $(NF-2)}'`
CXX=`pwd | awk -F/ '{print $(NF)}'`
#SSE=`pwd | awk -F/ '{print $NF}'`

# Job givaro with SSE option flag
# by default sse is enabled
#if [ "$SSE" == "withoutSSE" ]; then
#  GIVARO_SSEFLAG="--disable-sse --disable-sse2 --disable-sse3 --disable-ssse3 --disable-sse4.1 --disable-sse4.2 --disable-avx --disable-avx2 --disable-fma --disable-fma4"
#fi

JENKINS_DIR=${SOURCE_DIRECTORY%%/workspace/*}
LOCAL_DIR="$JENKINS_DIR"/local

# Where to install givaro binaries
# Keep default for local installation.
PREFIX_INSTALL="$LOCAL_DIR/$CXX"

# Add path to compilers (if needed)
export PATH="$PATH":"/usr/local/bin":"$LOCAL_DIR/$CXX/bin"
# Add specific locations (if needed)
export LD_LIBRARY_PATH="$LD_LIBRARY_PATH":"/usr/local/lib":"$LOCAL_DIR/$CXX/lib"
echo "LD_LIBRARY_PATH = ${LD_LIBRARY_PATH}"

# Where is GMP installed (compiled with cxx interface)
# Keep empty if in usual folders (i.e. /usr or /usr/local)
GMP_PATH="$LOCAL_DIR/$CXX"

# /!\ Warning /!\ This could be an issue if you changed
# the local installation directory
rm -rf "$PREFIX_INSTALL"/bin/givaro* "$PREFIX_INSTALL"/include/givaro* "$PREFIX_INSTALL"/include/gmp++ "$PREFIX_INSTALL"/include/recint "$PREFIX_INSTALL"/lib/libgivaro*


#================#
# Setup Variables#
#================#

if [ "$CXX" == "icpc" ]; then
     distribution=`uname -m`
     if [ "$distribution" == "i686" ]; then 	
	source /usr/local/bin/compilervars.sh ia32
     else
	source /usr/local/bin/compilervars.sh intel64
     fi
fi

# Particular case for Fedora23: g++=g++-5.3
if [[ "$ARCH" == "linbox-fedora-amd64" &&  "$CXX" == "g++-6" ]]; then
    CXX="g++"
fi
CC=`echo $CXX | sed  's/icpc/icc/;s/clang++/clang/;s/++/cc/'`

#==================================#
# Automated installation and tests #
#==================================#

echo "|=== JENKINS AUTOMATED SCRIPT ===| ./autogen.sh CXX=$CXX CC=$CC CXXFLAGS=$CXXFLAGS --prefix=$PREFIX_INSTALL --with-gmp=$GMP_PATH"
./autogen.sh CXX=$CXX CC=$CC CXXFLAGS=$CXXFLAGS --prefix="$PREFIX_INSTALL" --with-gmp="$GMP_PATH"
V="$?"; if test "x$V" != "x0"; then exit "$V"; fi

echo "|=== JENKINS AUTOMATED SCRIPT ===| make install"
make install
V="$?"; if test "x$V" != "x0"; then exit "$V"; fi

echo "|=== JENKINS AUTOMATED SCRIPT ===| make perfpublisher"
make perfpublisher

echo "|=== JENKINS AUTOMATED SCRIPT ===| make examples"
make examples
V="$?"; if test "x$V" != "x0"; then exit "$V"; fi
(cd examples && make clean)
