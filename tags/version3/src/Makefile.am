INCLUDES=-I$(top_srcdir)

SUBDIRS=kernel library 

lib_LTLIBRARIES= libgivaro.la

libgivaro_la_SOURCES = dummy.C

libgivaro_la_LIBADD= \
	kernel/gmp++/libgivgmpxx.la kernel/bstruct/libgivbstruct.la kernel/integer/libgivinteger.la kernel/memory/libgivmemory.la kernel/rational/libgivrational.la kernel/system/libgivsystem.la kernel/zpz/libgivzpz.la library/tools/libgivtools.la library/poly1/libgivpoly1.la

