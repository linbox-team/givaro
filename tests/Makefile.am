# Copyright(c)'1994-2009 by The Givaro group
# This file is part of Givaro.
# Givaro is governed by the CeCILL-B license under French law
# and abiding by the rules of distribution of free software.
# see the COPYRIGHT file for more details.

SUBDIRS=
check: $(BASE_TESTS)

# -DSTD_RECINT_SIZE=n creates numbers of size 2^(2^n)
# -DLOOPS=m changes the number of checks per test
OPTFLAGS= -DSTD_RECINT_SIZE=8
OPTLINKS=

AM_CXXFLAGS = $(GIVARO_CXXFLAGS) $(OPTFLAGS)
AM_CPPFLAGS = -I$(top_builddir) -I$(top_builddir)/src/kernel/system -I$(top_builddir)/src/kernel/memory -I$(top_builddir)/src/kernel/field -I$(top_builddir)/src/kernel/ring -I$(top_builddir)/src/kernel/integer -I$(top_builddir)/src/kernel/rational -I$(top_builddir)/src/kernel -I$(top_builddir)/src/library/poly1 -I$(top_builddir)/src/library/vector -I$(top_builddir)/src/library/matrix -I$(top_builddir)/src/kernel/bstruct -I$(top_builddir)/src/library/tools -I$(top_builddir)/src/kernel/recint -I$(top_srcdir)/src/kernel

LDADD = $(OPTLINKS) ../src/libgivaro.la $(GMP_LIBS)
AM_LDFLAGS=-static $(LDFLAGS)

PERFPUBLISHERFILE = tests-report.xml

EXTRA_DIST =            \
		test-fieldarith.h

BASIC_TESTS =               \
		test-regression     \
		test-extension      \
		test-random         \
		test-ifactor        \
		test-ffarith        \
		test-qfield         \
		test-ringarith      \
		test-mod            \
		test-modsqroot      \
		test-brillhart      \
		test-trunc          \
		test-matrix         \
		test-poly           \
		test-midmul         \
		test-crt            \
		test-geom           \
		test-integer        \
		test-conversion     \
		test-rational-class \
		test-ratrecon       \
		test-modularmulprecomp

# *_mg tests accept -DMG_DEFAULT=MG_ACTIVE to enable Montgomery representation
RECINT_TESTS=                    \
    test-rmint_arith_addmul_mg   \
    test-rmint_arith_addmul      \
    test-ruint_arith_addmul      \
    test-rmint_arith_add_mg      \
    test-rmint_arith_add         \
    test-ruint_arith_add         \
    test-rint_arith_add          \
    test-rmint_arith_sub_mg      \
    test-rmint_arith_sub         \
    test-ruint_arith_sub         \
    test-rint_arith_sub          \
    test-rmint_arith_mul_mg      \
    test-rmint_arith_mul         \
    test-ruint_arith_mul         \
    test-rint_arith_mul          \
    test-rmint_arith_mod_mg      \
    test-rmint_arith_mod         \
    test-ruint_arith_div         \
    test-rint_arith_mod          \
    test-rmint_inv_div_mg        \
    test-rmint_inv_div           \
    test-rmint_constructors_mg   \
    test-rmint_constructors      \
    test-rmint_neg_mg            \
    test-rmint_neg               \
    test-rmint_comparisons_mg    \
    test-rmint_comparisons       \
    test-rmint_complex_mg        \
    test-rmint_complex           \
    test-recint_rand_mg          \
    test-recint_rand             \
    test-recint_exp_mg           \
    test-recint_exp              \
    test-recint_extra_mg         \
    test-recint_extra            \
    test-recint_convert_mg       \
    test-recint_convert          \
    test-recint_cast             \
    test-rmint_montgomery        \
    test-ruint_arazi             \
    test-ruint_comparisons       \
    test-ruint_complex           \
    test-ruint_constructors      \
    test-ruint_fiddling          \
    test-ruint_lmul              \
    test-ruint_operators         \
    test-ruint_shift             \
    test-ruint_square

EXTRA_PROGRAMS = $(BASIC_TESTS) $(RECINT_TESTS)
CLEANFILES = $(EXTRA_PROGRAMS) $(PERFPUBLISHERFILE)
TESTS = $(EXTRA_PROGRAMS)

#=====
# Sources

test_regression_SOURCES   = test-regression.C
test_extension_SOURCES   = test-extension.C
test_random_SOURCES   = test-random.C
test_ifactor_SOURCES  = test-ifactor.C
test_ffarith_SOURCES  = test-ffarith.C test-fieldarith.h
test_qfield_SOURCES  = test-qfield.C test-fieldarith.h
test_ringarith_SOURCES= test-ringarith.C
test_modularmulprecomp_SOURCES= test-modularmulprecomp.C
test_mod_SOURCES      = test-mod.C
test_modsqroot_SOURCES= test-modsqroot.C
test_brillhart_SOURCES= test-brillhart.C
test_trunc_SOURCES    = test-trunc.C
test_poly_SOURCES    = test-poly.C
test_midmul_SOURCES    = test-midmul.C
test_matrix_SOURCES    = test-matrix.C
test_crt_SOURCES      = test-crt.C
test_geom_SOURCES     = test-geom.C
test_integer_SOURCES  = test-integer.C
test_rational_class_SOURCES = test-rational-class.C
test_ratrecon_SOURCES = test-ratrecon.C
test_conversion_SOURCES = test-conversion.C

test_rmint_arith_addmul_mg_SOURCES   = test-rmint_arith.C
test_rmint_arith_addmul_SOURCES      = test-rmint_arith.C
test_ruint_arith_addmul_SOURCES      = test-ruint_arith.C
test_rmint_arith_add_mg_SOURCES      = test-rmint_arith.C
test_rmint_arith_add_SOURCES         = test-rmint_arith.C
test_ruint_arith_add_SOURCES         = test-ruint_arith.C
test_rint_arith_add_SOURCES          = test-rint_arith.C
test_rmint_arith_sub_mg_SOURCES      = test-rmint_arith.C
test_rmint_arith_sub_SOURCES         = test-rmint_arith.C
test_ruint_arith_sub_SOURCES         = test-ruint_arith.C
test_rint_arith_sub_SOURCES          = test-rint_arith.C
test_rmint_arith_mul_mg_SOURCES      = test-rmint_arith.C
test_rmint_arith_mul_SOURCES         = test-rmint_arith.C
test_ruint_arith_mul_SOURCES         = test-ruint_arith.C
test_rint_arith_mul_SOURCES          = test-rint_arith.C
test_rmint_arith_mod_mg_SOURCES      = test-rmint_arith.C
test_rmint_arith_mod_SOURCES         = test-rmint_arith.C
test_ruint_arith_div_SOURCES         = test-ruint_arith.C
test_rint_arith_mod_SOURCES          = test-rint_arith.C

test_rmint_inv_div_mg_SOURCES        = test-rmint_inv_div.C
test_rmint_inv_div_SOURCES           = test-rmint_inv_div.C
test_rmint_constructors_mg_SOURCES   = test-rmint_constructors.C
test_rmint_constructors_SOURCES      = test-rmint_constructors.C
test_rmint_neg_mg_SOURCES            = test-rmint_neg.C
test_rmint_neg_SOURCES               = test-rmint_neg.C
test_rmint_comparisons_mg_SOURCES    = test-rmint_comparisons.C
test_rmint_comparisons_SOURCES       = test-rmint_comparisons.C
test_rmint_complex_mg_SOURCES        = test-rmint_complex.C
test_rmint_complex_SOURCES           = test-rmint_complex.C
test_recint_rand_mg_SOURCES          = test-recint_rand.C
test_recint_rand_SOURCES             = test-recint_rand.C
test_recint_exp_mg_SOURCES           = test-recint_exp.C
test_recint_exp_SOURCES              = test-recint_exp.C
test_recint_extra_mg_SOURCES         = test-recint_extra.C
test_recint_extra_SOURCES            = test-recint_extra.C
test_recint_convert_mg_SOURCES       = test-recint_convert.C
test_recint_convert_SOURCES          = test-recint_convert.C

test_recint_cast_SOURCES         = test-recint_cast.C
test_rmint_montgomery_SOURCES    = test-rmint_montgomery.C
test_ruint_arazi_SOURCES         = test-ruint_arazi.C
test_ruint_comparisons_SOURCES   = test-ruint_comparisons.C
test_ruint_complex_SOURCES       = test-ruint_complex.C
test_ruint_constructors_SOURCES  = test-ruint_constructors.C
test_ruint_fiddling_SOURCES      = test-ruint_fiddling.C
test_ruint_lmul_SOURCES          = test-ruint_lmul.C
test_ruint_operators_SOURCES     = test-ruint_operators.C
test_ruint_shift_SOURCES         = test-ruint_shift.C
test_ruint_square_SOURCES        = test-ruint_square.C

#=====
# Additionnal flags

test_rmint_arith_addmul_mg_CXXFLAGS  = $(AM_CXXFLAGS) -DMG_DEFAULT=MG_ACTIVE -DNOT_IN_PLACE -DRI_OP=addmul -DGMP_OP=mpz_addmul -DGMP_OPUI=mpz_addmul_ui
test_rmint_arith_addmul_CXXFLAGS     = $(AM_CXXFLAGS) -DNOT_IN_PLACE -DRI_OP=addmul -DGMP_OP=mpz_addmul -DGMP_OPUI=mpz_addmul_ui
test_ruint_arith_addmul_CXXFLAGS     = $(AM_CXXFLAGS) -DNOT_IN_PLACE -DRI_OP=addmul -DGMP_OP=mpz_addmul -DGMP_OPUI=mpz_addmul_ui
test_rmint_arith_add_mg_CXXFLAGS     = $(AM_CXXFLAGS) -DMG_DEFAULT=MG_ACTIVE -DRI_OP=add -DGMP_OP=mpz_add -DGMP_OPUI=mpz_add_ui
test_rmint_arith_add_CXXFLAGS        = $(AM_CXXFLAGS) -DRI_OP=add -DGMP_OP=mpz_add -DGMP_OPUI=mpz_add_ui
test_ruint_arith_add_CXXFLAGS        = $(AM_CXXFLAGS) -DRI_OP=add -DGMP_OP=mpz_add -DGMP_OPUI=mpz_add_ui
test_rint_arith_add_CXXFLAGS         = $(AM_CXXFLAGS) -DRI_OP=add -DGMP_OP=mpz_add -DGMP_OPUI=mpz_add_ui
test_rmint_arith_sub_mg_CXXFLAGS     = $(AM_CXXFLAGS) -DMG_DEFAULT=MG_ACTIVE -DRI_OP=sub -DGMP_OP=mpz_sub -DGMP_OPUI=mpz_sub_ui
test_rmint_arith_sub_CXXFLAGS        = $(AM_CXXFLAGS) -DRI_OP=sub -DGMP_OP=mpz_sub -DGMP_OPUI=mpz_sub_ui
test_ruint_arith_sub_CXXFLAGS        = $(AM_CXXFLAGS) -DRI_OP=sub -DGMP_OP=mpz_sub -DGMP_OPUI=mpz_sub_ui
test_rint_arith_sub_CXXFLAGS         = $(AM_CXXFLAGS) -DRI_OP=sub -DGMP_OP=mpz_sub -DGMP_OPUI=mpz_sub_ui
test_rmint_arith_mul_mg_CXXFLAGS     = $(AM_CXXFLAGS) -DMG_DEFAULT=MG_ACTIVE -DRI_OP=mul -DGMP_OP=mpz_mul -DGMP_OPUI=mpz_mul_ui
test_rmint_arith_mul_CXXFLAGS        = $(AM_CXXFLAGS) -DRI_OP=mul -DGMP_OP=mpz_mul -DGMP_OPUI=mpz_mul_ui
test_ruint_arith_mul_CXXFLAGS        = $(AM_CXXFLAGS) -DRI_OP=mul -DGMP_OP=mpz_mul -DGMP_OPUI=mpz_mul_ui
test_rint_arith_mul_CXXFLAGS         = $(AM_CXXFLAGS) -DRI_OP=mul -DGMP_OP=mpz_mul -DGMP_OPUI=mpz_mul_ui
test_rmint_arith_mod_mg_CXXFLAGS     = $(AM_CXXFLAGS) -DMG_DEFAULT=MG_ACTIVE -DRI_OP=mod -DGMP_OP=mpz_mod -DGMP_OPUI=mpz_mod_ui
test_rmint_arith_mod_CXXFLAGS        = $(AM_CXXFLAGS) -DRI_OP=mod -DGMP_OP=mpz_mod -DGMP_OPUI=mpz_mod_ui
test_ruint_arith_div_CXXFLAGS        = $(AM_CXXFLAGS) -DNOT_IN_PLACE -DRI_OP=div_q -DGMP_OP=mpz_fdiv_q -DGMP_OPUI=mpz_fdiv_q_ui
test_rint_arith_mod_CXXFLAGS         = $(AM_CXXFLAGS) -DARITH_MOD -DNOT_IN_PLACE -DRI_OP=div_r -DGMP_OP=mpz_tdiv_r

test_rmint_inv_div_mg_CXXFLAGS        = $(AM_CXXFLAGS) -DMG_DEFAULT=MG_ACTIVE
test_rmint_constructors_mg_CXXFLAGS   = $(AM_CXXFLAGS) -DMG_DEFAULT=MG_ACTIVE
test_rmint_neg_mg_CXXFLAGS            = $(AM_CXXFLAGS) -DMG_DEFAULT=MG_ACTIVE
test_rmint_comparisons_mg_CXXFLAGS    = $(AM_CXXFLAGS) -DMG_DEFAULT=MG_ACTIVE
test_rmint_complex_mg_CXXFLAGS        = $(AM_CXXFLAGS) -DMG_DEFAULT=MG_ACTIVE
test_recint_rand_mg_CXXFLAGS          = $(AM_CXXFLAGS) -DMG_DEFAULT=MG_ACTIVE
test_recint_exp_mg_CXXFLAGS           = $(AM_CXXFLAGS) -DMG_DEFAULT=MG_ACTIVE
test_recint_extra_mg_CXXFLAGS         = $(AM_CXXFLAGS) -DMG_DEFAULT=MG_ACTIVE
test_recint_convert_mg_CXXFLAGS       = $(AM_CXXFLAGS) -DMG_DEFAULT=MG_ACTIVE

#=====
# Perfpublisher script interaction - AB 2014/11/17
perfpublisher:
	+./perfpublisher.sh "$(PERFPUBLISHERFILE)" "$(EXTRA_PROGRAMS)" "$(CXX)"
	mv $(PERFPUBLISHERFILE) ..
	make clean
	mv ../$(PERFPUBLISHERFILE) .


#=====
# for compilation of new tests
GIVARO_BIN=@bindir@
define comp_new_test
	$(AM_V_CXX)$(CXXCOMPILE) $(OPTFLAGS) -c -o $@.$(OBJEXT) $< `$(GIVARO_BIN)/givaro-config --cflags`
	$(AM_V_CXXLD)$(CXXLINK) $@.$(OBJEXT) `$(GIVARO_BIN)/givaro-config --libs` $(LDADD)
endef

%:%.C
	$(comp_new_test)

%:%.cpp
	$(comp_new_test)
