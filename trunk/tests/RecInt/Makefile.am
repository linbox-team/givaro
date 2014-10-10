# Copyright(c)'1994-2009 by The Givaro group
# This file is part of Givaro.
# Givaro is governed by the CeCILL-B license under French law
# and abiding by the rules of distribution of free software.
# see the COPYRIGHT file for more details.
# --- Makefile.am stolen and modified by Alexis Breust

all: $(EXTRA_PROGRAMS)

# -DSTD_RECINT_SIZE=n creates numbers of size 2^(2^n)
# -DLOOPS=m changes the number of checks per test
OPTFLAGS = -DSTD_RECINT_SIZE=8
OPTLINKS=

AM_CXXFLAGS = @TESTS_CFLAGS@ $(OPTFLAGS) $(GMP_CFLAGS) -I$(top_srcdir) -I$(top_srcdir)/src/kernel/system -I$(top_srcdir)/src/kernel/recint

LDADD = $(OPTLINKS) -L${top_srcdir}/src -L${top_srcdir}/src/.libs -lgivaro $(GMP_LIBS) $(LDFLAGS)
AM_LDFLAGS = -static

#=====
# Programs
# MG_PROGS accept -DMG_DEFAULT=MG_ACTIVE
# to enable Montgomery representations by default

EXTRA_PROGRAMS=             \
    rmint_arith_addmul_mg   \
    rmint_arith_addmul      \
    ruint_arith_addmul      \
    rmint_arith_add_mg      \
    rmint_arith_add         \
    ruint_arith_add         \
    rmint_arith_sub_mg      \
    rmint_arith_sub         \
    ruint_arith_sub         \
    rmint_arith_mul_mg      \
    rmint_arith_mul         \
    ruint_arith_mul         \
    rmint_arith_mod_mg      \
    rmint_arith_mod         \
    ruint_arith_div         \
    rmint_inv_div_mg        \
    rmint_inv_div           \
    rmint_constructors_mg   \
    rmint_constructors      \
    rmint_neg_mg            \
    rmint_neg               \
    rmint_comparisons_mg    \
    rmint_comparisons       \
    rmint_complex_mg        \
    rmint_complex           \
    recint_rand_mg          \
    recint_rand             \
    recint_exp_mg           \
    recint_exp              \
    recint_extra_mg         \
    recint_extra            \
    recint_convert_mg       \
    recint_convert          \
    recint_cast             \
    rmint_montgomery        \
    ruint_arazi             \
    ruint_comparisons       \
    ruint_complex           \
    ruint_constructors      \
    ruint_fiddling          \
    ruint_lmul              \
    ruint_operators         \
    ruint_shift             \
    ruint_square

CLEANFILES=$(EXTRA_PROGRAMS)
TESTS=$(EXTRA_PROGRAMS)

#=====
# Sources

rmint_arith_addmul_mg_SOURCES   = rmint_arith.C
rmint_arith_addmul_SOURCES      = rmint_arith.C
ruint_arith_addmul_SOURCES      = ruint_arith.C
rmint_arith_add_mg_SOURCES      = rmint_arith.C
rmint_arith_add_SOURCES         = rmint_arith.C
ruint_arith_add_SOURCES         = ruint_arith.C
rmint_arith_sub_mg_SOURCES      = rmint_arith.C
rmint_arith_sub_SOURCES         = rmint_arith.C
ruint_arith_sub_SOURCES         = ruint_arith.C
rmint_arith_mul_mg_SOURCES      = rmint_arith.C
rmint_arith_mul_SOURCES         = rmint_arith.C
ruint_arith_mul_SOURCES         = ruint_arith.C
rmint_arith_mod_mg_SOURCES      = rmint_arith.C
rmint_arith_mod_SOURCES         = rmint_arith.C
ruint_arith_div_SOURCES         = ruint_arith.C
    
rmint_inv_div_mg_SOURCES        = rmint_inv_div.C
rmint_inv_div_SOURCES           = rmint_inv_div.C
rmint_constructors_mg_SOURCES   = rmint_constructors.C
rmint_constructors_SOURCES      = rmint_constructors.C
rmint_neg_mg_SOURCES            = rmint_neg.C
rmint_neg_SOURCES               = rmint_neg.C
rmint_comparisons_mg_SOURCES    = rmint_comparisons.C
rmint_comparisons_SOURCES       = rmint_comparisons.C
rmint_complex_mg_SOURCES        = rmint_complex.C
rmint_complex_SOURCES           = rmint_complex.C
recint_rand_mg_SOURCES          = recint_rand.C
recint_rand_SOURCES             = recint_rand.C
recint_exp_mg_SOURCES           = recint_exp.C
recint_exp_SOURCES              = recint_exp.C
recint_extra_mg_SOURCES         = recint_extra.C
recint_extra_SOURCES            = recint_extra.C
recint_convert_mg_SOURCES       = recint_convert.C
recint_convert_SOURCES          = recint_convert.C

recint_cast_SOURCES         = recint_cast.C
rmint_montgomery_SOURCES    = rmint_montgomery.C
ruint_arazi_SOURCES         = ruint_arazi.C
ruint_comparisons_SOURCES   = ruint_comparisons.C
ruint_complex_SOURCES       = ruint_complex.C
ruint_constructors_SOURCES  = ruint_constructors.C
ruint_fiddling_SOURCES      = ruint_fiddling.C
ruint_lmul_SOURCES          = ruint_lmul.C
ruint_operators_SOURCES     = ruint_operators.C
ruint_shift_SOURCES         = ruint_shift.C
ruint_square_SOURCES        = ruint_square.C

#=====
# Additionnal flags

rmint_arith_addmul_mg_CXXFLAGS  = $(AM_CXXFLAGS) -DMG_DEFAULT=MG_ACTIVE -DNOT_IN_PLACE -DRI_OP=addmul -DGMP_OP=mpz_addmul -DGMP_OPUI=mpz_addmul_ui
rmint_arith_addmul_CXXFLAGS     = $(AM_CXXFLAGS) -DNOT_IN_PLACE -DRI_OP=addmul -DGMP_OP=mpz_addmul -DGMP_OPUI=mpz_addmul_ui
ruint_arith_addmul_CXXFLAGS     = $(AM_CXXFLAGS) -DNOT_IN_PLACE -DRI_OP=addmul -DGMP_OP=mpz_addmul -DGMP_OPUI=mpz_addmul_ui
rmint_arith_add_mg_CXXFLAGS     = $(AM_CXXFLAGS) -DMG_DEFAULT=MG_ACTIVE -DRI_OP=add -DGMP_OP=mpz_add -DGMP_OPUI=mpz_add_ui
rmint_arith_add_CXXFLAGS        = $(AM_CXXFLAGS) -DRI_OP=add -DGMP_OP=mpz_add -DGMP_OPUI=mpz_add_ui
ruint_arith_add_CXXFLAGS        = $(AM_CXXFLAGS) -DRI_OP=add -DGMP_OP=mpz_add -DGMP_OPUI=mpz_add_ui
rmint_arith_sub_mg_CXXFLAGS     = $(AM_CXXFLAGS) -DMG_DEFAULT=MG_ACTIVE -DRI_OP=sub -DGMP_OP=mpz_sub -DGMP_OPUI=mpz_sub_ui
rmint_arith_sub_CXXFLAGS        = $(AM_CXXFLAGS) -DRI_OP=sub -DGMP_OP=mpz_sub -DGMP_OPUI=mpz_sub_ui
ruint_arith_sub_CXXFLAGS        = $(AM_CXXFLAGS) -DRI_OP=sub -DGMP_OP=mpz_sub -DGMP_OPUI=mpz_sub_ui
rmint_arith_mul_mg_CXXFLAGS     = $(AM_CXXFLAGS) -DMG_DEFAULT=MG_ACTIVE -DRI_OP=mul -DGMP_OP=mpz_mul -DGMP_OPUI=mpz_mul_ui
rmint_arith_mul_CXXFLAGS        = $(AM_CXXFLAGS) -DRI_OP=mul -DGMP_OP=mpz_mul -DGMP_OPUI=mpz_mul_ui
ruint_arith_mul_CXXFLAGS        = $(AM_CXXFLAGS) -DRI_OP=mul -DGMP_OP=mpz_mul -DGMP_OPUI=mpz_mul_ui
rmint_arith_mod_mg_CXXFLAGS     = $(AM_CXXFLAGS) -DMG_DEFAULT=MG_ACTIVE -DRI_OP=mod -DGMP_OP=mpz_mod -DGMP_OPUI=mpz_mod_ui
rmint_arith_mod_CXXFLAGS        = $(AM_CXXFLAGS) -DRI_OP=mod -DGMP_OP=mpz_mod -DGMP_OPUI=mpz_mod_ui
ruint_arith_div_CXXFLAGS        = $(AM_CXXFLAGS) -DNOT_IN_PLACE -DRI_OP=div_q -DGMP_OP=mpz_fdiv_q -DGMP_OPUI=mpz_fdiv_q_ui

rmint_inv_div_mg_CXXFLAGS        = $(AM_CXXFLAGS) -DMG_DEFAULT=MG_ACTIVE
rmint_constructors_mg_CXXFLAGS   = $(AM_CXXFLAGS) -DMG_DEFAULT=MG_ACTIVE
rmint_neg_mg_CXXFLAGS            = $(AM_CXXFLAGS) -DMG_DEFAULT=MG_ACTIVE
rmint_comparisons_mg_CXXFLAGS    = $(AM_CXXFLAGS) -DMG_DEFAULT=MG_ACTIVE
rmint_complex_mg_CXXFLAGS        = $(AM_CXXFLAGS) -DMG_DEFAULT=MG_ACTIVE
recint_rand_mg_CXXFLAGS          = $(AM_CXXFLAGS) -DMG_DEFAULT=MG_ACTIVE
recint_exp_mg_CXXFLAGS           = $(AM_CXXFLAGS) -DMG_DEFAULT=MG_ACTIVE
recint_extra_mg_CXXFLAGS         = $(AM_CXXFLAGS) -DMG_DEFAULT=MG_ACTIVE
recint_convert_mg_CXXFLAGS       = $(AM_CXXFLAGS) -DMG_DEFAULT=MG_ACTIVE

