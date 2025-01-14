ACLOCAL_AMFLAGS = -I m4

# Why configure does not do that ??
#AM_FCFLAGS=-I/usr/include #-I/home/abeelen/Bureau/Poker/Poker/Healpix/Healpix_2.15a/src/f90/mod

AM_FCFLAGS = -I$(HEALPIX_INC)

bin_PROGRAMS = poker_mbb poker_count_tasks poker_give_map_b

pkglib_LTLIBRARIES = libpoker.la

# healpix is not compiled with -fPIC so I can not add it here (where it should be)
libpoker_la_LIBADD= $(LIBS) $(FLIBS) 

libpoker_la_SOURCES = \
poker_types.F90 \
poker_input.F90 \
poker_tools.F90

LDADD=libpoker.la -lhealpix $(LIBS) $(FLIBS) -L$(HEALPIX_LIB)

poker_mbb_SOURCES = \
	poker_mbb.F90

poker_count_tasks_SOURCES = \
	poker_count_tasks.F90 

poker_give_map_b_SOURCES = \
	poker_give_map_b.F90
