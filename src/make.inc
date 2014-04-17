ESTER_ROOT=$(CURDIR)/..

USE_PGPLOT=0

CC=g++
CCFLAGS=-O3 -fPIC
FC=gfortran
FCFLAGS=-O3 -fPIC -fno-align-commons -fdefault-real-8 -fdefault-double-8 -ffixed-line-length-132
LD=g++
LDFLAGS=

BLAS_LIB=-lblas
LAPACK_LIB=-llapack

LIBS=$(LAPACK_LIB) $(BLAS_LIB)  -L/usr/local/Cellar/gfortran/4.8.2/gfortran/lib -lgfortran

ifeq ($(USE_PGPLOT),1)
	LIBS:=-lcpgplot $(LIBS)
endif

BLAS_H=$(ESTER_ROOT)/src/include/cblas.h
LAPACK_H=$(ESTER_ROOT)/src/include/lapack.h
