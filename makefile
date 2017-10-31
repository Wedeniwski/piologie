#################################
##
## Piologie V 1.3.2
## multi-precision arithmetic
## Makefile
##
## (c) 1996-2001 HiPiLib
##               www.hipilib.de
##
## Sebastian Wedeniwski
## 10/21/2001
##

#
# include the configuration file
#

include config

#
# delete a file:
#
# DELETE = del       # Windows / DOS / OS/2
DELETE = rm -f     # UNIX / MKS-Tools / GNU-Tools

#
# copy a file:
#
# COPY = copy       # Windows / DOS / OS/2
COPY = cp         # UNIX / MKS-Tools / GNU-Tools


PIOLOGIE_OBJS   = digit.$(OBJ) natural.$(OBJ) integer.$(OBJ) rational.$(OBJ)
SUPPLEMENT_OBJS = nmbrthry.$(OBJ) pi.$(OBJ) ispprime.$(OBJ)
ALL_OBJS        = $(PIOLOGIE_OBJS) $(SUPPLEMENT_OBJS)

all: lib

objs: $(ALL_OBJS)

lib: $(PIOLOGIE_OBJS)
	$(DELETE) piologie.$(LIB_SUFFIX)
	$(MAKE_KEYLIB)

keylib: digit.$(SUFFIX) digit.h test.$(SUFFIX)
	$(COMPILE_ONLY) $(CFLAGS) $(CFLAGS_KEY) digit.$(SUFFIX)
	$(MAKE) objs
	$(DELETE) piologie.$(LIB_SUFFIX)
	$(MAKE_KEYLIB)
	$(COMPILE_ONLY) $(CFLAGS) $(CFLAGS_KEY) test.$(SUFFIX)
	$(LINK_TEST)

piologie.$(LIB_SUFFIX): $(ALL_OBJS)
	$(DELETE) $@
	$(MAKE_LIB)

digit.$(SUFFIX):
	$(COPY) digit.cpp $@

digit.$(OBJ): digit.$(SUFFIX) digit.h
	$(COMPILE_ONLY) $(CFLAGS) digit.$(SUFFIX)

natural.$(SUFFIX):
	$(COPY) natural.cpp $@

natural.$(OBJ): natural.$(SUFFIX) natural.h digit.h
	$(COMPILE_ONLY) $(CFLAGS) natural.$(SUFFIX)

integer.$(SUFFIX):
	$(COPY) integer.cpp $@

integer.$(OBJ): integer.$(SUFFIX) integer.h natural.h digit.h
	$(COMPILE_ONLY) $(CFLAGS) integer.$(SUFFIX)

rational.$(SUFFIX):
	$(COPY) rational.cpp $@

rational.$(OBJ): rational.$(SUFFIX) rational.h integer.h natural.h digit.h
	$(COMPILE_ONLY) $(CFLAGS) rational.$(SUFFIX)

nmbrthry.$(SUFFIX):
	$(COPY) nmbrthry.cpp $@

nmbrthry.$(OBJ): nmbrthry.$(SUFFIX) modulo.h natural.h digit.h
	$(COMPILE_ONLY) $(CFLAGS) nmbrthry.$(SUFFIX)

pi.$(SUFFIX):
	$(COPY) pi.cpp $@

pi.$(OBJ): pi.$(SUFFIX) pi.h nmbrthry.h integer.h natural.h digit.h
	$(COMPILE_ONLY) $(CFLAGS) pi.$(SUFFIX)


test.$(SUFFIX):
	$(COPY) test.cpp $@

test: test.$(SUFFIX) piologie.$(LIB_SUFFIX)
	$(COMPILE_ONLY) $(CFLAGS) $@.$(SUFFIX)
	$(LINK)

cfrac.$(SUFFIX):
	$(COPY) cfrac.cpp $@

cfrac: cfrac.$(SUFFIX) piologie.$(LIB_SUFFIX)
	$(COMPILE_ONLY) $(CFLAGS) $@.$(SUFFIX)
	$(LINK)

check.$(SUFFIX):
	$(COPY) check.cpp $@

check: check.$(SUFFIX) piologie.$(LIB_SUFFIX)
	$(COMPILE_ONLY) $@.$(SUFFIX)
	$(LINK)

constant.$(SUFFIX):
	$(COPY) constant.cpp $@

constant: constant.$(SUFFIX) piologie.$(LIB_SUFFIX)
	$(COMPILE_ONLY) $(CFLAGS) $@.$(SUFFIX)
	$(LINK)

zeta.$(SUFFIX):
	$(COPY) zeta.cpp $@

zeta: zeta.$(SUFFIX) piologie.$(LIB_SUFFIX)
	$(COMPILE_ONLY) $(CFLAGS) $@.$(SUFFIX)
	$(LINK)
	$(MAKE) makesh

makesh.$(SUFFIX):
	$(COPY) makesh.cpp $@

makesh: makesh.$(SUFFIX)
	$(COMPILE_ONLY) $(CFLAGS) $@.$(SUFFIX)
	$(LINK)


dvi: piologie.dvi manual.dvi
	cd doc
	cd src
	nuweb -o piologie.w
	latex piologie.tex
	nuweb -o piologie.w
	makeindex piologie.idx
	latex piologie.tex
	latex manual.tex
	latex manual.tex
	cd ..
	cd ..

ps: piologie.ps manual.ps dvi
	cd doc
	cd src
	dvips piologie.dvi
	dvips manual.dvi
	cd ..
	cd ..

clean:
	$(DELETE) *.$(OBJ)


#
# Supported Compiler
#

ap: clean
	$(COPY) config.ap config; make

borland: clean
	$(COPY) config.i386_bc config; make

dec: clean
	$(COPY) config.dec config; make

edg: clean
	$(COPY) config.edg config; make

gnu: clean
	$(COPY) config.gnu config; make

gnu28: clean
	$(COPY) config.gnu28 config; make

gnu_mips4: clean
	$(COPY) config.gnu_mips4 config; make

gnu_sparc: clean
	$(COPY) config.gnu_sparc config; make

gnu28_sparc: clean
	$(COPY) config.gnu28_sparc config; make

hp: clean
	$(COPY) config.hp config; make

hpa: clean
	$(COPY) config.hpa config; make

i386_ibm: clean
	$(COPY) config.i386_ibm config
	make

kai: clean
	$(COPY) config.kai config; make

os390: clean
	$(COPY) config.os390 config
	make

ppc_ibm: clean
	$(COPY) config.ppc_ibm config
	make

sgi: clean
	$(COPY) config.sgi config; make

sgi_8000: clean
	$(COPY) config.sgi_8000 config; make

sun: clean
	$(COPY) config.sun config; make

visual: clean
	$(COPY) config.i386_vc config; make

watcom: clean
	$(COPY) config.i386_wc config; make
