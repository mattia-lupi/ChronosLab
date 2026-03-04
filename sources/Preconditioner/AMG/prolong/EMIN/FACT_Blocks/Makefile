CXX           = g++

CXXFLAGS      = -O2 -std=c++11 -fopenmp $(INCLUDE)
PCXXFLAGS     = -O1 -pg -std=c++11 $(INCLUDE)

CPLUS_INCLUDE_PATH = /usr/include/x86_64-linux-gnu/c++/4.8

CXX           = CPLUS_INCLUDE_PATH=$(CPLUS_INCLUDE_PATH) g++

INCLUDE       = -I. -I./include -I./ICHOL

DFLAGS        = -O0 -std=c++11 -g -fopenmp -Wall $(INCLUDE)

DEST	      = .

EXTHDRS	      = ICHOL/ICHOL_wrapper.h ICHOL/precision.h include/icholRF_apply.h\
		include/readMat.h include/wrCSRmat.h

HDRS	      =

LDFLAGS	      = -O2 -fopenmp -L. -L/usr/local/lib -L./lapack -L./ICHOL

LPFLAGS	      = -O1 -pg -L. -L/usr/local/lib -L./lapack -L./ICHOL

LDDFLAGS      = -O0 -g -fopenmp -L. -L/usr/local/lib -L./lapack -L./ICHOL

#LIBS          = -lICHOL -llapacke -llapack -lblas -lgfortran
LIBS          = -lICHOL -llapacke -llapack -lblas -lgfortran

LINKER	      = $(CXX)

MAKEFILE      = Makefile

OBJS	      = icholRF_apply.o main.o readMat.o

PRINT	      = pr

CAT	      = cat

PROGRAM	      = a.out

SRCS	      = icholRF_apply.cpp main.cpp readMat.cpp

all:		$(PROGRAM)

$(PROGRAM):     $(OBJS) 
		$(LINKER) -o $(PROGRAM) $(LDFLAGS) $(OBJS) $(LIBS)

prof:;          @make "CXXFLAGS=$(PCXXFLAGS)" "LDFLAGS=$(LPFLAGS)" $(PROGRAM)

debug:;         @make "CXXFLAGS=$(DFLAGS)" "LDFLAGS=$(LDDFLAGS)" $(PROGRAM)

clean:;		rm -f $(OBJS) $(PROGRAM).cat $(PROGRAM).ftn

depend:;	mkmf -f $(MAKEFILE) PROGRAM=$(PROGRAM) DEST=$(DEST)

index:;		ctags -wx $(HDRS) $(SRCS)

install:	$(PROGRAM)
		install -s $(PROGRAM) $(DEST)

program:        $(PROGRAM)

profile:        $(PROFILE)

tags:           $(HDRS) $(SRCS); ctags $(HDRS) $(SRCS)

update:		$(DEST)/$(PROGRAM)
###
icholRF_apply.o: include/icholRF_apply.h ICHOL/precision.h
main.o: include/readMat.h include/wrCSRmat.h ICHOL/ICHOL_wrapper.h \
	ICHOL/precision.h include/icholRF_apply.h
