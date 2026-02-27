CXX           = g++

CXXFLAGS      = -O2 -std=c++11 -fopenmp $(INCLUDE)
PCXXFLAGS     = -O1 -pg -std=c++11 $(INCLUDE)

CPLUS_INCLUDE_PATH = /usr/include/x86_64-linux-gnu/c++/4.8

CXX           = CPLUS_INCLUDE_PATH=$(CPLUS_INCLUDE_PATH) g++

INCLUDE       = -I. -I./include

DFLAGS        = -O0 -std=c++11 -g -fopenmp -Wall $(INCLUDE)

DEST	      = .

EXTHDRS	      =

HDRS	      =

LDFLAGS	      = -O2 -fopenmp -L. -L/usr/local/lib

LPFLAGS	      = -O1 -pg -L. -L/usr/local/lib

LDDFLAGS      = -O0 -g -fopenmp -L. -L/usr/local/lib

#LIBS          = -llapacke -llapack -lblas -lgfortran
LIBS          = -llapacke -llapack -lblas -lgfortran

LINKER	      = $(CXX)

MAKEFILE      = Makefile

OBJS	      = main.o

PRINT	      = pr

CAT	      = cat

PROGRAM	      = a.out

SRCS	      = main.cpp

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
