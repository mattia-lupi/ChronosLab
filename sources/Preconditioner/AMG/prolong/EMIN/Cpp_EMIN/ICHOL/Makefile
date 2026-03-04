FC            = gfortran
FFLAGS        = -O2 -fopenmp -fPIC
DFFLAGS       = -fopenmp -fPIC -fbacktrace
ARFLAGS	      = 

CFLAGS	      =

DEST	      = .

EXTHDRS	      =

HDRS	      = CPT_TriMatProd.h CPT_TriMatProd_int.h Complete_Matrix_int.h\
		ICHOL_wrapper.h IRheapsort_int.h IcholRF_apply_int.h\
		Ichol_RootFree_shift_int.h RI_SortSplit_int.h SPRS_3MP_int.h\
		Swap_int.h count_rowterms_compl_int.h free_fptr.h\
		mkiat_FTglo_compl_int.h mkiat_FTloc_compl_int.h\
		mvcoef_compl_int.h updt_list_int.h

INSTALL	      = cp

LIBRARY	      = libICHOL.a

MAKEFILE      = Makefile

OBJS	      = A0_precision.o A_resize.o C2FOR.o CPT_TriMatProd.o Complete_Matrix.o\
		ICHOL_wrapper.o IcholRF_apply.o Ichol_RootFree_shift.o\
		RI_SortSplit.o SPRS_3MP.o Swap.o count_rowterms_compl.o\
		free_fptr.o mkiat_FTglo_compl.o mkiat_FTloc_compl.o\
		mvcoef_compl.o updt_list.o IRheapsort.o

PRINT	      = pr

SHELL	      = /bin/sh

SRCS	      = A0_precision.F90 A_resize.f90 C2FOR.f90 CPT_TriMatProd.f90 Complete_Matrix.f90\
		ICHOL_wrapper.f90 IcholRF_apply.f90 Ichol_RootFree_shift.f90\
		RI_SortSplit.f90 SPRS_3MP.f90 Swap.f90 count_rowterms_compl.f90\
		free_fptr.f90 mkiat_FTglo_compl.f90 mkiat_FTloc_compl.f90\
		mvcoef_compl.f90 updt_list.f90 IRheapsort.F90

.SUFFIXES: .prj .f90 .F90

.f90.o:
	$(FC) $(FFLAGS) -c $<

.F90.o:
	$(FC) $(FFLAGS) -c $<

SYSHDRS	      =

all:		$(LIBRARY)

$(LIBRARY):	$(OBJS)
		@echo "Loading $(LIBRARY) ..."
		@ar $(ARFLAGS) cr $(LIBRARY) $(OBJS)
		@ranlib $(LIBRARY)
		@echo "done"

debug:;		@make "FFLAGS=$(DFFLAGS) -g -C" $(LIBRARY)

clean:;		@rm -f $(OBJS) core

clobber:;	@rm -f $(OBJS) $(LIBRARY) core tags

depend:;	@mkmf -f $(MAKEFILE)

echo:;		@echo $(HDRS) $(SRCS)

extract:;	@ar x $(DEST)/$(LIBRARY)

index:;		@ctags -wx $(HDRS) $(SRCS)

install:	$(LIBRARY)
		@echo Installing $(LIBRARY) in $(DEST)
		@if [ "$(DEST)" != . -a "$(DEST)" != "" ]; then \
		(rm -f $(DEST)/$(LIBRARY); $(INSTALL) $(LIBRARY) $(DEST)); fi

print:;		@$(PRINT) $(HDRS) $(SRCS)

tags:           $(HDRS) $(SRCS); @ctags $(HDRS) $(SRCS)
