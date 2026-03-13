CXX           = g++

CXXFLAGS      = -O2 -std=c++11 -fopenmp $(INCLUDE)
PCXXFLAGS     = -O1 -pg -std=c++11 $(INCLUDE)

#CPLUS_INCLUDE_PATH = /usr/include/x86_64-linux-gnu/c++/4.8

CXX           = CPLUS_INCLUDE_PATH=$(CPLUS_INCLUDE_PATH) g++

INCLUDE       = -I. -I./include -I./ICHOL

DFLAGS        = -O0 -std=c++11 -g -fopenmp -Wall $(INCLUDE)

DEST	      = .

EXTHDRS	      = ICHOL/ICHOL_wrapper.h ICHOL/precision.h include/DEFL_PCG.h\
		include/DEFL_PCG_matfree.h include/DebEnv.h\
		include/EMIN_ImpProl.h include/EMIN_matfree.h\
		include/KP_spmat.h include/LinvP_spmat.h include/NLSP_PCG.h\
		include/Orth_Q.h include/Prol_add_Cnodes.h\
		include/Transp_Patt.h include/UinvDP_spmat.h include/ZT_mult.h\
		include/Z_mult.h include/apply_perm.h include/copy_Prol.h\
		include/count_rowterms.h include/cpt_Trace_Acc.h\
		include/cpt_ichol_K.h include/ddot_par.h include/dnrm2_par.h\
		include/gather_B_QR.h include/gather_B_Z.h\
		include/gather_B_dump.h include/gather_K.h include/gather_col.h\
		include/gather_f.h include/get_diagpos.h\
		include/icholRF_apply.h include/ichol_matfree.h\
		include/ichol_mf.h include/inl_blas1.h include/load_Jacobi.h\
		include/mkiat_Tglo.h include/mkiat_Tloc.h include/mult_K_col.h\
		include/mult_K_matfree.h include/mvjcol.h include/parm_EMIN.h\
		include/print_Q.h include/print_Z.h include/readMat.h\
		include/solve_RRQR.h include/spmatv_blk.h include/wrCSRmat.h

HDRS	      =

LDFLAGS	      = -O2 -fopenmp -L. -L/usr/local/lib -L../lapack -L./ICHOL

LPFLAGS	      = -O1 -pg -L. -L/usr/local/lib -L../lapack -L./ICHOL

LDDFLAGS      = -O0 -g -fopenmp -L. -L/usr/local/lib -L../lapack -L./ICHOL

LIBRARIES     = -lICHOL -llapacke -llapack -lblas -lgfortran

LIBS          = $(LIBRARIES)

LINKER	      = $(CXX)

MAKEFILE      = Makefile

OBJS	      = DEFL_PCG.o DEFL_PCG_matfree.o DebEnv.o EMIN_ImpProl.o\
		EMIN_matfree.o KP_spmat.o LinvP_spmat.o NLSP_PCG.o Orth_Q.o\
		Prol_add_Cnodes.o Transp_Patt.o UinvDP_spmat.o ZT_mult.o\
		Z_mult.o apply_perm.o copy_Prol.o count_rowterms.o\
		cpt_Trace_Acc.o cpt_ichol_K.o ddot_par.o dnrm2_par.o\
		gather_B_QR.o gather_B_Z.o gather_B_dump.o gather_K.o\
		gather_col.o gather_f.o icholRF_apply.o ichol_matfree.o\
		load_Jacobi.o main.o mkiat_Tglo.o mkiat_Tloc.o mult_K_col.o\
		mult_K_matfree.o mvjcol.o print_Q.o print_Z.o readMat.o\
		solve_RRQR.o spmatv.o spmatv_blk.o wrCSRmat.o

PRINT	      = pr

CAT	      = cat

PROGRAM	      = a.out

SRCS	      = DEFL_PCG.cpp DEFL_PCG_matfree.cpp DebEnv.cpp EMIN_ImpProl.cpp\
		EMIN_matfree.cpp KP_spmat.cpp LinvP_spmat.cpp NLSP_PCG.cpp\
		Orth_Q.cpp Prol_add_Cnodes.cpp Transp_Patt.cpp UinvDP_spmat.cpp\
		ZT_mult.cpp Z_mult.cpp apply_perm.cpp copy_Prol.cpp\
		count_rowterms.cpp cpt_Trace_Acc.cpp cpt_ichol_K.cpp\
		ddot_par.cpp dnrm2_par.cpp gather_B_QR.cpp gather_B_Z.cpp\
		gather_B_dump.cpp gather_K.cpp gather_col.cpp gather_f.cpp\
		icholRF_apply.cpp ichol_matfree.cpp load_Jacobi.cpp main.cpp\
		mkiat_Tglo.cpp mkiat_Tloc.cpp mult_K_col.cpp mult_K_matfree.cpp\
		mvjcol.cpp print_Q.cpp print_Z.cpp readMat.cpp solve_RRQR.cpp\
		spmatv.cpp spmatv_blk.cpp wrCSRmat.cpp
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
DEFL_PCG.o: include/parm_EMIN.h include/spmatv_blk.h include/apply_perm.h \
	include/Orth_Q.h include/ddot_par.h include/dnrm2_par.h \
	include/icholRF_apply.h
DEFL_PCG_matfree.o: include/parm_EMIN.h include/mult_K_matfree.h \
	include/apply_perm.h include/Orth_Q.h include/ddot_par.h \
	include/dnrm2_par.h include/ichol_mf.h include/ichol_matfree.h \
	include/KP_spmat.h include/LinvP_spmat.h include/UinvDP_spmat.h \
	include/wrCSRmat.h
DebEnv.o: include/DebEnv.h ICHOL/precision.h
EMIN_ImpProl.o: include/parm_EMIN.h include/wrCSRmat.h include/Transp_Patt.h \
	include/EMIN_matfree.h include/gather_K.h include/cpt_ichol_K.h \
	include/copy_Prol.h include/gather_B_QR.h include/gather_B_Z.h \
	include/gather_f.h include/cpt_Trace_Acc.h include/print_Q.h \
	include/print_Z.h include/DEFL_PCG.h include/NLSP_PCG.h \
	include/Prol_add_Cnodes.h
EMIN_matfree.o: include/parm_EMIN.h include/wrCSRmat.h include/Transp_Patt.h \
	include/load_Jacobi.h include/copy_Prol.h include/gather_B_QR.h \
	include/gather_B_dump.h include/gather_f.h include/cpt_Trace_Acc.h \
	include/print_Q.h include/ichol_mf.h include/Prol_add_Cnodes.h \
	include/DEFL_PCG_matfree.h
LinvP_spmat.o: include/LinvP_spmat.h
NLSP_PCG.o: include/parm_EMIN.h include/spmatv_blk.h include/apply_perm.h \
	include/Z_mult.h include/ZT_mult.h include/ddot_par.h \
	include/dnrm2_par.h include/icholRF_apply.h
Transp_Patt.o: include/count_rowterms.h include/mkiat_Tloc.h \
	include/mkiat_Tglo.h include/mvjcol.h
UinvDP_spmat.o: include/UinvDP_spmat.h
count_rowterms.o: include/count_rowterms.h
cpt_ichol_K.o: include/DebEnv.h ICHOL/precision.h include/parm_EMIN.h \
	ICHOL/ICHOL_wrapper.h
gather_B_QR.o: include/inl_blas1.h ICHOL/precision.h include/DebEnv.h \
	include/parm_EMIN.h include/solve_RRQR.h
gather_B_Z.o: include/gather_f.h
gather_B_dump.o: include/DebEnv.h ICHOL/precision.h include/parm_EMIN.h
gather_K.o: include/get_diagpos.h include/gather_col.h
ichol_matfree.o: include/ichol_mf.h include/gather_col.h \
	include/get_diagpos.h ICHOL/ICHOL_wrapper.h ICHOL/precision.h
main.o: include/DebEnv.h ICHOL/precision.h include/readMat.h \
	include/wrCSRmat.h include/EMIN_ImpProl.h
mkiat_Tglo.o: include/mkiat_Tglo.h
mkiat_Tloc.o: include/mkiat_Tloc.h
mult_K_matfree.o: include/mult_K_col.h
