SUBDIRS = space hamiltonian 

DIST_SUBDIRS = space hamiltonian 

bin_PROGRAMS = real_spinless_bose_square_model \
               complex_spinless_bose_square_model \
			   complex_spinful_bose_square_model

real_spinless_bose_square_model_DEPENDENCIES = $(SRCDIR)/util/libGeneral.a $(SRCDIR)/tensor/libTensor.a $(SRCDIR)/dmrg/realdmrg/libRealDMRG.a $(top_srcdir)/model/space/libSpace.a $(top_srcdir)/model/hamiltonian/libHamiltonian.a
real_spinless_bose_square_model_SOURCES = real_spinless_bose_square_model.cc
real_spinless_bose_square_model_LDADD = -L$(SRCDIR)/util -L$(SRCDIR)/tensor -L$(SRCDIR)/dmrg/realdmrg -L$(top_srcdir)/model/space -L$(top_srcdir)/model/hamiltonian\
										-lSpace -lHamiltonian -lRealDMRG -lTensor -lGeneral @BLAS_LIBS@ @LAPACK_LIBS@ @MKL_LIBS@

complex_spinless_bose_square_model_DEPENDENCIES = $(SRCDIR)/util/libGeneral.a $(SRCDIR)/tensor/libTensor.a $(SRCDIR)/dmrg/complexdmrg/libComplexDMRG.a $(top_srcdir)/model/space/libSpace.a $(top_srcdir)/model/hamiltonian/libHamiltonian.a
complex_spinless_bose_square_model_SOURCES = complex_spinless_bose_square_model.cc
complex_spinless_bose_square_model_LDADD = -L$(SRCDIR)/util -L$(SRCDIR)/tensor -L$(SRCDIR)/dmrg/complexdmrg -L$(top_srcdir)/model/space -L$(top_srcdir)/model/hamiltonian\
										-lSpace -lHamiltonian -lComplexDMRG -lTensor -lGeneral @BLAS_LIBS@ @LAPACK_LIBS@ @MKL_LIBS@

complex_spinful_bose_square_model_DEPENDENCIES = $(SRCDIR)/util/libGeneral.a $(SRCDIR)/tensor/libTensor.a $(SRCDIR)/dmrg/complexdmrg/libComplexDMRG.a $(top_srcdir)/model/space/libSpace.a $(top_srcdir)/model/hamiltonian/libHamiltonian.a
complex_spinful_bose_square_model_SOURCES = complex_spinful_bose_square_model.cc
complex_spinful_bose_square_model_LDADD = -L$(SRCDIR)/util -L$(SRCDIR)/tensor -L$(SRCDIR)/dmrg/complexdmrg -L$(top_srcdir)/model/space -L$(top_srcdir)/model/hamiltonian\
									      -lSpace -lHamiltonian -lComplexDMRG -lTensor -lGeneral @BLAS_LIBS@ @LAPACK_LIBS@ @MKL_LIBS@
