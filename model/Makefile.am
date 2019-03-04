SUBDIRS = space hamiltonian 

DIST_SUBDIRS = space hamiltonian 

bin_PROGRAMS = real_spinless_bose_square_model \
               complex_spinless_bose_square_model \
			   complex_spinful_bose_square_model

#real_spinless_bose_square_model_DEPENDENCIES = $(SRCDIR)/util/General.a $(SRCDIR)/tensor/Tensor.a $(SRCDIR)/dmrg/realdmrg/RealDMRG.a $(top_srcdir)/model/space/Space.a $(top_srcdir)/model/hamiltonian/Hamiltonian.a
real_spinless_bose_square_model_SOURCES = real_spinless_bose_square_model.cc
real_spinless_bose_square_model_LDADD = -L$(SRCDIR)/util -L$(SRCDIR)/tensor -L$(SRCDIR)/dmrg/realdmrg -L$(top_srcdir)/model/space -L$(top_srcdir)/model/hamiltonian\
										-lSpace -lHamiltonian -lRealDMRG -lTensor -lGeneral @BLAS_LIBS@ @LAPACK_LIBS@ @MKL_LIBS@

#complex_spinless_bose_square_model_DEPENDENCIES = $(SRCDIR)/util/General.a $(SRCDIR)/tensor/Tensor.a $(SRCDIR)/dmrg/complexdmrg/ComplexDMRG.a $(top_srcdir)/model/space/Space.a $(top_srcdir)/model/hamiltonian/Hamiltonian.a
complex_spinless_bose_square_model_SOURCES = complex_spinless_bose_square_model.cc
complex_spinless_bose_square_model_LDADD = -L$(SRCDIR)/util -L$(SRCDIR)/tensor -L$(SRCDIR)/dmrg/complexdmrg -L$(top_srcdir)/model/space -L$(top_srcdir)/model/hamiltonian\
										-lSpace -lHamiltonian -lComplexDMRG -lTensor -lGeneral @BLAS_LIBS@ @LAPACK_LIBS@ @MKL_LIBS@

#complex_spinful_bose_square_model_DEPENDENCIES = $(SRCDIR)/util/General.a $(SRCDIR)/tensor/Tensor.a $(SRCDIR)/dmrg/complexdmrg/ComplexDMRG.a $(top_srcdir)/model/space/Space.a $(top_srcdir)/model/hamiltonian/Hamiltonian.a
complex_spinful_bose_square_model_SOURCES = complex_spinful_bose_square_model.cc
complex_spinful_bose_square_model_LDADD = -L$(SRCDIR)/util -L$(SRCDIR)/tensor -L$(SRCDIR)/dmrg/complexdmrg -L$(top_srcdir)/model/space -L$(top_srcdir)/model/hamiltonian\
									      -lSpace -lHamiltonian -lComplexDMRG -lTensor -lGeneral @BLAS_LIBS@ @LAPACK_LIBS@ @MKL_LIBS@
