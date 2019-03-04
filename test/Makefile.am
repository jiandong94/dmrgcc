SUBDIRS = 

bin_PROGRAMS = test_input \
               test_real_matrix \
			   test_real_matrix_block \
               test_real_tensor_lattice \
			   test_real_tensor_space \
               test_complex \
               test_complex_matrix \
               test_complex_matrix_block \
               test_complex_tensor_lattice \
               test_complex_tensor_space 

test_input_SOURCES = test_input.cc
test_real_matrix_SOURCES = test_real_matrix.cc
test_real_matrix_block_SOURCES = test_real_matrix_block.cc
test_real_tensor_lattice_SOURCES = test_real_tensor_lattice.cc
test_real_tensor_space_SOURCES = test_real_tensor_space.cc
test_complex_SOURCES = test_complex.cc
test_complex_matrix_SOURCES = test_complex_matrix.cc
test_complex_matrix_block_SOURCES = test_complex_matrix_block.cc
test_complex_tensor_lattice_SOURCES = test_complex_tensor_lattice.cc
test_complex_tensor_space_SOURCES = test_complex_tensor_space.cc

test_input_LDADD = -L$(SRCDIR)/util -lGeneral
test_real_matrix_LDADD = -L$(SRCDIR)/util -L$(SRCDIR)/tensor \
                         -lTensor -lGeneral @BLAS_LIBS@ @LAPACK_LIBS@ @MKL_LIBS@
test_real_matrix_block_LDADD = -L$(SRCDIR)/util -L$(SRCDIR)/tensor \
                               -lTensor -lGeneral @BLAS_LIBS@ @LAPACK_LIBS@ @MKL_LIBS@
test_real_tensor_lattice_LDADD = -L$(SRCDIR)/util -L$(SRCDIR)/tensor -L$(SRCDIR)/dmrg/realdmrg \
                                 -lRealDMRG -lTensor -lGeneral @BLAS_LIBS@ @LAPACK_LIBS@ @MKL_LIBS@
test_real_tensor_space_LDADD = -L$(SRCDIR)/util -L$(SRCDIR)/tensor -L$(SRCDIR)/dmrg/realdmrg -L$(top_srcdir)/model/space \
                               -lSpace -lRealDMRG -lTensor -lGeneral @BLAS_LIBS@ @LAPACK_LIBS@ @MKL_LIBS@
test_complex_LDADD = -L$(SRCDIR)/util -L$(SRCDIR)/tensor \
                     -lTensor -lGeneral @BLAS_LIBS@ @LAPACK_LIBS@ @MKL_LIBS@
test_complex_matrix_LDADD = -L$(SRCDIR)/util -L$(SRCDIR)/tensor \
                            -lTensor -lGeneral @BLAS_LIBS@ @LAPACK_LIBS@ @MKL_LIBS@
test_complex_matrix_block_LDADD = -L$(SRCDIR)/util -L$(SRCDIR)/tensor \
                                  -lTensor -lGeneral @BLAS_LIBS@ @LAPACK_LIBS@ @MKL_LIBS@
test_complex_tensor_lattice_LDADD = -L$(SRCDIR)/util -L$(SRCDIR)/tensor -L$(SRCDIR)/dmrg/complexdmrg \
                                 -lComplexDMRG -lTensor -lGeneral @BLAS_LIBS@ @LAPACK_LIBS@ @MKL_LIBS@
test_complex_tensor_space_LDADD = -L$(SRCDIR)/util -L$(SRCDIR)/tensor -L$(SRCDIR)/dmrg/complexdmrg -L$(top_srcdir)/model/space \
                               -lSpace -lComplexDMRG -lTensor -lGeneral @BLAS_LIBS@ @LAPACK_LIBS@ @MKL_LIBS@
