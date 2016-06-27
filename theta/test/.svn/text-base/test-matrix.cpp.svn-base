#include "interface/random.hpp"
#include "interface/matrix.hpp"
#include "test/utils.hpp"
#include "interface/exception.hpp"

#include <boost/test/unit_test.hpp>

using namespace theta;

BOOST_AUTO_TEST_SUITE(matrix_tests)

//test basic matrix properties.
BOOST_AUTO_TEST_CASE(matrix0){
   const size_t N=10;
   Matrix m(N,N);
   BOOST_REQUIRE(m.get_n_rows()==N);
   BOOST_REQUIRE(m.get_n_cols()==N);
   for(size_t i=0; i<N; i++){
      for(size_t j=0; j<N; j++){
         BOOST_REQUIRE(m(i,j)==0);
      }
   }
   m(1,0) = 1.0;
   BOOST_REQUIRE(m(1,0)==1.0);
   //copy constructor:
   Matrix mm(m);
   BOOST_REQUIRE(mm.get_n_rows()==N && mm.get_n_cols()==N);
   for(size_t i=0; i<N; i++){
      for(size_t j=0; j<N; j++){
         BOOST_REQUIRE(m(i,j)==mm(i,j));
      }
   }
}

//test matrix cholesky decomposition
BOOST_AUTO_TEST_CASE(matrix1){
   const unsigned int N=10;
   //create a positive definite matrix n:
   std::auto_ptr<RandomSource> rnd_src(new RandomSourceTaus());
   Random rnd(rnd_src);//no setSeed to keep it the same every time ...
   Matrix m(N,N);
   for(size_t i=0; i<N; i++){
      for(size_t j=0; j<N;j++){
        m(i,j) = rnd.uniform();
      }
   }
   Matrix n(N,N);
   for(size_t i=0; i<N; i++){
      for(size_t j=0; j<N;j++){
         double sum = 0;
         for(size_t k=0; k<N; k++){
	    sum += m(i,k) * m(j,k);
	 }
         n(i,j) = sum;
      }
   }
   Matrix n_old(n);
   n.cholesky_decomposition();
   Matrix l(n);
   Matrix lt(n);
   for(size_t i=0; i<N; i++){
      for(size_t j=i+1; j<N; j++){
         lt(j,i) = l(i,j) = 0.0;
      }
   }
   //compare l * lt to the saved n, n_old:
   for(size_t i=0; i<N; i++){
      for(size_t j=0; j<N; j++){
         double n_ij = 0;
	 for(size_t k=0; k<N; k++){
	    n_ij += l(i,k) * lt(k,j);
	 }
	 BOOST_REQUIRE(close_to(n_ij, n_old(i,j), 1));
      }
   }
   //test errors:
   n = n_old;
   n(0,0) = -1.0;
   bool exception = false;
   try{
      n.cholesky_decomposition();
   }
   catch(std::range_error & e){
      exception = true;
   }
   BOOST_REQUIRE(exception);
}

//test matrix inversion, using cholesky:
BOOST_AUTO_TEST_CASE(matrix2){
   const unsigned int N=10;
   //create a positive definite matrix n (see matrix 1 test case):
   std::auto_ptr<RandomSource> rnd_src(new RandomSourceTaus());
   Random rnd(rnd_src);//no setSeed to keep it the same every time ...
   Matrix m(N,N);
   for(size_t i=0; i<N; i++){
      for(size_t j=0; j<N;j++){
        m(i,j) = rnd.uniform();
      }
   }
   Matrix n(N,N);
   for(size_t i=0; i<N; i++){
      for(size_t j=0; j<N;j++){
         double sum = 0;
         for(size_t k=0; k<N; k++){
	    sum += m(i,k) * m(j,k);
	 }
         n(i,j) = sum;
      }
   }
   Matrix n_old(n);
   n.invert_cholesky();
   Matrix unity(N,N);
   for(size_t i=0; i<N; i++){
       for(size_t j=0; j<N; j++){
           for(size_t k=0; k<N; k++){
               unity(i,j) += n(i,k) * n_old(k,j);
           }
       }
   }
   //test matrix for unity:
   for(size_t i=0; i<N; i++){
       for(size_t j=0; j<N; j++){
           if(i==j) BOOST_CHECK(close_to(unity(i,j), 1, 100));
           else BOOST_CHECK(close_to(unity(i,j), 0, 100));
       }
   }
}


BOOST_AUTO_TEST_SUITE_END()
