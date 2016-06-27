#include "interface/histogram.hpp"
#include "interface/histogram-function.hpp"
#include "interface/phys.hpp"
#include "test/utils.hpp"
#include "interface/random.hpp"
#include "interface/exception.hpp"


#include <boost/test/unit_test.hpp>
#include <iostream>


using namespace std;
using namespace theta;

// general note: use odd bin numbers to test SSE implementation for which
// an odd number of bins is a special case.

BOOST_AUTO_TEST_SUITE(histogram_tests)

// test the utils fast functions
BOOST_AUTO_TEST_CASE(utils_test){
	double* x=0, * y=0, * y2=0;
	const unsigned int align = 16;
	for(int N = 1; N < 30; ++N){
		int N_alloc = N + (N % 2);
		assert(0==posix_memalign((void**)(&x), align, sizeof(double) * N_alloc));
    	assert(0==posix_memalign((void**)(&y), align, sizeof(double) * N_alloc));
    	assert(0==posix_memalign((void**)(&y2), align, sizeof(double) * N_alloc));

    	for(int j=0; j<N; ++j){
    		x[j] = N + j*j;
    		y2[j] = j;
    	}
    	utils::copy_fast(y, x, N);
    	for(int j=0; j<N; ++j){
			BOOST_REQUIRE(x[j] == y[j]);
		}
    	utils::add_fast(x, y2, N);
    	for(int j=0; j<N; ++j){
			BOOST_REQUIRE(close_to_relative(y[j] + y2[j], x[j]));
		}
    	utils::copy_fast(x, y, N);
    	double c = M_PI;
    	utils::add_fast_with_coeff(x, y2, c, N);
    	for(int j=0; j<N; ++j){
			BOOST_REQUIRE(close_to_relative(y[j] + c * y2[j], x[j]));
		}
    	utils::copy_fast(x, y, N);
    	utils::mul_fast(x, c, N);
    	for(int j=0; j<N; ++j){
			BOOST_REQUIRE(close_to_relative(c * y[j], x[j]));
		}

    	free(x);
    	free(y);
    	free(y2);
	}
}

//test constructors and copy assignment
BOOST_AUTO_TEST_CASE(ctest){
   //default construction:
    Histogram1D h_def;
    BOOST_CHECK(h_def.get_nbins()==0);
    BOOST_CHECK(h_def.get_data()==0);
    
   const size_t nbins = 101;
   Histogram1D m(nbins, -1, 1);
   BOOST_CHECK(m.get_nbins()==nbins);
   BOOST_CHECK(m.get_xmin()==-1);
   BOOST_CHECK(m.get_xmax()==1);
   BOOST_CHECK(m.get_sum()==0);
   //fill a bit:
   for(size_t i=0; i<nbins; i++){
       m.set(i, i*i);
   }

   // try different byte offsets as there are a lot of optimizations:
   for(size_t off=0; off < 67; ++off){
	   const size_t nb = 32 + off;
	   Histogram1D h0(nb, -1, 1);
	   BOOST_REQUIRE(h0.get_nbins()==nb);
	   for(size_t i=0; i<nb; i++){
		   h0.set(i, i*i - i);
	   }
	   Histogram1D h0_init(h0);
	   BOOST_CHECK(histos_equal(h0, h0_init));
	   Histogram1D h0_assign(100, -1, 1);
	   h0_assign = h0;
	   BOOST_CHECK(histos_equal(h0, h0_assign));
   }

   //copy constructor:
   Histogram1D mcopy(m);
   BOOST_REQUIRE(mcopy.get_data()!=m.get_data());
   BOOST_CHECK(histos_equal(m, mcopy));

   //copy assignment:
   Histogram1D m200(2 * nbins, -1, 1);
   BOOST_CHECK(m200.get_nbins()==2*nbins);
   m200 = m;
   BOOST_CHECK(m200.get_nbins()==nbins);
   BOOST_CHECK(histos_equal(m200, m));

   //copy assignment with empty histo:
   Histogram1D h_empty;
   h_empty = m;
   BOOST_CHECK(histos_equal(h_empty, m));
   
   Histogram1D h_empty2;
   BOOST_CHECK(h_empty2.get_nbins()==0);
   BOOST_CHECK(h_empty2.get_data()==0);
   h_empty = h_empty2;
   BOOST_CHECK(h_empty.get_nbins()==0);
   BOOST_CHECK(h_empty.get_data()==0);
   
   bool exception = false;
   try{
      Histogram1D m2(nbins, 1, 0);
   }
   catch(invalid_argument & ex){
      exception = true;
   }
   BOOST_CHECK(exception);
   
   exception = false;
   try{
      Histogram1D m2(nbins, -1, -1);
   }
   catch(invalid_argument & ex){
      exception = true;
   }
   BOOST_CHECK(exception);
}

BOOST_AUTO_TEST_CASE(add_coeff){
	const size_t nbins=100;
	Histogram1D m1(nbins, -1, 1);
	Histogram1D m2(nbins, -1, 1);
	Histogram1D m_expected(nbins, -1, 1);
	const double c = M_PI;
	for(size_t i=0; i<nbins; i++){
		double v1 = 0.33*i;
		double v2 = 10. - i;
		m1.set(i, v1);
		m2.set(i, v2);
		m_expected.set(i, v1 + c*v2);
	}
	{
		Histogram1D s(m1);
		s.add_with_coeff(c, m2);
		BOOST_CHECK(histos_equal(s, m_expected));
	}

    const double c2 = 4.3;
    Histogram1D m3(nbins, -1, 1);
    for(size_t i=0; i<nbins; i++){
    	double v = 2.*i*i;
    	m3.set(i, v);
    	double expected = m1.get(i) + c * m2.get(i) + c2 * v;
    	m_expected.set(i, expected);
    }
    {
    	Histogram1D s(m1);
    	s.add_with_coeff2(c, m2, c2, m3);
    	BOOST_CHECK(histos_equal(s, m_expected));
    	s = m1;
    	s.add_with_coeff(c, m2);
    	s.add_with_coeff(c2, m3);
    	BOOST_CHECK(histos_equal(s, m_expected));
    }
}


//test get, set, fill:
BOOST_AUTO_TEST_CASE(getset){
    const size_t nbins=100;
    Histogram1D m(nbins, -1, 1);
    volatile double sum = 0.0;
    for(size_t i=0; i<nbins; i++){
        double a = sqrt(i+0.0);
        sum += a;
        m.set(i, a);
    }

    BOOST_CHECK(close_to_relative(m.get_sum(),sum));

    for(size_t i=0; i<nbins; i++){
        double a = sqrt(i+0.0);
        BOOST_CHECK(m.get(i) == a);
    }

    //fill:
    double content = m.get(0);
    m.fill(-0.999, 1.7);
    content += 1.7;
    BOOST_CHECK(content==m.get(0));
    sum += 1.7;
    BOOST_CHECK(close_to_relative(m.get_sum(),sum));

    //fill in underflow, content should not change
    content = m.get(0);
    double delta = 10.032;
    m.fill(-1.001, delta);
    BOOST_CHECK(content==m.get(0));
    BOOST_CHECK(close_to_relative(m.get_sum(),sum));

    //fill in overflow, content should not change:
    content = m.get(nbins-1);
    delta = 7.032;
    m.fill(1.001, delta);
    BOOST_CHECK(content==m.get(nbins-1));
    BOOST_CHECK(close_to_relative(m.get_sum(),sum));
}

//test +=
BOOST_AUTO_TEST_CASE(test_plus){
    std::auto_ptr<RandomSource> rnd_src(new RandomSourceTaus());
    Random rnd(rnd_src);
    const size_t nbins = 101;
    Histogram1D m0(nbins, 0, 1);
    Histogram1D m1(m0);
    Histogram1D m_expected(m0);
    for(size_t i=0; i<nbins; ++i){
        volatile double g0 = rnd.get();
        m0.set(i, g0);
        volatile double g1 = rnd.get();
        m1.set(i, g1);
        g0 += g1;
        m_expected.set(i, g0);
    }
    //m0+=m1 should add the histogram m1 to m0, leaving m1 untouched ...
    Histogram1D m1_before(m1);
    Histogram1D m0_before(m0);
    m0+=m1;
    //the sum of weights could be slightly off, but check_histos_equal does handle this.
    BOOST_CHECK(histos_equal(m0, m_expected));
    BOOST_CHECK(histos_equal(m1, m1_before));
    //... and it should commute:
    m1+=m0_before;
    BOOST_CHECK(histos_equal(m1, m_expected));
    //BOOST_CHECKPOINT("test_plus m1, m0");
    BOOST_CHECK(histos_equal(m1, m0));
}

//test *=
BOOST_AUTO_TEST_CASE(test_multiply){
    std::auto_ptr<RandomSource> rnd_src(new RandomSourceTaus());
    Random rnd(rnd_src);
    const size_t nbins = 101;
    Histogram1D m0(nbins, 0, 1);
    Histogram1D m1(m0);
    Histogram1D m0m1(m0);
    Histogram1D m0factor_expected(m0);
    double factor = rnd.get();
    for(size_t i=0; i<nbins; i++){
        double g0 = rnd.get();
        m0.set(i, g0);
        m0factor_expected.set(i, g0*factor);
        double g1 = rnd.get();
        m1.set(i, g1);
        m0m1.set(i, g0*g1);
    }
    Histogram1D m0_before(m0);
    m0*=factor;
    BOOST_CHECK(histos_equal(m0, m0factor_expected));
    m1*=m0_before;
    BOOST_CHECK(histos_equal(m1, m0m1));
    bool exception = false;
    //check error behaviour:
    try{
       Histogram1D m2(nbins + 1, 0, 1);
       m0 *= m2;
    }
    catch(invalid_argument & ex){
       exception = true;
    }
    BOOST_REQUIRE(exception);
}


BOOST_AUTO_TEST_SUITE_END()
