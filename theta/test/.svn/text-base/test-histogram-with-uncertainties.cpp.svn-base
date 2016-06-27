#include "interface/histogram-with-uncertainties.hpp"
#include "test/utils.hpp"
#include "interface/random.hpp"
#include "interface/exception.hpp"


#include <boost/test/unit_test.hpp>
#include <iostream>


using namespace std;
using namespace theta;

namespace{

bool histos_equal(const Histogram1DWithUncertainties & h1, const Histogram1DWithUncertainties & h2){
    if(h1.get_nbins()!=h2.get_nbins()) return false;
    if(h1.get_xmin()!=h2.get_xmin()) return false;
    if(h1.get_xmax()!=h2.get_xmax()) return false;
    const size_t n = h1.get_nbins();
    for(size_t i=0; i<n; i++){
        if(h1.get_value(i)!=h2.get_value(i)) return false;
        if(!close_to_relative(h1.get_uncertainty(i), h2.get_uncertainty(i))) return false;
    }
    return true;
}


}

// general note: use odd bin numbers to test SSE implementation for which
// an odd number of bins is a special case.

BOOST_AUTO_TEST_SUITE(histogram_with_uncertainties_tests)

//test constructors and copy assignment
BOOST_AUTO_TEST_CASE(ctest){
   //default construction:
    Histogram1DWithUncertainties h_def;
    BOOST_CHECK(h_def.get_nbins()==0);
    
   const size_t nbins = 101;
   Histogram1DWithUncertainties m(nbins, -1, 1);
   BOOST_CHECK(m.get_nbins()==nbins);
   BOOST_CHECK(m.get_xmin()==-1);
   BOOST_CHECK(m.get_xmax()==1);
   // check that uncertainties are 0.0:
   for(size_t i=0; i<nbins; i++){
	   BOOST_CHECK(m.get_uncertainty(i)==0.0);
   }
   //fill a bit:
   for(size_t i=0; i<nbins; i++){
       m.set(i, i*i);
   }
   // check that uncertainties are still 0.0:
   for(size_t i=0; i<nbins; i++){
   	   BOOST_CHECK(m.get_uncertainty(i)==0.0);
   }
   const double dv = 15.8;
   m.set(0, m.get(0), dv);
   BOOST_CHECK(m.get_uncertainty(0)==dv);
   for(size_t i=1; i<nbins; i++){
	   BOOST_CHECK(m.get_uncertainty(i)==0.0);
       m.set(i, i*i, i);
   }

   //copy constructor:
   Histogram1DWithUncertainties mcopy(m);
   BOOST_CHECK(histos_equal(m, mcopy));

   //copy assignment:
   Histogram1DWithUncertainties m200(2 * nbins, -1, 1);
   BOOST_CHECK(m200.get_nbins()==2*nbins);
   m200 = m;
   BOOST_CHECK(m200.get_nbins()==nbins);
   BOOST_CHECK(histos_equal(m200, m));

   //copy assignment with empty histo:
   Histogram1DWithUncertainties h_empty;
   h_empty = m;
   BOOST_CHECK(histos_equal(h_empty, m));
   
   Histogram1DWithUncertainties h_empty2;
   BOOST_CHECK(h_empty2.get_nbins()==0);
   h_empty = h_empty2;
   BOOST_CHECK(h_empty.get_nbins()==0);
   
   bool exception = false;
   try{
      Histogram1DWithUncertainties m2(nbins, 1, 0);
   }
   catch(invalid_argument & ex){
      exception = true;
   }
   BOOST_CHECK(exception);
   
   exception = false;
   try{
      Histogram1DWithUncertainties m2(nbins, -1, -1);
   }
   catch(invalid_argument & ex){
      exception = true;
   }
   BOOST_CHECK(exception);
}



//test get, set:
BOOST_AUTO_TEST_CASE(getset){
    const size_t nbins=100;
    Histogram1DWithUncertainties m(nbins, -1, 1);
    for(size_t i=0; i<nbins; i++){
        double a = sqrt(i+0.0);
        double da = (i + 1) * 0.1;
        m.set(i, a, da);
    }

    for(size_t i=0; i<nbins; i++){
        double a = sqrt(i+0.0);
        double da = (i + 1) * 0.1;
        BOOST_CHECK(m.get_value(i) == a);
        BOOST_CHECK(close_to_relative(m.get_uncertainty(i), da));
    }
    
    //zero uncertainty case (optimised ...):
    Histogram1DWithUncertainties m0(nbins, -1, 1);
    for(size_t i=0; i<nbins; i++){
        double a = sqrt(i+0.3);
        m0.set(i, a, 0.0);
    }
    for(size_t i=0; i<nbins; i++){
        double a = sqrt(i+0.3);
        BOOST_CHECK(m0.get_value(i) == a);
        BOOST_CHECK(m0.get_uncertainty(i) == 0.0);
    }
    
}

//test +=
BOOST_AUTO_TEST_CASE(test_plus){
    std::auto_ptr<RandomSource> rnd_src(new RandomSourceTaus());
    Random rnd(rnd_src);
    const size_t nbins = 101;
    Histogram1DWithUncertainties m0(nbins, 0, 1);
    Histogram1DWithUncertainties m1(m0);
    Histogram1DWithUncertainties m_expected(m0);
    for(size_t i=0; i<nbins; ++i){
        //also include negative values to check error propagation ...
        volatile double g0 = rnd.get() - 0.5;
        m0.set(i, g0, 0.25*g0);
        volatile double g1 = rnd.get() - 0.5;
        m1.set(i, g1, 0.5*g1);
        double expected = g0 + g1;
        double error_expected = sqrt(0.25*g0*0.25*g0 + 0.5*g1*0.5*g1);
        m_expected.set(i, expected, error_expected);
    }
    m0 += m1;
    //the sum of weights could be slightly off, but check_histos_equal does handle this.
    BOOST_CHECK(histos_equal(m0, m_expected));
    
    
    // check += with rhs Histogram1D without uncertainties:
    Histogram1D h_nu(nbins, 0, 1);
    for(size_t i=0; i<nbins; ++i){
        h_nu.set(i, rnd.get());
    }
    Histogram1DWithUncertainties m0_copy(m0);
    m0 += h_nu;
    Histogram1D h_expected(nbins, 0, 1);
    for(size_t i=0; i<nbins; ++i){
        double val_expected = m0_copy.get_value(i) + h_nu.get(i);
        h_expected.set(i, val_expected);
        // uncertainty should not change:
        double unc_expected = m0_copy.get_uncertainty(i);
        BOOST_CHECK(m0.get_value(i)==val_expected);
        BOOST_CHECK(close_to_relative(m0.get_uncertainty(i), unc_expected));
    }
    BOOST_CHECK(histos_equal(m0.get_values_histogram(), h_expected));
    
    
    // check += where the rhs has no uncertainties, although it is a Histogram1DWithUncerytainties:
    m0 = m0_copy;
    Histogram1DWithUncertainties h_wu(h_nu);
    m0 += h_wu;
    BOOST_CHECK(histos_equal(m0.get_values_histogram(), h_expected));    
}

BOOST_AUTO_TEST_CASE(add_coeff){
	const size_t nbins=100;
	Histogram1DWithUncertainties m1(nbins, -1, 1);
	Histogram1DWithUncertainties m2(nbins, -1, 1);
	Histogram1DWithUncertainties m_expected(nbins, -1, 1);
	const double c = M_PI;
	for(size_t i=0; i<nbins; i++){
		double v1 = 0.33*i;
		double dv1 = fabs(v1)/4.;
		double v2 = 10. - i;
		double dv2 = fabs(v2)/8.;
		m1.set(i, v1, dv1);
		m2.set(i, v2, dv2);
		m_expected.set(i, v1 + c*v2, sqrt(pow(dv1, 2) + pow(fabs(c) * dv2, 2)));
	}
	{
		Histogram1DWithUncertainties s(m1);
		s.add_with_coeff(c, m2);
		BOOST_CHECK(histos_equal(s, m_expected));
	}
}


BOOST_AUTO_TEST_SUITE_END()
