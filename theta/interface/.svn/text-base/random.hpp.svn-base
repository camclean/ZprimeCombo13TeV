#ifndef RANDOM_HPP
#define RANDOM_HPP

#include "interface/decls.hpp"

#include <vector>
#include <memory>
#include <limits>
#include <boost/utility.hpp>
#include <boost/static_assert.hpp>

// the random number algorithms use "unsigned int" as 32 bit field. It is ok to be larger,
// but must not be smaller:
BOOST_STATIC_ASSERT(std::numeric_limits<unsigned int>::digits >= 32);

namespace theta {

/// abstract base class for pseudo random number generators
class RandomSource{
    friend class Random;
    public:
        virtual ~RandomSource(){}
    protected:
        /** \brief Fill the buffer with full 32 bit pseudorandom numbers.
         */
        virtual void fill(std::vector<unsigned int> & buffer) = 0;
        
        /** \brief Set the seed of the generator.
         */
        virtual void set_seed(unsigned int seed) = 0;
};

/** \brief Random number distribution generator
 *
 * This class is used to generate pseudo random number in user code. For the actual generation of "bare"
 * pseudo-randomness, it can use any class derived from RandomSource.
 */
class Random: private boost::noncopyable {
private:
    /* The Gauss Zigurrat method. */
    static const double ytab[128];
    /* tabulated values for 2^24 times x[i]/x[i+1],
    * used to accept for U*x[i+1]<=x[i] without any floating point operations */
    static const long int ktab[128];
    /* tabulated values of 2^{-24}*x[i] */
    static const double wtab[128];
    
    std::auto_ptr<RandomSource> rnd;
    
    //cached random data:
    std::vector<unsigned int> random_data;
    std::vector<unsigned int>::const_iterator current;
    std::vector<unsigned int>::const_iterator end;

    double gauss_zig(double);
    unsigned int poisson_root(double mean);
    
    void fill(){
        rnd->fill(random_data);
        current = random_data.begin();
    }
public:
    
    /** \brief Construct a generator with \c rnd_ as underlying source
     *
     * Takes ownership of rnd_ memory.
     */
    //While larger buffer sizes help very much for small buffers, tests have shown that
    // choosing a buffer sizes larger than around 100 does not increase performance significantly.
    explicit Random(std::auto_ptr<RandomSource> & rnd_): rnd(rnd_), random_data(100),
        current(random_data.end()), end(random_data.end()){
    }
    
    const RandomSource & get_source() const{
        return *rnd;
    }
    
    /** \brief get a random number distributed according to a normal distribution with the given standard deviation
     *
     * Uses the GNU Scientific Library "zigurrat" method. For details and references, see there.
     */
    double gauss(double sigma = 1.0){
        return gauss_zig(sigma);
    }
    
    /** \brief get a random number distributed according to a poisson distribution
     *
     * Uses the GNU Scientific Library "root" method. For details and references, see there.
     */
    unsigned int poisson(double mean){
        return poisson_root(mean);
    }
    
    /** \brief get a 32 bit random number from the generator
     */
    unsigned int get() {
       if(current==end) fill();
       return *(current++);
    }
    
    /** \brief Return a random integer between 0 and n-1 inclusively.
     */
    unsigned int get_uniform_int(unsigned int n);
    
    /// returns a uniform random number on [0,1)
    double uniform(){
        return get() / 4294967296.0;
    }
    
    /** \brief set the random number generator seed
     *
     * This calls RandomSource::set_seed
     */
    void set_seed(unsigned int n){
        rnd->set_seed(n);
        //invalidate the buffer:
        current = end;
    }
};

    /** \brief Tausworthe generator
     *
     * See Pierre L'Ecuyer: "Maximally Equidistributed Combined Tausworthe Generators", Math. Comp. 65, 1996 <br>
     * with seeding modifications described in Pierre L'Ecuyer: "Tables of Maximally Equidistributed Combined LFSR Generators", Math. Comp. 68, 1999.
     */
    class RandomSourceTaus: public RandomSource{
    private:
        unsigned int s1, s2, s3;
    protected:
        //@{
        /** \brief Implement the pure virtual methods from RandomSource
         */
        virtual void fill(std::vector<unsigned int> & buffer);
        virtual void set_seed(unsigned int);
        //@}
    public:
        /// Default constructor; uses the same seed each time
        RandomSourceTaus();
    };
    
    /** \brief The Mersenne Twister generator MT19937
     *
     * See Matsumoto, Makoto and Nishimura, Takuji: "Mersenne twister: a 623-dimensionally equidistributed uniform pseudo-random number generator"
     * ACM Trans. Model. Comput. Simul. 1, 1998
     *
     * This algorithm is approximately 50% slower than the \link RandomSourceTaus Tausworthe generator \endlink.
     */
    class RandomSourceMersenneTwister: public RandomSource{
    private:
        static const int N = 624;
        static const int M = 397;
        static const unsigned int MATRIX_A = 0x9908b0dfU;
        static const unsigned int UPPER_MASK = 0x80000000U;
        static const unsigned int LOWER_MASK = 0x7fffffffU;
        unsigned int mag01[2];//={0x0U, MATRIX_A};
        /* mag01[x] = x * MATRIX_A  for x=0,1 */
        unsigned long mt[N];
        int mti;
    protected:
        //@{
        /** \brief Implement the pure virtual methods from RandomSource
         */
        virtual void fill(std::vector<unsigned int> & buffer);
        virtual void set_seed(unsigned int);
        //@}
    public:
        /// Default constructor; uses the same seed each time
        RandomSourceMersenneTwister();

    };
}


#endif
