#include "interface/utils.hpp"

#include <cmath>
#include <limits>

#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>

#ifdef __clang__
#define BOOST_DETAIL_FENV_HPP
#include <fenv.h>
#endif

#include <boost/math/special_functions/gamma.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/filesystem.hpp>


namespace fs = boost::filesystem;


namespace{
//note: these functions are useful to have compatibility
// with both V2 and V3 of boost::filesystem
// as in V2, path::filename returns a string whereas in
// V3, path::filename returns a path.
std::string to_string(const std::string & s){
    return s;
}
std::string to_string(const fs::path & p){
    return p.string();
}

// note: newer boost::filesystem versions have fs::path::canonical, but that's not available on
// older boost versions, so at least resolve "." and ".." in guessing theta_dir:
fs::path resolve_dots(const fs::path & p){
    std::vector<fs::path::const_iterator> components;
    for(fs::path::iterator it=p.begin(); it!=p.end(); ++it){
        if(to_string(*it)=="."){
            continue;
        }
        if(to_string(*it)==".."){
            // remove last component:
            if(components.size() == 0){
                throw std::invalid_argument("to many ..");
            }
            components.pop_back();
            continue;
        }
        components.push_back(it);
    }
    fs::path result;
    for(std::vector<fs::path::const_iterator>::const_iterator cit=components.begin(); cit!=components.end(); ++cit){
        result /= **cit;
    }
    std::string ps = p.string();
    if(ps.size()>0 && ps[ps.size()-1]=='/'){
        result /= "/";
    }
    return result;
}

}



namespace theta{ namespace utils{


std::string replace_theta_dir(const std::string & path) {
   return boost::algorithm::replace_all_copy(path, "$THETA_DIR", theta_dir);
}

void fill_theta_dir(char** argv){
    if(argv!=0){
        fs::path guessed_self_path = argv[0];
        if(fs::is_regular_file(guessed_self_path)){
            fs::path cp = resolve_dots(fs::system_complete(guessed_self_path));
            theta_dir = cp.parent_path().parent_path().string();
            return;
        }
    }
    if(fs::is_symlink("/proc/self/exe")){
        char path[4096];
        ssize_t s = readlink("/proc/self/exe", path, 4096);
        if(s > 0 && s < 4096){
            path[s] = '\0';
            theta_dir = fs::path(path).parent_path().parent_path().string();
        }
    }
}

std::string theta_dir;

int roots_quad(double & x1, double & x2, double b, double c){
    double D = b*b - 4*c;
    if(D < 0){
        x1 = x2 = NAN;
        return 0;
    }
    if(D==0){
        x1 = x2 = -0.5 * b;
        return 1;
    }
    if(b==0.0){
        x1 = -0.5 * sqrt(D);
        x2 = 0.5 * sqrt(D);
    }
    else if(b > 0){
        x1 = 0.5 * (-b - sqrt(D));
        x2 = c / x1;
    }
    else{
        x2 = 0.5 * (-b + sqrt(D));
        x1 = c / x2;
    }
    return 2;
}

double lngamma(double x){
    return boost::math::lgamma(x);
}

/** The inverse Error function in form of the inverse
 * of the cumulative distribution function phi of the 
 * standard gaussian distribution.
 * This is some old fortran routine translation to C. Source:
 * http://lib.stat.cmu.edu/apstat/241
 * via:
 * http://home.online.no/~pjacklam/notes/invnorm/#Other_algorithms
 */
double phi_inverse(double p){
/* ALGORITHM AS241  APPL. STATIST. (1988) VOL. 37, NO. 3
 * Produces the normal deviate Z corresponding to a given lower
 * tail area of P; Z is accurate to about 1 part in 10**16. */
   const double split1 = 0.425, split2 = 5.0, const1 = 0.180625, const2 = 1.6;


/* Coefficients for P close to 0.5 */
   double a[8] = {3.3871328727963666080,
              1.3314166789178437745E+2,
              1.9715909503065514427E+3,
              1.3731693765509461125E+4,
              4.5921953931549871457E+4,
              6.7265770927008700853E+4,
              3.3430575583588128105E+4,
             2.5090809287301226727E+3};
   double b[8] = {0, 4.2313330701600911252E+1,
                   6.8718700749205790830E+2,
                   5.3941960214247511077E+3,
                   2.1213794301586595867E+4,
                   3.9307895800092710610E+4,
                   2.8729085735721942674E+4,
                   5.2264952788528545610E+3};
/* Coefficients for P not close to 0, 0.5 or 1. */
   double c[8] = {1.42343711074968357734E0,
              4.63033784615654529590E0,
              5.76949722146069140550E0,
              3.64784832476320460504E0,
              1.27045825245236838258E0,
              2.41780725177450611770E-1,
              2.27238449892691845833E-2,
              7.74545014278341407640E-4};
   double d[8] = {0, 2.05319162663775882187E0,
                 1.67638483018380384940E0,
                 6.89767334985100004550E-1,
                 1.48103976427480074590E-1,
                 1.51986665636164571966E-2,
                 5.47593808499534494600E-4,
                 1.05075007164441684324E-9 };

/*  Coefficients for P near 0 or 1. */
   double e[8] = {6.65790464350110377720E0,
              5.46378491116411436990E0,
              1.78482653991729133580E0,
              2.96560571828504891230E-1,
              2.65321895265761230930E-2,
              1.24266094738807843860E-3,
              2.71155556874348757815E-5,
              2.01033439929228813265E-7};
   double f[8] = {0, 5.99832206555887937690E-1,
              1.36929880922735805310E-1,
              1.48753612908506148525E-2,
              7.86869131145613259100E-4,
              1.84631831751005468180E-5,
              1.42151175831644588870E-7,
              2.04426310338993978564E-15};

   double q = p - 0.5;
   double r;

   if(fabs(q) <= split1){
      r = const1 - q*q;
      return q * (((((((a[7] * r + a[6]) * r + a[5]) * r + a[4]) * r + a[3])
              * r + a[2]) * r + a[1]) * r + a[0]) /
                 (((((((b[7] * r + b[6]) * r + b[5]) * r + b[4]) * r + b[3])
              * r + b[2]) * r + b[1]) * r + 1.0);
   }

   if(q <= 0.0)  r = p;
   else r = 1.0 - p;
   
   if(r <= 0.0) return std::numeric_limits<double>::signaling_NaN();
   r = sqrt(-log(r));
   double result;
   if(r <= split2){
      r -= const2;
      result = (((((((c[7] * r + c[6]) * r + c[5]) * r + c[4]) * r + c[3])
              * r + c[2]) * r + c[1]) * r + c[0]) /
                (((((((d[7] * r + d[6]) * r + d[5]) * r + d[4]) * r + d[3])
              * r + d[2]) * r + d[1]) * r + 1.0);
   }
   else{
      r -= split2;
      result = (((((((e[7] * r + e[6]) * r + e[5]) * r + e[4]) * r + e[3])
              * r + e[2]) * r + e[1]) * r + e[0]) /
                (((((((f[7] * r + f[6]) * r + f[5]) * r + f[4]) * r + f[3])
              * r + f[2]) * r + f[1]) * r + 1.0);
   }
   if(q < 0) return -result;
   else return result;
}

}}
