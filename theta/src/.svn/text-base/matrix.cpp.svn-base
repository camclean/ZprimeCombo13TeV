#include "interface/matrix.hpp"
#include "interface/exception.hpp"

#include <cmath>
#include <algorithm>
#include <limits>

using namespace theta;
using namespace std;

Matrix::Matrix(size_t r, size_t c) : rows(r), cols(c==0?r:c), elements(rows*cols, 0) {}

void Matrix::reset(size_t r, size_t c){
    rows = r;
    cols = c;
    elements.resize(r*c);
    std::fill(elements.begin(), elements.end(), 0);
}

void Matrix::invert_cholesky(){
    cholesky_decomposition();
    Matrix result(rows, cols);
    double max_tii, min_tii;
    max_tii = -numeric_limits<double>::infinity();
    min_tii = numeric_limits<double>::infinity();
    //column i and row j of the inverse:
    for(size_t i=0; i<rows; i++){
        double t_ii = (*this)(i,i);
        if(t_ii==0.0){
           throw range_error("Matrix inversion not possible: division by zero");
        }
        max_tii = std::max(max_tii, fabs(t_ii));
        min_tii = std::min(min_tii, fabs(t_ii));
        //the diagonals of the inverse are the inverse of the diagonals:
        result(i,i) = 1.0/t_ii;
        for(size_t j=i+1; j<cols; j++){
            //inverse(j,i)
            result(j,i) = 0;
            for(size_t k=0; k<j; k++){
                result(j,i) += (*this)(j,k) * result(k,i);
            }
            result(j,i) *= -1.0/(*this)(j,j);
        }
    }
    
    if(max_tii / min_tii > 1E12){
        throw range_error("Matrix inversion: problem is very bad (largest / smallest eigenvalue > 1E12)");
    }
    
    //now, result=L^-1 where L is the lower triangular matrix from the cholesky decomposition,
    //i.e. (originally) this = L * L^t.
    //So (this)^-1 = (L^-1)^t * L^-1
    for(size_t i=0; i<rows; i++){
        for(size_t j=0; j<rows; j++){
            (*this)(i,j) = 0.0;
            for(size_t k=0; k<rows; k++){
                (*this)(i,j) += result(k, i) * result(k,j);
            }
        }
    }
}


void Matrix::cholesky_decomposition(){
   if(cols!=rows) throw range_error("cholesky: not square");
   Matrix & m = *this;
   if(m(0,0) < 0.0) throw range_error("cholesky: not positive definite");
   double l_00 = m(0,0) = sqrt(m(0,0));
   
   if(rows > 1){
      double m_10 = m(1,0);
      double m_11 = m(1,1);
      double l_10 = m_10 / l_00;
      double diag = m_11 - l_10 * l_10;
      if(diag < 0) throw range_error("cholesky: not positive definite");
      m(1,1) = sqrt(diag);
      m(1,0) = l_10;
   }
   for(size_t k=2; k<rows; k++){
      double m_kk = m(k,k);
      for(size_t i=0; i<k; i++){
         double sum = 0;
	 double m_ki = m(k,i);
	 double m_ii = m(i,i);
	 if(i>0){
	    for(size_t j=0; j<i; j++){
	       sum += m(i,j)*m(k,j);
	    }
	 }
	 m_ki = (m_ki - sum) / m_ii;
	 m(k,i) = m_ki;
      }
      double sum = 0;
      for(size_t i=0; i<k; i++){
         sum += m(k,i) * m(k,i);
      }
      double diag = m_kk - sum;
      double l_kk = sqrt(diag);
      if(diag < 0) throw range_error("cholesky: not positive definite");
      m(k,k) = l_kk;
   }
   
   for(size_t i=1; i<cols; i++){
      for(size_t j=0; j<i; j++){
         m(j,i) = m(i,j);
      }
   }
}
