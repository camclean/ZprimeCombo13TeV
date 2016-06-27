#ifndef MATRIX_H
#define MATRIX_H

#include <vector>
#include <cstddef>
#include "interface/exception.hpp"

namespace theta {

/** \brief A matrix of doubles
 *
 * Beyond mere data, this matrix class provides algorithms commonly required for
 * positive definite matrices: cholesky decomposition and matrix inversion based
 * on cholesky decomposition. A usage example for these transformation is the random number generation
 * of a multivariate normal distribution with a non-trivial covariance matrix.
 */
class Matrix {
private:
    size_t rows, cols;
    std::vector<double> elements;

public:
    /** \brief Construct a new Matrix with the given dimensions.
     * 
     * The matrix will be initialized with zero.
     * 
     * \param rows The number of rows.
     * \param cols The number of columns.
     * */
    explicit Matrix(size_t rows = 0, size_t cols = 0);


    /** \brief resize and clear the matrix.
     *
     * Resizes the matrix and fills it with 0.
     */
    void reset(size_t rows, size_t cols);

    /** \brief Detemines the cholesky decomposition of this matrix in-place.
     * 
     * Throws a range_error if this Matrix is not square, not symmetric, or not positive definite.
     */
    void cholesky_decomposition();

    /** \brief Determines the inverse of this Matrix in-place, using a cholesky decomposition.
     * 
     * Throws the same exceptions as \sa cholesky_decomposition
     */
    void invert_cholesky();

    /** \brief Access matrix elements.
     * 
     *  \param row Row index, 0 &lt;= row &lt; getRows()
     *  \param col Column index, 0 &lt;=
     *  \return A reference to the specified matrix element, which can be used as r-value also. 
     */
    double& operator()(size_t row, size_t col) {
        return elements[cols * row + col];
    }

    /** \brief Access matrix elements, const overload.
     * 
     *  \param row Row index, 0 &lt;= row &lt; getRows()
     *  \param col Column index, 0 &lt;=
     *  \return The specified matrix element. 
     */
    double operator()(size_t row, size_t col) const {
        return elements[cols * row + col];
    }
    
    /** \brief Add other matrix to this
     *
     * \c other must have the same number of rows and columns. Otherwise,
     * an invalid_argument exception is thrown.
     */
    const Matrix & operator+=(const Matrix & other){
        if(rows!=other.rows || cols!=other.cols) throw std::invalid_argument("Matrix::operator+= called for matrices of different dimensions");
        for(size_t i=0; i<elements.size(); ++i){
            elements[i] += other.elements[i];
        }
        return *this;
    }
    
    /** \brief Add two matrices
     *
     * \c other must have the same number of rows and columns. Otherwise,
     * an invalid_argument exception is thrown.
     */
    Matrix operator+(const Matrix & other) const{
        Matrix result(*this);
        result += other;
        return result;
    }

    /** \brief The number of rows.
     * 
     * \return The number of rows this Matrix was constructed with.
     * \sa getCols
     */
    size_t get_n_rows() const {
        return rows;
    }

    /** \brief The number of columns.
     * 
     * \return The number of columns this Matrix was constructed with.
     * \sa getRows
     */
    size_t get_n_cols() const {
        return cols;
    }
};

}

#endif
