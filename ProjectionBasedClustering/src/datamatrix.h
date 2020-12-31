/*Copyright (C) Kristian Nybo, Jarkko Venna
 *
 *This software is released under the GNU Lesser General Public
 *License. See the included file LICENSE for details.*/

#ifndef DATAMATRIX_HH
#define DATAMATRIX_HH

#include <Rcpp.h>
#include "exception.h"
#include <iostream>

namespace dredviz {
class DataMatrix
{
private:

  double **data;
  size_t rows;
  size_t cols;

public:

  /* Creates an UNINITIALIZED 1*1 matrix. */

  DataMatrix ();


  /* Creates an UNINITIALIZED row * col matrix. */

  DataMatrix (size_t row, size_t col);

  /* Creates a matrix out of an R object */
  DataMatrix (const Rcpp::NumericMatrix mat);

  DataMatrix (const DataMatrix &);

  virtual ~ DataMatrix ();

  size_t getRows () const
  {
    return this->rows;
  }
  size_t getCols () const
  {
    return this->cols;
  }

  double &operator () (size_t i, size_t j)
  {
    return data[i][j];
  }

  double operator () (size_t i, size_t j) const
  {
    return data[i][j];
  }


  DataMatrix & operator+= (const DataMatrix & other)
  {
    for (size_t i = 0; i < other.getRows (); i++)
      for (size_t j = 0; j < other.getCols (); j++)
        (*this) (i, j) += other (i, j);

    return *this;
  }


  DataMatrix & operator-= (const DataMatrix & other)
  {
    for (size_t i = 0; i < other.getRows (); i++)
      for (size_t j = 0; j < other.getCols (); j++)
        (*this) (i, j) -= other (i, j);

    return *this;
  }


  DataMatrix & operator= (const DataMatrix & other);

  virtual double getMin () const;
  virtual double getMax () const;
  virtual double getAverage () const;

  virtual void scale (double factor);
  void normalize (double desiredMaxElement);

  /* Calculates the dot product of *this and other. */

  double elementwiseProduct (const DataMatrix & other) const;


  friend std::ostream & operator<< (std::ostream & s, const DataMatrix & m);

  // Exceptions:

class OutOfBoundsExc:public Exception
  {
  public:
    OutOfBoundsExc (std::string errMsg):Exception (errMsg)
    {
    }
  };

class ZeroDimExc:public Exception
  {
  public:
    ZeroDimExc (std::string errMsg):Exception (errMsg)
    {
    }
  };
};
}

#endif
