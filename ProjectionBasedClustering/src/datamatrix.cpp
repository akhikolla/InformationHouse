/*Copyright (C) Kristian Nybo, Jarkko Venna
 *
 *This software is released under the GNU Lesser General Public
 *License. See the included file LICENSE for details.*/

#include "datamatrix.h"

#include <sstream>
#include <vector>
#include <string>
#include <stdexcept>

namespace dredviz {
DataMatrix::DataMatrix ():rows (1), cols (1)
{
  data = NULL;
  data = new double *[rows];

  for (size_t i = 0; i < rows; i++)
    data[i] = new double[cols];

  /*for(size_t i = 0; i < rows; i++)
     for(size_t j = 0; j < cols; j++)
     data[i][j] = 0.0; */
}


DataMatrix::DataMatrix (size_t rows, size_t cols):rows (rows), cols (cols)
{
  data = NULL;
  if (0 == rows * cols)
    throw (ZeroDimExc ("Tried to create a zero-dimensional matrix."));

  data = new double *[rows];

  for (size_t i = 0; i < rows; i++)
    data[i] = new double[cols];

  /*for(size_t i = 0; i < rows; i++)
     for(size_t j = 0; j < cols; j++)
     data[i][j] = 0.0; */
}

DataMatrix::DataMatrix (const Rcpp::NumericMatrix mat):
    rows (mat.nrow()),
    cols (mat.ncol())
{
  data = NULL;
  data = new double *[rows];

  for (size_t i = 0; i < rows; i++)
    data[i] = new double[cols];

  for (size_t i = 0; i < rows; i++)
    for (size_t j = 0; j < cols; j++)
      (*this) (i, j) = mat (i, j);
}

DataMatrix::DataMatrix (const DataMatrix & m):
    rows (m.getRows ()),
    cols (m.getCols ())
{
  data = NULL;
  data = new double *[rows];

  for (size_t i = 0; i < rows; i++)
    data[i] = new double[cols];


  for (size_t i = 0; i < rows; i++)
    for (size_t j = 0; j < cols; j++)
      (*this) (i, j) = m (i, j);
}


DataMatrix::~DataMatrix ()
{
  if (data == NULL)
    return;
  for (size_t i = 0; i < rows; i++)
    delete[]data[i];

  delete[]data;
}



DataMatrix & DataMatrix::operator= (const DataMatrix & m)
{
  if (this == &m)
    return (*this);

  /* Let's avoid reallocation if the matrices have the same dimensions. */

  if (!(rows == m.getRows () && cols == m.getCols ()))
  {
    for (size_t i = 0; i < rows; i++)
      delete[]data[i];

    delete[]data;


    rows = m.getRows ();
    cols = m.getCols ();
    data = NULL;
    data = new double *[rows];

    for (size_t i = 0; i < rows; i++)
      data[i] = new double[cols];
  }


  for (size_t i = 0; i < rows; i++)
    for (size_t j = 0; j < cols; j++)
      (*this) (i, j) = m (i, j);

  return *this;
}

/*
double& DataMatrix::operator() (size_t i, size_t j)
{
   if(i > (rows - 1) || j > (cols - 1))
    {
      std::ostringstream ost;
      ost << "Tried to access element (" << i << ", " << j << ") in a "
	  << rows << "x" << cols << " matrix.";
	  throw(OutOfBoundsExc(ost.str()));
      }

    return data[i][j];
}*/

/*
double DataMatrix::operator() (size_t i, size_t j) const
{
   if(i > (rows - 1) || j > (cols - 1))
    {
      std::ostringstream ost;
      ost << "Tried to access element (" << i << ", " << j << ") in a "
	  << rows << "x" << cols << " matrix.";
      throw(OutOfBoundsExc(ost.str()));
    }

  return data[i][j];
}*/


std::ostream & operator<< (std::ostream & s, const DataMatrix & m)
{
  for (size_t i = 0; i < m.getRows (); i++)
  {
    for (size_t j = 0; j < m.getCols (); j++)
    {
      s << m (i, j) << " ";
    }
    s << "\n";
  }

  return s;
}


void
DataMatrix::scale (double factor)
{
  for (size_t i = 0; i < rows; i++)
    for (size_t j = 0; j < cols; j++)
      (*this) (i, j) *= factor;
}


double
DataMatrix::getMax () const
{
  double max = (*this) (0, 0);

  for (size_t i = 0; i < rows; i++)
    for (size_t j = 0; j < cols; j++)
      if ((*this) (i, j) > max)
        max = (*this) (i, j);

  return max;
}


double
DataMatrix::getMin () const
{
  double min = (*this) (0, 0);

  for (size_t i = 0; i < rows; i++)
    for (size_t j = 0; j < cols; j++)
      if ((*this) (i, j) < min)
        min = (*this) (i, j);

  return min;
}


double
DataMatrix::elementwiseProduct (const DataMatrix & other) const
{
  double result = 0.0;

  for (size_t i = 0; i < rows; i++)
    for (size_t j = 0; j < cols; j++)
      result += (*this) (i, j) * other (i, j);

  return result;
}


double
DataMatrix::getAverage () const
{
  double average = 0.0;

  for (size_t i = 0; i < rows; i++)
    for (size_t j = 0; j < cols; j++)
      average += (*this) (i, j);

  return average / (rows * cols);
}


void
DataMatrix::normalize (double desiredMaxDistance)
{
  scale (desiredMaxDistance / getMax ());
}
}
