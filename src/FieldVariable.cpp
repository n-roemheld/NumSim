#include "FieldVariable.h"
#include <iostream>
#include <math.h>

  //!constructor
  FieldVariable::FieldVariable(std::array< int, 2 > size, std::array< double, 2 > origin, std::array< double, 2 > meshWidth) :
    Array2D(size), origin_(origin), meshWidth_(meshWidth)
    {};

  //!get the value at the Cartesian coordinate (x,y). The value is linearly interpolated between stored points.
  // works only for points within the domain
  double FieldVariable::interpolateAt(double x, double y) const
    {
      double dist_i = (x - origin_[0])/meshWidth_[0];
      double dist_j = (y - origin_[1])/meshWidth_[1];

      // find smaller neighbouring indices
      int i_left = int((dist_i));
      int j_lower = int((dist_j));

      // manipulate indices to be within domain
      if ( i_left < 0)
      {
        i_left = 0;
      }
      if ( j_lower < 0)
      {
        j_lower = 0;
      }

      // compute weights for bilinear interpolation
      double alpha_x = 1-(dist_i - i_left);
      double alpha_y = 1-(dist_j - j_lower);

      // find bigger neighbouring indices
      double i_right;
      double j_upper;

      // manipulate indices to be within domain
      if (i_left == size_[0]-1)
      {
        i_right = i_left;
      }
      else
      {
        i_right = i_left+1;
      }
      if (j_lower == size_[1]-1)
      {
        j_upper = j_lower;
      }
      else
      {
        j_upper = j_lower +1;
      }

      // compute bilinear interpolation
      double left_lower_val = Array2D::operator()(i_left, j_lower);
      double left_upper_val = Array2D::operator()(i_left, j_upper);
      double right_lower_val = Array2D::operator()(i_right, j_lower);
      double right_upper_val = Array2D::operator()(i_right, j_upper);

      double mean_left = alpha_y * left_lower_val + (1-alpha_y) * left_upper_val;
      double mean_right = alpha_y * right_lower_val + (1-alpha_y) * right_upper_val;

      return alpha_x*mean_left + (1-alpha_x)*mean_right;


    };

    void FieldVariable::setToZero()
    {
      for(int i = 0; i < size_[0]; i++)
      {
        for(int j = 0; j < size_[1]; j++)
        {
          Array2D::operator()(i,j) = 0;
        }
      }
    }

    // const void* FieldVariable::data()
    // {
    //   return &data_;
    // }

    void* FieldVariable::data()
    {
      return &data_;
    }
