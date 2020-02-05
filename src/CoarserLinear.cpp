#include "CoarserLinear.h"
#include "CoarserDefault.h"

// void CoarserLinear::restrict(std::shared_ptr<MGGrid> mggf, std::shared_ptr<MGGrid> mggc) // pointer to fine and coarse grid  // make input const?
// {
//     for (int j = mggc->pJBegin(); j < mggc->pJEnd(); j++) // coarse grid indices
//     {
//         for (int i = mggc->pIBegin(); i < mggc->pIEnd(); i++)
//         {
//             int i_f = (i-mggc->pIBegin())*2 + mggf->pIBegin(); // fine grid indices (bottom left corner)
//             int j_f = (j-mggc->pJBegin())*2 + mggf->pJBegin();
//             mggc->rhs(i,j) = .25 * (mggf->resVec(i_f,j_f) + mggf->resVec(i_f+1,j_f) + mggf->resVec(i_f,j_f+1) + mggf->resVec(i_f+1,j_f+1));
//             // mggc->p must be zeros
//             mggc->p(i,j) = 0;
//             // mggc->resVec must be initialized, values don't matter
//         }
//     }
// }


void CoarserLinear::restrict(std::shared_ptr<MGGrid> mggf, std::shared_ptr<MGGrid> mggc) // pointer to fine and coarse grid  // make input const?
{
    for (int j = mggc->pJBegin(); j < mggc->pJEnd(); j++) // coarse grid indices
    {
        for (int i = mggc->pIBegin(); i < mggc->pIEnd(); i++)
        {
            int i_f = (i-mggc->pIBegin())*2 + mggf->pIBegin(); // fine grid indices (bottom left corner)
            int j_f = (j-mggc->pJBegin())*2 + mggf->pJBegin();


            mggc->rhs(i,j) = 1/64 * ( 9*mggf->resVec(i_f,j_f) + 9*mggf->resVec(i_f+1,j_f) + 9*mggf->resVec(i_f,j_f+1) + 9*mggf->resVec(i_f+1,j_f+1) 
            + 3*mggf->resVec(i_f-1,j_f) + 3*mggf->resVec(i_f-1,j_f+1) + 3*mggf->resVec(i_f+2,j_f) + 3*mggf->resVec(i_f+2,j_f+1) 
            + 3*mggf->resVec(i_f,j_f-1) + 3*mggf->resVec(i_f+1,j_f-1) + 3*mggf->resVec(i_f,j_f+2) + 3*mggf->resVec(i_f+1,j_f+2)
            + mggf->resVec(i_f-1,j_f-1) + mggf->resVec(i_f+2,j_f-1) + mggf->resVec(i_f-1,j_f+2) + mggf->resVec(i_f+2,j_f+2)
            );
            // mggc->p must be zeros (should be zero already)
            mggc->p(i,j) = 0;
            // mggc->resVec must be initialized, values don't matter
        }
    }
}

// void CoarserLinear::interpolate(std::shared_ptr<MGGrid> mggc, std::shared_ptr<MGGrid> mggf)
// {
//   for (int j = mggc->pJBegin(); j < mggc->pJEnd()-1; j++) // coarse grid indices
//   {
//       for (int i = mggc->pIBegin(); i < mggc->pIEnd()-1; i++)
//       {
//         int i_f = (i-mggc->pIBegin())*2 + mggf->pIBegin(); // fine grid indices (bottom left corner)
//         int j_f = (j-mggc->pJBegin())*2 + mggf->pJBegin();
//         double left_low = mggc->p(i,j);
//         double right_low = mggc->p(i+1,j);
//         double left_high = mggc->p(i,j+1);
//         double right_high = mggc->p(i+1,j+1);
//         mggf->resVec(i_f+1,j_f+1) = 0.75*0.75*left_low + 0.25*0.75*right_low + 0.75*0.25*left_high + 0.25*0.25*right_high;
//         mggf->resVec(i_f+2,j_f+1) = 0.25*0.75*left_low + 0.75*0.75*right_low + 0.25*0.25*left_high + 0.75*0.25*right_high;
//         mggf->resVec(i_f+1,j_f+2) = 0.75*0.25*left_low + 0.25*0.25*right_low + 0.75*0.75*left_high + 0.25*0.75*right_high;
//         mggf->resVec(i_f+2,j_f+2) = 0.25*0.25*left_low + 0.75*0.25*right_low + 0.25*0.75*left_high + 0.75*0.75*right_high;
//       }
//     }
// }

void CoarserLinear::interpolate(std::shared_ptr<MGGrid> mggc, std::shared_ptr<MGGrid> mggf)
{
  for (int j = mggc->pJBegin(); j < mggc->pJEnd(); j++) // coarse grid indices
  {
      for (int i = mggc->pIBegin(); i < mggc->pIEnd(); i++)
      {
        int i_f = (i-mggc->pIBegin())*2 + mggf->pIBegin(); // fine grid indices (bottom left corner)
        int j_f = (j-mggc->pJBegin())*2 + mggf->pJBegin();
        double center = mggc->p(i,j);
        double left = mggc->p(i-1,j);
        double right = mggc->p(i+1,j);
        double lower = mggc->p(i,j-1);
        double upper = mggc->p(i,j+1);
        double left_lower = mggc->p(i-1,j-1);
        double left_upper = mggc->p(i-1,j+1);
        double right_lower = mggc->p(i+1,j-1);
        double right_upper = mggc->p(i+1,j+1);
        mggf->resVec(i_f  ,j_f  ) = .75*.75 * center + .75*.25 * (left +lower) + .25*.25 * left_lower; // left_lower_f
        mggf->resVec(i_f  ,j_f+1) = .75*.75 * center + .75*.25 * (left +upper) + .25*.25 * left_upper; // left_upper_f
        mggf->resVec(i_f+1,j_f  ) = .75*.75 * center + .75*.25 * (right+lower) + .25*.25 * right_lower; // right_lower_f
        mggf->resVec(i_f+1,j_f+1) = .75*.75 * center + .75*.25 * (right+upper) + .25*.25 * right_upper; // right_upper_f
      }
  }
  // Set boundary values on fine grid?!
}
