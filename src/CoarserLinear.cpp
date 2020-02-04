#include "CoarserLinear.h"
#include "CoarserDefault.h"

void CoarserLinear::restrict(std::shared_ptr<MGGrid> mggf, std::shared_ptr<MGGrid> mggc)
{

}

void CoarserLinear::interpolate(std::shared_ptr<MGGrid> mggc, std::shared_ptr<MGGrid> mggf)
{
  for (int j = mggc->pJBegin(); j < mggc->pJEnd()-1; j++) // coarse grid indices
  {
      for (int i = mggc->pIBegin(); i < mggc->pIEnd()-1; i++)
      {
        int i_f = (i-mggc->pIBegin())*2 + mggf->pIBegin(); // fine grid indices (bottom left corner)
        int j_f = (j-mggc->pJBegin())*2 + mggf->pJBegin();
        int left_low = mggc->p(i,j);
        int right_low = mggc->p(i+1,j);
        int left_high = mggc->p(i,j+1);
        int right_high = mggc->p(i+1,j+1);
        mggf->resVec(i_f+1,j_f+1) = 0.75*0.75*left_low + 0.25*0.75*right_low + 0.75*0.25*left_high + 0.25*0.25*right_high;
        mggf->resVec(i_f+2,j_f+1) = 0.25*0.75*left_low + 0.75*0.75*right_low + 0.25*0.25*left_high + 0.75*0.25*right_high;
        mggf->resVec(i_f+1,j_f+2) = 0.75*0.25*left_low + 0.25*0.25*right_low + 0.75*0.75*left_high + 0.25*0.75*right_high;
        mggf->resVec(i_f+2,j_f+2) = 0.25*0.25*left_low + 0.75*0.25*right_low + 0.25*0.75*left_high + 0.75*0.75*right_high;
      }
    }
}
