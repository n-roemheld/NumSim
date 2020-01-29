#include "CoarserDefault.h"

void CoarserDefault::restrict(std::shared_ptr<MGGrid> mggf, std::shared_ptr<MGGrid> mggc) // pointer to fine and coarse grid  // make input const?
{
    for (int j = mggc->pJBegin(); j < mggc->pJEnd(); j++) // coarse grid indices
    {
        for (int i = mggc->pIBegin(); i < mggc->pIEnd(); i++)
        {
            int i_f = (i-mggc->pIBegin())*2 + mggf->pIBegin(); // fine grid indices (bottom left corner)
            int j_f = (j-mggc->pJBegin())*2 + mggf->pJBegin();
            mggc->rhs(i,j) = .25 * (mggf->resVec(i_f,j_f) + mggf->resVec(i_f+1,j_f) + mggf->resVec(i_f,j_f+1) + mggf->resVec(i_f+1,j_f+1));
            // mggc->p must be zeros
            // mggc->resVec must be initialized, values don't matter
        }
    }
}

void CoarserDefault::interpolate(std::shared_ptr<MGGrid> mggc, std::shared_ptr<MGGrid> mggf)
{
    for (int j = mggc->pJBegin(); j < mggc->pJEnd(); j++) // coarse grid indices
    {
        for (int i = mggc->pIBegin(); i < mggc->pIEnd(); i++)
        {
            int i_f = (i-mggc->pIBegin())*2 + mggf->pIBegin(); // fine grid indices (bottom left corner)
            int j_f = (j-mggc->pJBegin())*2 + mggf->pJBegin();
            // reusing resVec on fine grid to store corrector vector (interpolated solution from coarse grid)
            mggf->resVec(i_f,j_f) = mggc->resVec(i,j);
            mggf->resVec(i_f+1,j_f) = mggc->resVec(i,j);
            mggf->resVec(i_f,j_f+1) = mggc->resVec(i,j);
            mggf->resVec(i_f+1,j_f+1) = mggc->resVec(i,j);

            // mggf->p must not be changed
            // mggf->rhs must not be changed
        }
    }
};
