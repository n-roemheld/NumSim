#include "Coarser2.h"

void Coarser2::restrict(std::shared_ptr<MGGrid> mggf, std::shared_ptr<MGGrid> mggc) // pointer to fine and coarse grid  // make input const?
{
    for (int j = mggc->pJBegin(); j < mggc->pJEnd(); j++) // coarse grid indices
    {
        for (int i = mggc->pIBegin(); i < mggc->pIEnd(); i++)
        {
            int i_f = (i-mggc->pIBegin())*2 + mggf->pIBegin(); // fine grid indices (bottom left corner)
            int j_f = (j-mggc->pJBegin())*2 + mggf->pJBegin();
            mggc->rhs(i,j) = 1/16 * ( 4*mggf->resVec(i_f,j_f) 
            + 2 * ( mggf->resVec(i_f-1,j_f) + mggf->resVec(i_f+1,j_f) + mggf->resVec(i_f,j_f-1) + mggf->resVec(i_f,j_f+1) ) 
            + ( mggf->resVec(i_f-1,j_f-1) + mggf->resVec(i_f-1,j_f+1) + mggf->resVec(i_f+1,j_f-1) + mggf->resVec(i_f+1,j_f+1) )
            );
            // mggc->p must be zeros
            mggc->p(i,j) = 0;
            // mggc->resVec must be initialized, values don't matter
        }
    }
}

void Coarser2::interpolate(std::shared_ptr<MGGrid> mggc, std::shared_ptr<MGGrid> mggf)
{
    for (int j = mggc->pJBegin(); j < mggc->pJEnd(); j++) // coarse grid indices
    {
        for (int i = mggc->pIBegin(); i < mggc->pIEnd(); i++)
        {
            int i_f = (i-mggc->pIBegin())*2 + mggf->pIBegin(); // fine grid indices (bottom left corner)
            int j_f = (j-mggc->pJBegin())*2 + mggf->pJBegin();
            // reusing resVec on fine grid to store corrector vector (interpolated solution from coarse grid)
            mggf->resVec(i_f,j_f) = mggc->p(i,j);
            mggf->resVec(i_f+1,j_f) = .5 * (mggc->p(i,j) + mggc->p(i+1,j));
            mggf->resVec(i_f,j_f+1) = .5 * (mggc->p(i,j) + mggc->p(i,j+1));
            mggf->resVec(i_f+1,j_f+1) = .25 * (mggc->p(i,j) + mggc->p(i+1,j) + mggc->p(i,j+1) + mggc->p(i+1,j+1));

            // mggf->p must not be changed
            // mggf->rhs must not be changed
        }
    }
};
