#include "SmootherGaussSeidelUnique.h"

SmootherGaussSeidelUnique::SmootherGaussSeidelUnique(int numberOfIterationsPre, int numberOfIterationsPost) :
    Smoother(numberOfIterationsPre, numberOfIterationsPost)
{};

void SmootherGaussSeidelUnique::presmooth(std::shared_ptr<MGGrid> mgg)
{
    smooth(mgg, numberOfIterationsPre_); 
};

void SmootherGaussSeidelUnique::postsmooth(std::shared_ptr<MGGrid> mgg)
{
    smooth(mgg, numberOfIterationsPost_);
};

void SmootherGaussSeidelUnique::smooth(std::shared_ptr<MGGrid> mgg, int numberOfIterations)
{
     std::array<double,2> mW = mgg->meshWidth();
	double dx = mW[0];
	double dy = mW[1];

    for(int it = 0; it < numberOfIterations; it++)
    {
        setBoundaryValues(mgg);
        for(int j = mgg->pJBegin(); j < mgg->pJEnd(); j++)
        {
            for(int i = mgg->pIBegin(); i < mgg->pIEnd(); i++)
            {
                if (i == mgg->pIBegin() && j == mgg->pJBegin())
                {
                    mgg->p(i,j) = 0;
                }
                else
                {
                    mgg->p(i,j) = (dx*dx*dy*dy)/(2*(dx*dx+dy*dy))
                    * ( (mgg->p(i-1,j) + mgg->p(i+1,j)) / (dx*dx)
                    + (mgg->p(i,j-1) + mgg->p(i,j+1)) / (dy*dy)
                    - mgg->rhs(i,j) );
                }
                
            }
        }
    }
    setBoundaryValues(mgg);

};
