#include "SmootherJacobi.h"

SmootherJacobi::SmootherJacobi(int numberOfIterationsPre, int numberOfIterationsPost) :
    Smoother(numberOfIterationsPre, numberOfIterationsPost)
{};

void SmootherJacobi::presmooth(std::shared_ptr<MGGrid> mgg)
{
    smooth(mgg, numberOfIterationsPre_);
    // std::array<double,2> mW = mgg->meshWidth();
	// double dx = mW[0];
	// double dy = mW[1];

    // FieldVariable p_old = mgg->p(); // copy by value?

    // for(int it = 0; it < numberOfIterationsPre_; it++)
    // {
    //     setBoundaryValues(mgg);
    //     for(int j = mgg->pJBegin(); j < mgg->pJEnd(); j++)
    //     {
    //         for(int i = mgg->pJBegin(); i < mgg->pJEnd(); i++)
    //         {
    //             mgg->p(i,j) = (dx*dx*dy*dy)/(2*(dx*dx+dy*dy))
	// 			* ( (p_old(i-1,j) + p_old(i+1,j)) / (dx*dx)
	// 			+ (p_old(i,j-1) + p_old(i,j+1)) / (dy*dy)
	// 			- mgg->rhs(i,j) );
    //         }
    //     }
    // }
};

void SmootherJacobi::postsmooth(std::shared_ptr<MGGrid> mgg)
{
    smooth(mgg, numberOfIterationsPost_);
//    std::array<double,2> mW = mgg->meshWidth();
// 	double dx = mW[0];
// 	double dy = mW[1];

//     FieldVariable p_old = mgg->p(); // copy by reference (?)

//     for(int it = 0; it < numberOfIterationsPost_; it++)
//     {
//         setBoundaryValues(mgg);
//         for(int j = mgg->pJBegin(); j < mgg->pJEnd(); j++)
//         {
//             for(int i = mgg->pJBegin(); i < mgg->pJEnd(); i++)
//             {
//                 mgg->p(i,j) = (dx*dx*dy*dy)/(2*(dx*dx+dy*dy))
// 				* ( (p_old(i-1,j) + p_old(i+1,j)) / (dx*dx)
// 				+ (p_old(i,j-1) + p_old(i,j+1)) / (dy*dy)
// 				- mgg->rhs(i,j) );
//             }
//         }
//     }
};

void SmootherJacobi::smooth(std::shared_ptr<MGGrid> mgg, int numberOfIterations)
{
    std::array<double,2> mW = mgg->meshWidth();
	double dx = mW[0];
	double dy = mW[1];

    FieldVariable p_old = mgg->p(); // copy by reference (?)

    for(int it = 0; it < numberOfIterations; it++)
    {
        setBoundaryValues(mgg);
        for(int j = mgg->pJBegin(); j < mgg->pJEnd(); j++)
        {
            for(int i = mgg->pJBegin(); i < mgg->pJEnd(); i++)
            {
        //         mgg->p(i,j) = (dx*dx*dy*dy)/(2*(dx*dx+dy*dy))
				// * ( (p_old(i-1,j) + p_old(i+1,j)) / (dx*dx)
				// + (p_old(i,j-1) + p_old(i,j+1)) / (dy*dy)
				// - mgg->rhs(i,j) );
// GS
        mgg->p(i,j) = (dx*dx*dy*dy)/(2*(dx*dx+dy*dy))
* ( (mgg->p(i-1,j) + mgg->p(i+1,j)) / (dx*dx)
+ (mgg->p(i,j-1) + mgg->p(i,j+1)) / (dy*dy)
- mgg->rhs(i,j) );
            }
        }
    }
    setBoundaryValues(mgg);

};
