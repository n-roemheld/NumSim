#include "Smoother.h"

Smoother::Smoother(int numberOfIterationsPre, int numberOfIterationsPost) :
    numberOfIterationsPre_(numberOfIterationsPre), numberOfIterationsPost_(numberOfIterationsPost)
{};
    
void Smoother::setBoundaryValues(std::shared_ptr<MGGrid> mgg)
{
    // lower p ghost layer
		int j = mgg->pJBegin()-1;
		for(int i = mgg->pIBegin(); i < mgg->pIEnd(); i++)
		{
			mgg->p(i,j) = mgg->p(i,j+1);
		};
		// upper p ghost layer
		j = mgg->pJEnd();
		for(int i = mgg->pIBegin(); i < mgg->pIEnd(); i++)
		{
			mgg->p(i,j) = mgg->p(i,j-1);
		};
		// left p ghost layer
		int i = mgg->pIBegin()-1;
		for(int j = mgg->pJBegin()-1; j < mgg->pJEnd()+1; j++)
		{
			mgg->p(i,j) = mgg->p(i+1,j);
		}
		// right p ghost layer
		i = mgg->pIEnd();
		for(int j = mgg->pJBegin()-1; j < mgg->pJEnd()+1; j++)
		{
			mgg->p(i,j) = mgg->p(i-1,j);
		}
};
