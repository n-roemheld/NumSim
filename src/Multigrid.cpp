#include "Multigrid.h"
#include "GaussSeidel.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <math.h>

Multigrid:: Multigrid(std::shared_ptr<Discretization> discretization, std::shared_ptr<Smoother> sm, std::shared_ptr<Coarser> coa, std::shared_ptr<EndSolver> es, Cycle cycle, double epsilon, int maximumNumberOfIterations) :
    PressureSolver(discretization, epsilon, maximumNumberOfIterations), smoother_(sm), coarser_(coa), endSolver_(es), cycle_(cycle) // epsilon und maximumNumberOfIterations einfach auf 0 gesetzt, da nicht gebraucht
{};

void Multigrid::solve()
{
  // std::unique_ptr<GaussSeidel> pressureSolver_ = std::make_unique<GaussSeidel>(discretization_, epsilon_,
  //  maximumNumberOfIterations_);
  //  // discretization_->p(5,5) = -0.7;
  //
  //  pressureSolver_->solve();
   // discretization_->p(5,5) = -7;

    std::shared_ptr<MGGrid> mgg = std::make_shared<MGGrid>(discretization_->nCells(), discretization_->meshWidth(), discretization_->p(), discretization_->rhs());

    // for testing:
    for (int j = 0; j <= mgg->pJEnd(); j++)
    {
      for (int i = 0; i <= mgg->pIEnd(); i++)
      {
        mgg->p(i,j) = 0;
      }
    }

    if(cycle_.recursive)
    {
      int it = 0;
      double res_squared = 2*epsilon_*epsilon_;
      while(it < maximumNumberOfIterations_ && res_squared > epsilon_*epsilon_)
      {
        MGCycle(cycle_.maxLevel, mgg);
        res_squared = compute_res(mgg);
        it++;
      }
      // std::cout << it << std::endl;
    }
    else
    {
        // TODO: iterative
        MGLoop(cycle_.maxLevel, mgg);
    }

    for (int j = discretization_->pJBegin(); j < discretization_->pJEnd(); j++)
    {
      for (int i = discretization_->pIBegin(); i < discretization_->pIEnd(); i++)
      {
        discretization_->p(i,j) = mgg->p().operator()(i,j);
        discretization_->rhs(i,j) = mgg->rhs().operator()(i,j);
      }
    }

};

void Multigrid::MGCycle(int level, std::shared_ptr<MGGrid> mgg)
{
    // std::cout << "Level1: " << level << std::endl;
    if (mgg->nCells()[0] == 2 || mgg->nCells()[1] == 2)
    {
        
    }
    else if (level == 0)
    {
        endSolver_->solve(mgg);
    }
    else
    {
        // std::array< int, 2 > nCellsmgg = mgg->nCells();
        // std::cout << "MGG.nCells: " << nCellsmgg[0] << "," << nCellsmgg[1] << std::endl;
            // writeToConsole(mgg, "before pre");
        smoother_->presmooth(mgg);
            // writeToConsole(mgg, "after pre");
        computeResVec(mgg);
            // writeToConsole(mgg, "after compute resVec");
        std::shared_ptr<MGGrid> mggc = std::make_shared<MGGrid>(mgg->nCells(), mgg->meshWidth());
        coarser_->restrict(mgg, mggc); // mggCoarse is set complete, also p
            // writeToConsole(mggc, "after restrict");
        // std::array< int, 2 > nCellsmggc = mggc->nCells();
        // std::cout << "MGGc.nCells: " << nCellsmggc[0] << "," << nCellsmggc[1] << std::endl;
        for(int i = 0; i < cycle_.gamma.at(level); i++)
        {
            MGCycle(level-1, mggc);
        }
            // writeToConsole(mggc, "after recursive");
        coarser_->interpolate(mggc, mgg);
            // writeToConsole(mgg, "after interpolate");
        for(int j = mgg->pJBegin(); j < mgg->pJEnd(); j++)
        {
            for(int i = mgg->pIBegin(); i < mgg->pIEnd(); i++)
            {
                mgg->p(i,j) = mgg->p(i,j) +0.5 * mgg->resVec(i,j); // adding theta ???
            }
        }
            // writeToConsole(mgg, "after adding");
        smoother_->postsmooth(mgg);
            // writeToConsole(mgg, "after post");
    }
    //  std::cout << "Level2: " << level << std::endl;
};

void Multigrid::MGLoop(int maxLevel, std::shared_ptr<MGGrid> mgg)
{
    std::vector<std::shared_ptr<MGGrid>> grids;
    grids.push_back(mgg);
    for (int l = 0; l < maxLevel-1; l++)
    {
        smoother_->presmooth(grids.at(l));
        computeResVec(grids.at(l));
        std::shared_ptr<MGGrid> mggc = std::make_shared<MGGrid>(grids.at(l)->nCells(), grids.at(l)->meshWidth());
        coarser_->restrict(mgg, mggc); // mggCoarse is set complete, also p
        grids.push_back(mggc);
    }
    endSolver_->solve(mgg);
    for (int l = maxLevel-2; l >= 0; l--)
    {
        coarser_->interpolate(grids.at(l-1), grids.at(l));
        for(int j = grids.at(l)->pJBegin(); j < grids.at(l)->pJEnd(); j++)
        {
            for(int i = grids.at(l)->pIBegin(); i < grids.at(l)->pIEnd(); i++)
            {
                grids.at(l)->p(i,j) = grids.at(l)->p(i,j) + grids.at(l)->resVec(i,j); // adding theta ???
            }
        }
        smoother_->postsmooth(grids.at(l));
    }
}


void Multigrid::computeResVec(std::shared_ptr < MGGrid> mgg)
{
    std::array<double,2> mW = mgg->meshWidth();
	double dx = mW[0];
	double dy = mW[1];
    for(int j = mgg->pJBegin(); j < mgg->pJEnd(); j++)
    {
        for(int i = mgg->pIBegin(); i < mgg->pIEnd(); i++)
        {
            mgg->resVec(i,j) = mgg->rhs(i,j) - (mgg->p(i+1,j)-2*mgg->p(i,j)+mgg->p(i-1,j))/(dx*dx) - (mgg->p(i,j+1)-2*mgg->p(i,j)+mgg->p(i,j-1))/(dy*dy);
        }
    }
};

void Multigrid::setBoundaryValuesMGGrid(std::shared_ptr<MGGrid> mgg)
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

double Multigrid::compute_res(std::shared_ptr<MGGrid> mgg)
{
	//Array2D res_vec = Array2D(discretization_->nCells());
	double res = 0;
	FieldVariable p = mgg->p();
	std::array<double,2> mW = mgg->meshWidth();
	double dx = mW[0];
	double dy = mW[1];
	for(int j = mgg->pJBegin(); j < mgg->pJEnd(); j++)
	{
		for (int i = mgg->pIBegin(); i < mgg->pIEnd(); i++)
		{
			res += pow((p(i+1,j)-2*p(i,j)+p(i-1,j))/(dx*dx) + (p(i,j+1)-2*p(i,j)+p(i,j-1))/(dy*dy) - mgg->rhs(i,j),2);
		};
	};

	// average residuum with respect to number of cells
	res /= (mgg->nCells()[0]*mgg->nCells()[1]);
	return res;
};

void Multigrid::writePtoConsole(std::shared_ptr<MGGrid> mgg, std::string word)
{
    std::stringstream file;
    // write mesh width
    file << word << std::endl;
//   file << "nCells: " << mgg->nCells()[0] << "x" << mgg->nCells()[1]
//     << ", dx: " << mgg->meshWidth()[0] << ", dy: " << mgg->meshWidth()[1] << std::endl << std::endl;

  const int fieldWidth = 9;   // number of characters to use for a single value

  // write p
  // ---------
  // write header lines
  file << "p (" << mgg->p().size()[0] << "x" << mgg->p().size()[1] << "): " << std::endl
    << std::string(fieldWidth, ' ') << "|";
  for (int i = mgg->pIBegin()-1; i < mgg->pIEnd()+1; i++)
  {
    file << std::setw(fieldWidth) << i;
  }
  file << std::endl << std::string(fieldWidth*(mgg->p().size()[0]+2)+1, '-') << std::endl;

  // write p values
  for (int j = mgg->pJEnd(); j >= mgg->pJBegin()-1; j--)
  {
    file << std::setw(fieldWidth) << j << "|";
    for (int i = mgg->pIBegin()-1; i < mgg->pIEnd()+1; i++)
    {
      file << std::setw(fieldWidth) << std::setprecision(fieldWidth-6) << mgg->p(i,j);
    }
    file << std::endl;
  }
  file << std::endl;
  std::cout << file.str() << std::endl;
}

void Multigrid::writeRHStoConsole(std::shared_ptr<MGGrid> mgg, std::string word)
{
    std::stringstream file;
    // write mesh width
    file << word << std::endl;
//   file << "nCells: " << mgg->nCells()[0] << "x" << mgg->nCells()[1]
//     << ", dx: " << mgg->meshWidth()[0] << ", dy: " << mgg->meshWidth()[1] << std::endl << std::endl;

  const int fieldWidth = 9;   // number of characters to use for a single value

  // write rhs
  // ---------
  // write header lines
  file << "rhs (" << mgg->rhs().size()[0] << "x" << mgg->rhs().size()[1] << "): " << std::endl
    << std::string(fieldWidth, ' ') << "|";
  for (int i = mgg->pIBegin()-1; i < mgg->pIEnd()+1; i++)
  {
    file << std::setw(fieldWidth) << i;
  }
  file << std::endl << std::string(fieldWidth*(mgg->rhs().size()[0]+2)+1, '-') << std::endl;

  // write rhs values
  for (int j = mgg->pJEnd(); j >= mgg->pJBegin()-1; j--)
  {
    file << std::setw(fieldWidth) << j << "|";
    for (int i = mgg->pIBegin()-1; i < mgg->pIEnd()+1; i++)
    {
      file << std::setw(fieldWidth) << std::setprecision(fieldWidth-6) << mgg->rhs(i,j);
    }
    file << std::endl;
  }
  file << std::endl;
  std::cout << file.str() << std::endl;
}

void Multigrid::writeREStoConsole(std::shared_ptr<MGGrid> mgg, std::string word)
{
    std::stringstream file;
    // write mesh width
    file << word << std::endl;
//   file << "nCells: " << mgg->nCells()[0] << "x" << mgg->nCells()[1]
//     << ", dx: " << mgg->meshWidth()[0] << ", dy: " << mgg->meshWidth()[1] << std::endl << std::endl;

  const int fieldWidth = 9;   // number of characters to use for a single value

  // write resVec
  // ---------
  // write header lines
  file << "resVec (" << mgg->resVec().size()[0] << "x" << mgg->resVec().size()[1] << "): " << std::endl
    << std::string(fieldWidth, ' ') << "|";
  for (int i = mgg->pIBegin()-1; i < mgg->pIEnd()+1; i++)
  {
    file << std::setw(fieldWidth) << i;
  }
  file << std::endl << std::string(fieldWidth*(mgg->resVec().size()[0]+2)+1, '-') << std::endl;

  // write resVec values
  for (int j = mgg->pJEnd(); j >= mgg->pJBegin()-1; j--)
  {
    file << std::setw(fieldWidth) << j << "|";
    for (int i = mgg->pIBegin()-1; i < mgg->pIEnd()+1; i++)
    {
      file << std::setw(fieldWidth) << std::setprecision(fieldWidth-6) << mgg->resVec(i,j);
    }
    file << std::endl;
  }
  file << std::endl;
  std::cout << file.str() << std::endl;
}

void Multigrid::writeToConsole(std::shared_ptr<MGGrid> mgg, std::string word)
{
    writePtoConsole(mgg, word);
    writeRHStoConsole(mgg, "");
    writeREStoConsole(mgg, "");

}
