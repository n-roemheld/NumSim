#pragma once

#include <array>
#include <memory>

#include "FieldVariable.h"

class MGGrid
{
public:
    MGGrid(std::array< int, 2> nCells, std::array< double, 2> meshWidth, std::shared_ptr<FieldVariable> p, std::shared_ptr<FieldVariable> rhs);
    MGGrid(std::array< int, 2> nCells, std::array< double, 2> meshWidth);
    
    const std::array< double, 2 > 	meshWidth() const;

    const std::array< int, 2 > 	nCells() const;

    const FieldVariable & 	p() const;

    const FieldVariable &   rhs() const;

    const FieldVariable &   resVec() const;

    double 	p(int i, int j) const;

    double& 	p(int i, int j);

    double& 	rhs(int i, int j);

    double&     resVec(int i, int j);

    int 	pIBegin() const;

    int 	pIEnd() const;

    int 	pJBegin() const;

    int 	pJEnd() const;

protected:
    const std::array< int, 2 > nCells_;
    const std::array< double, 2 > meshWidth_;
    FieldVariable p_;
    FieldVariable rhs_;
    FieldVariable resVec_;
};