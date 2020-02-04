#pragma once

#include "Coarser.h"

class CoarserLinear : public Coarser
{
public:
  CoarserLinear() : Coarser()
  {};

  // restricts the current MGGrid to the coarser MGGrid and sets also nCells and meshWidth
  void restrict(std::shared_ptr<MGGrid> mggf, std::shared_ptr<MGGrid> mggc);

  // interpolates the coarse MGGrid to the finer MGGrid in p
  void interpolate(std::shared_ptr<MGGrid> mggc, std::shared_ptr<MGGrid> mggf);
};