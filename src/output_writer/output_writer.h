#pragma once

#include "Discretization.h"
#include "Partitioning.h"

#include <memory>

/** Inteface class for writing simulation data output.
 */
class OutputWriter
{
public:
  //! constructor
  //! @param discretization shared pointer to the discretization object that will contain all the data to be written to the file
  OutputWriter(std::shared_ptr<Discretization> discretization, const Partitioning &parti);

  //! write current velocities to file, filename is output_<count>.vti
  virtual void writeFile(double currentTime) = 0;

protected:

  std::shared_ptr<Discretization> discretization_;  //< a shared pointer to the discretization which contains all data that will be written to the file
  // const Partitioning partitioning_;                 //< the partitioning object that knowns about the domain decomposition, only significant when executing in parallel
  Partitioning partitioning_;                 //< the partitioning object that knowns about the domain decomposition, only significant when executing in parallel

  int fileNo_;   //< a counter that increments for every file, this number is part of the file name of output files
};
