#include "Computation.h"

class ComputationParallel : public Computation 
{
  public:
	void initialize (int argc, char *argv[]);
  
  private:
	void computeTimeStepWidth ();
  
  
  
}