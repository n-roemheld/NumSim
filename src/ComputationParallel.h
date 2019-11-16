#include "Computation.h"

class ComputationParallel : public Computation 
{
  public:
	void initialize (int argc, char *argv[]);
  
  private:
	void computeTimeStepWidth ();
    
	void computePreliminaryVelocities ();
    
    void velocity_communication(double (*u_or_f)(int, int), double (*v_or_g)(int, int));

  
  
}