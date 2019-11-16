#include "Computation.h"

class ComputationParallel : public Computation
{
  public:
	void initialize (int argc, char *argv[]);

  private:
	void computeTimeStepWidth ();

	void computePreliminaryVelocities();

  void applyBoundaryValues();

  void computeVelocities();

  void velocity_communication();



};
