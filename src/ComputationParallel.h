#include "Computation.h"
#include "output_writer/output_writer_paraview_parallel.h"
#include "output_writer/output_writer_text_parallel.h"

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
