#include "Computation.h"
#include "output_writer/output_writer_paraview_parallel.h"
#include "output_writer/output_writer_text_parallel.h"


class ComputationParallel : public Computation
{
  public:
	void initialize (int argc, char *argv[]);
  void runSimulation ();

  private:
	void computeTimeStepWidth ();

	void computePreliminaryVelocities();

  void applyBoundaryValues();

  void computeVelocities();

  void preliminaryVelocity_communication();
  void finalVelocity_communication();

  std::unique_ptr< OutputWriterParaviewParallel > outputWriterParaviewParallel_;
	std::unique_ptr< OutputWriterTextParallel > outputWriterTextParallel_;

  
	void gdbParallelDebuggingBarrier();

};
