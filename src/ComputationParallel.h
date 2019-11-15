#include "Computation.h"

class ComputationParallel : public Computation 
{
  public:
	void initialize (int argc, char *argv[]);
  
  private:
	void computeTimeStepWidth ();
	void computePreliminaryVelocities ();
    void send_boundary_vertical(int direction, int j_fixed, int i_begin, int i_end, int target_rank, double (*fVar)(int, int), bool do_nothing);
    void send_boundary_horizontal(int direction, int i_fixed, int j_begin, int j_end, int target_rank, double (*fVar)(int, int), bool do_nothing);

  
  
  
}