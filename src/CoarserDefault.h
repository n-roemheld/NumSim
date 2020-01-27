#include "Coarser.h"

class CoarserDefault : public Coarser
{

    public:
    CoarserDefault();

    // restricts the current MGGrid to the coarser MGGrid and sets also nCells and meshWidth 
    std::shared_ptr<MGGrid> restrict(std::shared_ptr<MGGrid> mgg);

    // interpolates the coarse MGGrid to the finer MGGrid in p
    std::shared_ptr<FieldVariable> interpolate(std::shared_ptr<MGGrid> mggCoarse);
};