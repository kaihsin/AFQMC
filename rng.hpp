#ifndef RNG_HPP
#define RNG_HPP
#include "mersennetwister.h"
class RNG:public MtRng64
{
  public:
    RNG(unsigned long seed) { MtRng64::init(seed); }
    double operator () () { return MtRng64::getReal2();};
    ~RNG(){ };
};

#endif
