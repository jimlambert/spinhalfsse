#ifndef RANDOM_HELPERS_H
#define RANDOM_HELPERS_H

#include <random>

#include "lattice.h"

namespace spinhalfsse
{
  class RandomHelpers
  {
    private:

      std::random_device _rd;
      std::mt19937 _mteng{_rd()};
      std::uniform_real_distribution<> _rdist{0.0, 0.99999};
      std::uniform_int_distribution<> _rleg{0,3};

    public:

      double random() { return _rdist(_mteng); }
      int randleg() { return _rleg(_mteng); }
  };
}

#endif // RANDOM_HELPERS_H

