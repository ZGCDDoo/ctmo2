#pragma once

#include <boost/random.hpp>

namespace Utilities
{

    using EngineTypeMt19937_t = boost::mt19937;
    using EngineTypeFibonacci3217_t = boost::lagged_fibonacci3217;
    using UniformDistribution_t = boost::uniform_real<double>;
    using UniformRngMt19937_t = boost::variate_generator<EngineTypeMt19937_t &, UniformDistribution_t>;
    using UniformRngFibonacci3217_t = boost::variate_generator<EngineTypeFibonacci3217_t &, UniformDistribution_t>;

    

} // namespace Utilities
