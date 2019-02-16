#pragma once

namespace MC
{

class ABC_MonteCarlo
{
  public:
    ABC_MonteCarlo() = default;
    ABC_MonteCarlo(const ABC_MonteCarlo &abc_monteCarlo) = default;
    ABC_MonteCarlo(ABC_MonteCarlo &&abc_monteCarlo) = default;

    ABC_MonteCarlo &operator=(const ABC_MonteCarlo &abc_monteCarlo) = delete;
    ABC_MonteCarlo &operator=(ABC_MonteCarlo &&abc_monteCarlo) = delete;

    virtual ~ABC_MonteCarlo() = 0;
    virtual void RunMonteCarlo() = 0;
}; // class ABC_MonteCarlo

ABC_MonteCarlo::~ABC_MonteCarlo() = default; // destructors must exist

} // namespace MC
