#pragma once

#include "../Utilities/Utilities.hpp"

namespace Models
{

class UTensor
{

  public:
    UTensor(const Json &UJson, const size_t &Nc) : NOrb_(UJson["NOrb"].get<size_t>()), Nc_(Nc)

    {
        std::cout << "Start of UTensor constructor" << std::endl;
        assert(UJson["UParameters"].size() == NOrb_ * (NOrb_ + 1) / 2);

        auxMu_ = ClusterMatrixCD_t(NOrb_, NOrb_).zeros();
        auxMu_ += UJson["mu"].get<double>() * ClusterMatrixCD_t(NOrb_, NOrb_).eye();

        for (size_t o1 = 0; o1 < NOrb_; o1++)
        {
            for (size_t o2 = o1; o2 < NOrb_; o2++)
            {

                const std::string o1o2Name = std::to_string(o1) + std::to_string(o2);
                std::cout << "o1o2Name = " << o1o2Name << std::endl;
                const Json &jj = UJson["UParameters"][o1o2Name];

                const double U = jj["U"].get<double>();

                double J_H = 0.0;
                double UPrime = 0.0;
                if (o1 != o2)
                {
                    J_H = jj["J_H"].get<double>();
                    UPrime = U - 2.0 * J_H;
                }

                UVec_.push_back(U);
                JHVec_.push_back(J_H);
                UPrimeVec_.push_back(UPrime);
                // auxMu_
            }
        }
        std::cout << "End of UTensor constructor" << std::endl;
    }

    ~UTensor()
    {
    }

    std::vector<double> UVec() const { return UVec_; };
    std::vector<double> UPrimeVec() const { return UPrimeVec_; };
    std::vector<double> JHVec() const { return JHVec_; };

    ClusterMatrixCD_t auxMu() const { return auxMu_; };

  protected:
    std::vector<double> UVec_;
    std::vector<double> UPrimeVec_;
    std::vector<double> JHVec_;
    ClusterMatrixCD_t auxMu_;
    const size_t NOrb_;
    const size_t Nc_;
};

} // namespace Models
