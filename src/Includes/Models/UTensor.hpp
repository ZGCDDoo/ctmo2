#pragma once

#include "../Utilities/Utilities.hpp"

namespace Models
{

class UTensor
{

  public:
    UTensor(const Json &UJson) : NOrb_(tJson["NOrb"].get<size_t>())

    {
        assert(UJson["UParameters"].size() == NOrb_ * (NOrb_ + 1) / 2);

        for (size_t o1 = 0; o1 < NOrb_; o1++)
        {
            for (size_t o2 = o1; o2 < NOrb_; o2++)
            {

                const std::string o1o2Name = std::to_string(o1) + std::to_string(o2);
                std::cout << "o1o2Name = " << o1o2Name << std::endl;
                const Json &jj = UJson["UParameters"][o1o2Name];

                const double U = jj["U"].get<double>();
                const double J_H = jj["J_H"].get<double>();

                UVec_.push_back(U);
                JHVec_.push_back(J_H);
                UPrimeVec_.push_back(U - 2.0 * J_H);
            }
        }
    }

    ~UTensor()
    {
    }

    std::vector<double> UVec() const { return UVec_; };
    std::vector<double> UPrimeVec() const { return UPrimeVec_; };
    std::vector<double> JHVec() const { return JHVec_; };

  protected:
    std::vector<double> UVec_;
    std::vector<double> UPrimeVec_;
    std::vector<double> JHVec_;
    const size_t NOrb_;
};

} // namespace Models
