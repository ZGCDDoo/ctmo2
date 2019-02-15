//
// Created by charles-david on 29/01/19.
//

#ifndef CTMO_CONFIGPARSER_HPP
#define CTMO_CONFIGPARSER_HPP

#include <snappy.h>
#include "ctmo/ImpuritySolver/VerticesSimple.hpp"

namespace Diagrammatic
{
class ConfigParser
{
  public:
    const size_t NUM_MAX_CONFIGS = 5000;
    const std::string NAME_PREFIX_BATCH = "configsBatch";

    explicit ConfigParser() : configs_() {}

    void SaveConfig(const Vertices &vertices_, const double &logDeterminant, const Sign_t &sign)
    {
        const size_t ss = vertices_.size();
        std::vector<int> auxSpins(ss);
        std::vector<double> taus(ss);

        for (size_t ii = 0; ii < ss; ++ii)
        {
            const auto vPartUp = vertices_.atUp(ii);
            auxSpins.at(ii) = (vPartUp.aux() == AuxSpin_t::Up) ? 1 : -1;
            taus.at(ii) = vPartUp.tau();
        }

        configs_[std::to_string(configId_)]["taus"] = taus;
        configs_[std::to_string(configId_)]["auxSpins"] = auxSpins;
        configs_[std::to_string(configId_)]["logDeterminant"] = logDeterminant;
        configs_[std::to_string(configId_)]["sign"] = sign;

        ++configId_;
        if (configId_ == NUM_MAX_CONFIGS)
        {
            const std::string batchStr = configs_.dump();
            std::string batchStrCompressed;
            snappy::Compress(batchStr.data(), batchStr.size(), &batchStrCompressed);
            const std::string fileName = NAME_PREFIX_BATCH + std::to_string(batchId_) + ".snappy";
            std::ofstream fout(fileName);
            fout << batchStrCompressed;
            fout.close();

            configs_.clear();
            configId_ = 0;
            ++batchId_;
        }
    }

  private:
    Json configs_;
    size_t configId_{0};
    size_t batchId_{0};
};
} // namespace Diagrammatic

#endif // CTMO_CONFIGPARSER_HPP
