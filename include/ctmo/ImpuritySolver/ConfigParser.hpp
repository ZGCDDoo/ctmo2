//
// Created by charles-david on 29/01/19.
//

#ifndef CTMO_CONFIGPARSER_HPP
#define CTMO_CONFIGPARSER_HPP

#include <snappy.h>
#include "ctmo/ImpuritySolver/VerticesSimple.hpp"
#include "ctmo/deps/nlohmann_json/json.hpp"

namespace Diagrammatic
{
class ConfigParser
{
public:
    const size_t NUM_MAX_CONFIGS = 500;
    const std::string NAME_PREFIX = "configs_";
    const std::string NAME_SUFFIX = "_ctmo.bin";

    explicit ConfigParser(const size_t &rank) : configs_(), rank_{rank}
    {}

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
//            if(tries_ != 0)
//            {
//                return;
//            }
            const std::vector<std::uint8_t> msgpackConfigs = nlohmann::json::to_msgpack(configs_);
            const std::string msgpackConfigsStr(msgpackConfigs.begin(), msgpackConfigs.end());
            std::string msgpackConfigsStrCompressed;
            snappy::Compress(msgpackConfigsStr.data(), msgpackConfigsStr.size(), &msgpackConfigsStrCompressed);
            const std::string fileName = NAME_PREFIX + std::to_string(rank_) + NAME_SUFFIX;

            const std::uint32_t sizeStr = msgpackConfigsStrCompressed.size();
            std::cout << "sizeStr = " << sizeStr << std::endl;
            std::vector<char> dataBin;
            dataBin.insert(dataBin.end(), reinterpret_cast<const char *>(&sizeStr),
                           reinterpret_cast<const char *>(&sizeStr) + sizeof(std::uint32_t));
            dataBin.insert(dataBin.end(), msgpackConfigsStrCompressed.begin(), msgpackConfigsStrCompressed.end());

            std::ofstream fout(fileName, std::ios::out | std::ios::app | std::ios::binary | std::ios::ate);
            fout.write(dataBin.data(), dataBin.size());
            fout.close();

            configs_.clear();
            configId_ = 0;
//            tries_++;
        }
    }


private:
    Json configs_;
    size_t configId_{0};
    const size_t rank_;
//    int tries_{0};
};
} // namespace Diagrammatic

#endif // CTMO_CONFIGPARSER_HPP
