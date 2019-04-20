//
// Created by charles-david on 29/01/19.
//

#ifndef CTMO_CONFIGPARSER_HPP
#define CTMO_CONFIGPARSER_HPP

#include <snappy.h>
#include "ctmo/ImpuritySolver/VerticesSimple.hpp"
#include "ctmo/deps/nlohmann_json/json.hpp"
#include "ctmo/Foundations/Logging.hpp"

#ifdef HAVE_POSTGRES
#include "ctmo/Foundations/PgDriver.hpp"
#endif

#include <csignal>

namespace Diagrammatic
{
class ConfigParser
{
public:
    const std::string NAME_PREFIX = "configs_";
    const std::string NAME_SUFFIX = "_ctmo.bin";
    const size_t NUM_MAX_CONFIGS = 500;
    const size_t NUM_MAX_BATCH_SAVED = 10;

    const static int SIMULATION_ID_DEFAULT = -9999;

    explicit ConfigParser(const size_t &rank, const int& simulation_id=SIMULATION_ID_DEFAULT)
            : configs_(), rank_{rank}, simulationId_{simulation_id}
    {}

    int batchesSaved() const
    { return batchesSaved_; }

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
//            if(batchesSaved_ == 21)
//            {   if(mpiUt::Tools::Rank() == mpiUt::Tools::master)
//                {
//                    const std::string jjStr = configs_.dump();
//                    std::ofstream foutTmp("batchConfig21.json");
//                    foutTmp << jjStr << std::endl;
//                    foutTmp.close();
//                }
//
//            }
            const std::vector<std::uint8_t> msgpackConfigs = nlohmann::json::to_msgpack(configs_);
            const std::string msgpackConfigsStr(msgpackConfigs.begin(), msgpackConfigs.end());
            std::string msgpackConfigsStrCompressed;
            snappy::Compress(msgpackConfigsStr.data(), msgpackConfigsStr.size(), &msgpackConfigsStrCompressed);
            const std::string fileName = NAME_PREFIX + std::to_string(rank_) + NAME_SUFFIX;

            const std::uint32_t sizeStr = msgpackConfigsStrCompressed.size();
            std::vector<char> dataBin;
            dataBin.insert(dataBin.end(), reinterpret_cast<const char *>(&sizeStr),
                           reinterpret_cast<const char *>(&sizeStr) + sizeof(std::uint32_t));
            dataBin.insert(dataBin.end(), msgpackConfigsStrCompressed.begin(), msgpackConfigsStrCompressed.end());

            std::ofstream fout(fileName, std::ios::out | std::ios::app | std::ios::binary | std::ios::ate);
            fout.write(dataBin.data(), dataBin.size());
            fout.close();

            configs_.clear();
            configId_ = 0;
            batchesSaved_++;
            if (batchesSaved_ > NUM_MAX_BATCH_SAVED)
            {
                Logging::Info("In SLMC, NUM_BATCH_SAVED attained, exiting");
                raise(SIGINT);
            }

#ifdef HAVE_POSTGRES
            DB::PgDriver* pgDriverPtr = DB::PgDriver::getInstance();
            pgDriverPtr->SaveBytea(msgpackConfigsStrCompressed, simulationId_);
#endif
        }
    }


private:
    Json configs_;
    size_t configId_{0};
    const size_t rank_;
    const int simulationId_;
    size_t batchesSaved_{0};

};
} // namespace Diagrammatic

#endif // CTMO_CONFIGPARSER_HPP
