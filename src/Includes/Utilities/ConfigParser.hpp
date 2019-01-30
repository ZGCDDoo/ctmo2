//
// Created by charles-david on 29/01/19.
//

#ifndef CTMO_CONFIGPARSER_HPP
#define CTMO_CONFIGPARSER_HPP

#include <snappy.h>
#include "Utilities.hpp"
#include "VerticesSimple.hpp"

namespace Diagrammatic
{
    class ConfigParser
    {
    public:

        const size_t NUM_MAX_CONFIGS = 50;
        const std::string NAME_PREFIX_BATCH = "configsBatch";

        explicit ConfigParser() : configs_()
        {}

        void SaveConfig(const Vertices &vertices_)
        {
            const size_t ss = vertices_.size();
            std::vector<double> auxSpins(ss);
            std::vector<double> taus(ss);

            for (size_t ii = 0; ii < ss; ++ii)
            {
                const auto vPartUp = vertices_.atUp(ii);
                auxSpins.at(ii) = (vPartUp.aux() == AuxSpin_t::Up) ? 1.0 : -1.0;
                taus.at(ii) = vPartUp.tau();
            }


            configs_[std::to_string(configId_)]["taus"] = taus;
            configs_[std::to_string(configId_)]["auxSpins"] = auxSpins;

            ++configId_;
            if (configId_ == NUM_MAX_CONFIGS)
            {
                const std::string fileName = NAME_PREFIX_BATCH + std::to_string(batchId_) + ".snappy";
                std::ofstream fout(fileName);
                fout << configs_.dump();
                fout.close();

                std::cout <<"Not using snappy yet, use snappy compress" << std::endl;

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
}

#endif //CTMO_CONFIGPARSER_HPP
