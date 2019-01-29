
#include "Logging.hpp"

namespace IO
{

using GreenSites_t = std::vector<std::vector<std::pair<size_t, size_t>>>;

// Return the greensites, given a file name for the model at hand,
// in the file, the greensites are given as a complex matrix of ints.
GreenSites_t BuildGreenSites(const std::string &fname)
{
    Logging::Trace("Start of IO::BuildGreenSites");
    GreenSites_t greenSites;
    ClusterMatrixCD_t inGreen;
    inGreen.load(fname);

    for (size_t i = 0; i < inGreen.n_rows; ++i)
    {
        std::vector<std::pair<size_t, size_t>> tmpvec;
        for (size_t j = 0; j < inGreen.n_rows; ++j)
        {
            tmpvec.emplace_back(inGreen(i, j).real(), inGreen(i, j).imag());
        }

        greenSites.push_back(tmpvec);
    }

    Logging::Trace("greensites = ");
    for (size_t i = 0; i < inGreen.n_rows; ++i)
    {
        for (size_t j = 0; j < inGreen.n_rows; ++j)
        {
            Logging::Trace(std::to_string(greenSites.at(i).at(j).first) + ", " + std::to_string(greenSites.at(i).at(j).second));
        }
    }

    Logging::Trace("End of IO::BuildGreenSites. ");
    return greenSites;
}

std::vector<std::pair<size_t, size_t>> BuildIndepSites(const GreenSites_t &greenSites)
{
    std::vector<std::pair<size_t, size_t>> indepSites;

    for (size_t i = 0; i < greenSites.size(); ++i)
    {
        for (size_t j = 0; j < greenSites.size(); ++j)
        {
            auto itsite = std::find(indepSites.begin(), indepSites.end(), greenSites.at(i).at(j));
            if (itsite == indepSites.end())
            {
                indepSites.emplace_back(greenSites.at(i).at(j));
            }
        }
    }

    return indepSites;
}

} // namespace IO