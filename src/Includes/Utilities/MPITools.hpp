#pragma once

#include "Utilities.hpp"

#ifdef HAVEMPI
#include <boost/mpi.hpp>
namespace mpi = boost::mpi;
#endif

namespace mpiUt
{

class Tools
{
  public:
    static const int master = 0;

    static int NWorkers()
    {
#ifdef HAVEMPI
        mpi::communicator world;
        return world.size();
#endif
#ifndef HAVEMPI
        return 1;
#endif
    }

    static int Rank()
    {
#ifdef HAVEMPI
        mpi::communicator world;
        return world.rank();
#endif
#ifndef HAVEMPI
        return 0;
#endif
    }

    static void Print(const std::string &message)
    {
#ifdef HAVEMPI
        mpi::communicator world;
        if (world.rank() == 0)
        {
            std::cout << "Rank " << std::to_string(Rank()) << ": " << message << std::endl;
        }
#endif

#ifndef HAVEMPI
        std::cout << message << std::endl;
#endif
    }

    static std::string SaveUpdStats(const std::vector<UpdStats_t> &updStatsVec)
    {

        std::vector<std::string> keys;
        UpdStats_t updStatsResult(updStatsVec.at(0));

        for (UpdStats_t::iterator it = updStatsResult.begin(); it != updStatsResult.end(); ++it)
        {
            keys.push_back(it->first);
        }

        for (size_t i = 1; i < updStatsVec.size(); i++)
        {
            for (const std::string &key : keys)
            {
                const std::valarray<size_t> tmp = updStatsVec.at(i).at(key);
                updStatsResult[key] += tmp;
            }
        }

        Json jjout(updStatsResult);
        for (const std::string key : keys)
        {

            updStatsResult[key][0] = double(updStatsResult[key][0]) / double(NWorkers());
            updStatsResult[key][1] = double(updStatsResult[key][1]) / double(NWorkers());
            size_t nbProposed = updStatsResult[key][0];
            size_t nbAccepted = updStatsResult[key][1];

            jjout[key] = {{"Proposed", nbProposed}, {"Accepted", nbAccepted}};
        }

        return jjout.dump(4);
    }

    static std::vector<cd_t> CubeCDToVecCD(const ClusterCubeCD_t &cubeCD)
    {
        // Print("start CubeCDToVecCD");

        std::vector<cd_t> vecCD;
        const size_t nrows = cubeCD.n_rows;
        const size_t ncols = cubeCD.n_cols;
        const size_t nslices = cubeCD.n_slices;
        vecCD.resize(cubeCD.n_elem, cd_t(0.0));

        for (size_t ii = 0; ii < nrows; ii++)
        {
            for (size_t jj = 0; jj < ncols; jj++)
            {
                for (size_t kk = 0; kk < nslices; kk++)
                {
                    const size_t index = ii + jj * nrows + (nrows * ncols) * kk;
                    vecCD.at(index) = cubeCD(ii, jj, kk);
                }
            }
        }

        // Print("End CubeCDToVecCD");
        return vecCD;
    }

    static ClusterCubeCD_t VecCDToCubeCD(std::vector<cd_t> &vecCD, const size_t &nrows, const size_t &ncols, const size_t &nslices)
    {
        // Print("start VecCDToCubeCD");

        ClusterCubeCD_t cubeCD(nrows, ncols, nslices);
        cubeCD.zeros();

        for (size_t ii = 0; ii < nrows; ii++)
        {
            for (size_t jj = 0; jj < ncols; jj++)
            {
                for (size_t kk = 0; kk < nslices; kk++)
                {
                    const size_t index = ii + jj * nrows + (nrows * ncols) * kk;
                    cubeCD(ii, jj, kk) = vecCD.at(index);
                }
            }
        }

        // Print("End VecCDToCubeCD");

        return cubeCD;
    }
};

} // namespace mpiUt