#pragma once

#include "IO.hpp"
#include "MPITools.hpp"
#include "../IS/ISResult.hpp"

#ifdef HAVEMPI
#include <boost/mpi.hpp>
namespace mpi = boost::mpi;
#endif

namespace mpiUt
{

class IOResult
{
  public:
    static void SaveISResults(const std::vector<Result::ISResult> &isResultVec, const IO::Base_IOModel &ioModel, const double &beta)
    {
        const int nworkers = Tools::NWorkers();
        assert(nworkers == int(isResultVec.size()));
        const size_t n_cols = isResultVec.at(0).n_cols_;
        const size_t n_rows = isResultVec.at(0).n_rows_;
        const size_t fillingSize = ioModel.fillingSites().size();
        std::valarray<cd_t> greenTabResultUp(cd_t(0.0, 0.0), n_rows * n_cols);
        std::valarray<cd_t> greenTabResultDown(cd_t(0.0, 0.0), n_rows * n_cols);
        std::valarray<double> fillingResultUp(0.0, fillingSize);
        std::valarray<double> fillingResultDown(0.0, fillingSize);

        //Average the greens of matsubara, and the fillingsSigma_.
        for (int i = 0; i < nworkers; i++)
        {
            greenTabResultUp += isResultVec.at(i).greenTabUp_;

#ifdef AFM
            greenTabResultDown += isResultVec.at(i).greenTabDown_;
#endif
            fillingResultUp += isResultVec.at(i).fillingUp_;
            fillingResultDown += isResultVec.at(i).fillingDown_;
        }

        greenTabResultUp /= double(nworkers);
        greenTabResultDown /= double(nworkers);
        fillingResultUp /= double(nworkers);
        fillingResultDown /= double(nworkers);
        //convert the greens to ClusterMatrixCD_t
        ClusterMatrixCD_t greenUp(n_rows, n_cols);
        ClusterMatrixCD_t greenDown(n_rows, n_cols);
        for (size_t j = 0; j < n_cols; j++)
        {
            for (size_t i = 0; i < n_rows; i++)
            {
                greenUp(i, j) = greenTabResultUp[i + n_rows * j];
#ifdef AFM
                greenDown(i, j) = greenTabResultDown[i + n_rows * j];
#endif
            }
        }
        //greenUp.print();
        //greenDown.print();
        const size_t PRECISION_OUT = 14;

        SaveTabular("greenUp", greenUp, beta, PRECISION_OUT, false);

        SaveTabular("greenDown", greenDown, beta, PRECISION_OUT, false);

        //Average the obsScale_
        SaveFillingMatrixs(fillingResultUp, fillingResultDown, ioModel);
        StatsJsons(isResultVec);
    }

    static void StatsJsons(const std::vector<Result::ISResult> &isResultVec)
    {
        Json jjResult(isResultVec.at(0).obsScal_);

        const size_t nworkers = Tools::NWorkers();
        const size_t jjSize = jjResult.size();

        //START STATS===================================================================
        ClusterMatrix_t statsMat(nworkers, jjSize); //Each row contains the data for each jj
        std::vector<std::string> keys;

        for (Json::iterator it = jjResult.begin(); it != jjResult.end(); ++it)
        {
            keys.push_back(it.key());
        }

        for (size_t i = 0; i < nworkers; i++)
        {
            std::map<std::string, double> tmpMap = isResultVec.at(i).obsScal_;
            for (size_t j = 0; j < jjSize; j++)
            {
                statsMat(i, j) = tmpMap[keys.at(j)];
            }
        }

        ClusterMatrix_t stddevs = arma::stddev(statsMat, 1, 0); // size = 1 x statsMat.n_cols
        ClusterMatrix_t means = arma::mean(statsMat, 0);

        //END STATS===================================================================

        jjResult.clear();
        double sqrtNWorkers = std::sqrt(double(nworkers));
        for (size_t i = 0; i < keys.size(); i++)
        {
            jjResult[keys.at(i)] = {means(0, i), stddevs(0, i) / sqrtNWorkers}; // the ith key contains a vector of size 2, mean and stddev
        }

        jjResult["NWorkers"] = {nworkers, 0.0};

        std::ofstream fout("Obs.json");
        fout << std::setw(4) << jjResult << std::endl;
        fout.close();
    }

    static void SaveTabular(const std::string &fname, const ClusterMatrixCD_t &greenTab, const double &beta,
                            const size_t &precision = 10, const bool &saveArma = false)
    {

        std::ofstream fout;
        fout.open(fname + std::string(".dat"), std::ios::out);

        for (size_t nn = 0; nn < greenTab.n_rows; nn++)
        {
            const double iwn = (2.0 * nn + 1.0) * M_PI / beta;
            fout << std::setprecision(precision) << iwn << " ";

            for (Site_t ii = 0; ii < greenTab.n_cols; ii++)
            {

                fout << std::setprecision(precision) << greenTab(nn, ii).real()
                     << " "
                     << std::setprecision(precision) << greenTab(nn, ii).imag()
                     << " ";
            }
            fout << "\n";
        }

        fout.close();

        if (saveArma)
        {
            assert(greenTab.save(fname + std::string(".arma"), arma::arma_ascii));
        }
    }

    static void SaveFillingMatrixs(std::valarray<double> &fillingResultUp, std::valarray<double> &fillingResultDown, const IO::Base_IOModel &ioModel)
    {
        ClusterMatrix_t nUpMatrix(ioModel.Nc, ioModel.Nc);
        nUpMatrix.zeros();
        ClusterMatrix_t nDownMatrix = nUpMatrix;

        for (size_t kk = 0; kk < ioModel.fillingSites().size(); kk++)
        {

            for (size_t ii = 0; ii < nUpMatrix.n_rows; ii++)
            {
                std::pair<size_t, size_t> fSites = ioModel.indepSites().at(ioModel.FindIndepSiteIndex(ii, ii));
                assert(fSites.first == fSites.second);
                auto fillingSites = ioModel.fillingSites();
                const size_t fsitesfirst = fSites.first;
                auto fIt = std::find(fillingSites.begin(), fillingSites.end(), fsitesfirst);

                if (fIt == fillingSites.end())
                {
                    throw std::runtime_error("Bad index in fIt nMatrix, SaveFillingMatrix");
                }

                size_t fsite = std::distance(fillingSites.begin(), fIt);
                nUpMatrix(ii, ii) = fillingResultUp[fsite];
                nDownMatrix(ii, ii) = fillingResultDown[fsite];
            }
        }
        assert(nUpMatrix.save("nUpMatrix.dat", arma::raw_ascii));
        assert(nDownMatrix.save("nDownMatrix.dat", arma::raw_ascii));
    }
};

} //namespace mpiUt