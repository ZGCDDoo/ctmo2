#pragma once
#include "../Utilities/Utilities.hpp"

namespace Models
{
template <size_t TNX, size_t TNY>
class ABC_H0
{

  public:
    static const size_t Nc;
    static const size_t Nx;
    static const size_t Ny;
    static const size_t n_rows;
    static const size_t n_cols;

    ABC_H0(const ABC_H0 &abc_h0) = default;
    ABC_H0(const Json &tJson) : RSites_(Nc),
                                KWaveVectors_(Nc)

    {
        assert(TNX == TNY);

        for (size_t i = 0; i < TNX; i++)
        {
            for (size_t j = 0; j < TNY; j++)
            {
                const size_t index = i + TNY * j;
                RSites_.at(index) = {static_cast<double>(i), static_cast<double>(j)};
                KWaveVectors_.at(index) = {static_cast<double>(i) * 2.0 * M_PI / static_cast<double>(TNX), static_cast<double>(j) * 2.0 * M_PI / static_cast<double>(TNY)};
            }
        }

        ReadInHoppings(tJson);
    }

    virtual ~ABC_H0() = 0;

    // double t() const { return t_; };
    // double tPrime() const { return tPrime_; };
    // double tPrimePrime() const { return tPrimePrime_; };
    ClusterSites_t RSites() const { return RSites_; };
    ClusterSites_t KWaveVectors() const { return KWaveVectors_; };

    void ReadInHoppings(const Json &tJson)
    {
        const size_t NOrb = tJson["NOrb"].get<size_t>();
        assert(tJson["tParameters"].size() == NOrb * (NOrb + 1) / 2);

        for (size_t o1 = 0; o1 < NOrb; o1++)
        {
            for (size_t o2 = o1; o2 < NOrb; o2++)
            {

                const std::string o1o2Name = std::to_string(o1) + std::to_string(o2);
                const Json &jj = jj["tParameters"][o1o2Name];

                tIntraOrbitalVec_.push_back(jj[o1o2Name]["tIntra"].get<double>());
                txVec_.push_back(jj[o1o2Name]["tx"].get<double>());
                tyVec_.push_back(jj[o1o2Name]["ty"].get<double>());
                txyVec_.push_back(jj[o1o2Name]["tx=y"].get<double>());
                tx_yVec_.push_back(jj[o1o2Name]["tx=-y"].get<double>());
                t2xVec_.push_back(jj[o1o2Name]["t2x"].get<double>());
                t2yVec_.push_back(jj[o1o2Name]["t2y"].get<double>());
            }
        }
    }

    double Eps0k(const double &kx, const double &ky, const size_t &NIndepOrbIndex)
    {
        const double eps0k = tIntraOrbitalVec_.at(NIndepOrbIndex) +
                             2.0 * txVec_.at(NIndepOrbIndex) * std::cos(kx) +
                             2.0 * tyVec_.at(NIndepOrbIndex) * std::cos(ky) +
                             2.0 * txyVec_.at(NIndepOrbIndex) * std::cos(kx + ky) +
                             2.0 * tx_yVec_.at(NIndepOrbIndex) * std::cos(kx - ky) +
                             2.0 * t2xVec_.at(NIndepOrbIndex) * std::cos(2.0 * kx) +
                             2.0 * t2yVec_.at(NIndepOrbIndex) * std::cos(2.0 * ky);

        return eps0k;
    }

    // ClusterMatrixCD_t
    // operator()(const double &kTildeX, const double &kTildeY) //return t(ktilde)
    // {
    //     const cd_t im = cd_t(0.0, 1.0);
    //     const SiteVector_t ktilde = {kTildeX, kTildeY};
    //     ClusterMatrixCD_t HoppingKTilde(Nc, Nc);
    //     HoppingKTilde.zeros();

    //     for (size_t i = 0; i < Nc; i++)
    //     {
    //         for (size_t j = 0; j < Nc; j++)
    //         {
    //             for (const SiteVector_t &K : this->KWaveVectors_)
    //             {
    //                 HoppingKTilde(i, j) += std::exp(im * dot(K + ktilde, RSites_.at(i) - RSites_[j])) * Eps0k(K(0) + kTildeX, K(1) + kTildeY);
    //             }
    //         }
    //     }
    //     return (HoppingKTilde / static_cast<double>(Nc));
    // }

    // void SaveTKTildeAndHybFM()
    // {
    //     //check if  file exists:
    //     using boost::filesystem::exists;
    //     if ((exists("tktilde.arma") && exists("tloc.arma")) && exists("hybFM.arma"))
    //     {
    //         ClusterMatrixCD_t tmp;
    //         tmp.load("tloc.arma");
    //         if (tmp.n_cols == Nc)
    //         {
    //             return;
    //         }
    //     }
    //     std::cout << "Calculating tktilde, tloc and hybFM. " << std::endl;

    //     const size_t kxtildepts = 2.0 * M_PI / 0.009 / std::min(Nx, Ny);

    //     ClusterCubeCD_t tKTildeGrid(Nc, Nc, kxtildepts * kxtildepts);
    //     tKTildeGrid.zeros();
    //     ClusterMatrixCD_t tLoc(Nc, Nc);
    //     tLoc.zeros();

    //     size_t sliceindex = 0;
    //     for (size_t kx = 0; kx < kxtildepts; kx++)
    //     {
    //         const double kTildeX = static_cast<double>(kx) / static_cast<double>(kxtildepts) * 2.0 * M_PI / static_cast<double>(Nx);
    //         for (size_t ky = 0; ky < kxtildepts; ky++)
    //         {
    //             const double kTildeY = static_cast<double>(ky) / static_cast<double>(kxtildepts) * 2.0 * M_PI / static_cast<double>(Ny);
    //             tKTildeGrid.slice(sliceindex) = (*this)(kTildeX, kTildeY);
    //             tLoc += tKTildeGrid.slice(sliceindex);
    //             sliceindex++;
    //         }
    //     }

    //     tKTildeGrid.save("tktilde.arma");
    //     tLoc /= static_cast<double>(tKTildeGrid.n_slices);
    //     tLoc.save("tloc.arma", arma::arma_ascii);

    //     //First moment of hyb
    //     ClusterMatrixCD_t hybFM(Nc, Nc);
    //     hybFM.zeros();

    //     const size_t Nkpts = tKTildeGrid.n_slices;
    //     for (size_t nn = 0; nn < Nkpts; nn++)
    //     {
    //         hybFM += tKTildeGrid.slice(nn) * tKTildeGrid.slice(nn);
    //     }
    //     hybFM /= Nkpts;
    //     hybFM -= tLoc * tLoc;
    //     hybFM.save("hybFM.arma", arma::arma_ascii);
    // }

  protected:
    ClusterSites_t RSites_;
    ClusterSites_t KWaveVectors_;

    const std::vector<double> tIntraOrbitalVec_;
    const std::vector<double> txVec_;
    const std::vector<double> tyVec_;
    const std::vector<double> txyVec_;
    const std::vector<double> tx_yVec_; //along x=-y diagonal
    const std::vector<double> t2xVec_;
    const std::vector<double> t2yVec_;
};

template <size_t TNX, size_t TNY>
ABC_H0<TNX, TNY>::~ABC_H0() {} //destructors must exist

template <size_t TNX, size_t TNY>
const size_t ABC_H0<TNX, TNY>::Nx = TNX;

template <size_t TNX, size_t TNY>
const size_t ABC_H0<TNX, TNY>::Ny = TNY;

template <size_t TNX, size_t TNY>
const size_t ABC_H0<TNX, TNY>::Nc = TNX *TNY;

template <size_t TNX, size_t TNY>
const size_t ABC_H0<TNX, TNY>::n_rows = TNX *TNY;

template <size_t TNX, size_t TNY>
const size_t ABC_H0<TNX, TNY>::n_cols = TNX *TNY;

} // namespace Models
