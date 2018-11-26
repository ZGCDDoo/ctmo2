#pragma once

#include "../Utilities/Utilities.hpp"
#include "../Utilities/Integrator.hpp"
#include "../Utilities/Conventions.hpp"
#include "../Utilities/Fourier_DCA.hpp"

namespace Models
{

template <typename TH0>
class HybFMAndTLoc
{

  public:
    static void CalculateHybFMAndTLoc(const TH0 &h0)
    {
        std::cout << "Start of CalculateHybFMAndTLoc: This should only be for DCA !!!! " << std::endl;

        Conventions::MapSS_t mapNames = Conventions::BuildFileNameConventions();
        const std::string tlocFName = mapNames["tlocFile"]; //tloc File Name
        const std::string hybFMFName = mapNames["hybFMFile"];

        using boost::filesystem::exists;
        if ((exists(tlocFName)) && exists(hybFMFName))
        {
            ClusterMatrixCD_t tmp;
            tmp.load(tlocFName);
            if (tmp.n_cols == TH0::Nc)
            {
                return;
            }
        }

        //Get TLoc = Int[t(ktilde)]
        std::cout << "Calculating tLoc" << std::endl;
        TKTildeK tktildeK(h0);
        const ClusterMatrixCD_t tlocK = Integrator::CubatureKTildeDCA(tktildeK);

        //Get Int[ t(ktilde)^2 ]
        std::cout << "Calculating hybFM" << std::endl;
        TKTildeSquaredK tktildesquaredK(h0);
        const ClusterMatrixCD_t tktildeSquaredKIntegrated = Integrator::CubatureKTildeDCA(tktildesquaredK);

        const ClusterMatrixCD_t hybFMK = tktildeSquaredKIntegrated - tlocK * tlocK;

        tlocK.save(tlocFName, arma::arma_binary);
        hybFMK.save(hybFMFName, arma::arma_binary);
        std::cout << "End of CalculateHybFMAndTLoc" << std::endl;
    }

    struct TKTildeK
    {

        TKTildeK(const TH0 &h0) : h0_(h0){};
        size_t n_rows() const { return h0_.n_rows(); };
        size_t n_cols() const { return h0_.n_cols(); };

        const size_t Nx = TH0::Nx;
        const size_t Ny = TH0::Ny;
        const size_t Nz = TH0::Nz;
        const size_t Nc = TH0::Nc;

        TH0 h0_;

        ClusterMatrixCD_t operator()(const double &kTildeX, const double &kTildeY, const double &kTildeZ) //return t(ktilde)^2
        {
            return (FourierDCA::RtoK(h0_(kTildeX, kTildeY, kTildeZ), h0_.RSites(), h0_.KWaveVectors()));
        }
    };

    struct TKTildeSquaredK
    {

        TKTildeSquaredK(const TH0 &h0) : h0_(h0){};
        size_t n_rows() const { return h0_.n_rows(); };
        size_t n_cols() const { return h0_.n_cols(); };

        const size_t Nx = TH0::Nx;
        const size_t Ny = TH0::Ny;
        const size_t Nz = TH0::Nz;
        const size_t Nc = TH0::Nc;

        TH0 h0_;

        ClusterMatrixCD_t operator()(const double &kTildeX, const double &kTildeY, const double &kTildeZ) //return t(ktilde)^2
        {
            const ClusterMatrixCD_t tmp = h0_(kTildeX, kTildeY, kTildeZ);
            return (FourierDCA::RtoK(tmp * tmp, h0_.RSites(), h0_.KWaveVectors()));
        }
    };
};

} // namespace Models
