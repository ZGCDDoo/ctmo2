#pragma once

#include "../Models/ABC_H0.hpp"
#include "../Utilities/Integrator.hpp"
#include "../Utilities/Conventions.hpp"
#include "../Utilities/Fourier_DCA.hpp"
#include "../Utilities/Logging.hpp"

namespace Models
{

class HybFMAndTLoc
{

  public:
    static void CalculateHybFMAndTLoc(const Models::ABC_H0 &h0)
    {
        Logging::Debug("Start of CalculateHybFMAndTLoc: This should only be for DCA ! ");

        Conventions::MapSS_t mapNames = Conventions::BuildFileNameConventions();
        const std::string tlocFName = mapNames.at("tlocFile"); // tloc File Name
        const std::string hybFMFName = mapNames.at("hybFMFile");

        using boost::filesystem::exists;
        if ((exists(tlocFName)) && exists(hybFMFName))
        {
            ClusterMatrixCD_t tmp;
            tmp.load(tlocFName);
            if (tmp.n_cols == h0.Nc)
            {
                return;
            }
        }

        // Get TLoc = Int[t(ktilde)]
        Logging::Warn("Calculating tLoc");

        TKTildeK tktildeK(h0);
        const ClusterMatrixCD_t tlocK = Integrator::CubatureKTildeDCA(tktildeK);

        // Get Int[ t(ktilde)^2 ]
        Logging::Warn("Calculating hybFM");
        TKTildeSquaredK tktildesquaredK(h0);
        const ClusterMatrixCD_t tktildeSquaredKIntegrated = Integrator::CubatureKTildeDCA(tktildesquaredK);

        const ClusterMatrixCD_t hybFMK = tktildeSquaredKIntegrated - tlocK * tlocK;

        tlocK.save(tlocFName, arma::arma_binary);
        hybFMK.save(hybFMFName, arma::arma_binary);
        Logging::Debug("End of CalculateHybFMAndTLoc");
    }

    struct TKTildeK
    {

        explicit TKTildeK(const Models::ABC_H0 &h0) : h0_(h0), Nx(h0.Nx), Ny(h0.Ny), Nz(h0.Nz), Nc(h0.Nc){};
        size_t n_rows() const { return h0_.n_rows(); };
        size_t n_cols() const { return h0_.n_cols(); };

        Models::ABC_H0 h0_;
        const size_t Nx;
        const size_t Ny;
        const size_t Nz;
        const size_t Nc;

        ClusterMatrixCD_t operator()(const double &kTildeX, const double &kTildeY, const double &kTildeZ) // return t(ktilde)^2
        {
            return (FourierDCA::RtoK(h0_(kTildeX, kTildeY, kTildeZ), h0_.RSites(), h0_.KWaveVectors()));
        }
    };

    struct TKTildeSquaredK
    {

        explicit TKTildeSquaredK(const Models::ABC_H0 &h0) : h0_(h0), Nx(h0.Nx), Ny(h0.Ny), Nz(h0.Nz), Nc(h0.Nc){};
        size_t n_rows() const { return h0_.n_rows(); };
        size_t n_cols() const { return h0_.n_cols(); };

        Models::ABC_H0 h0_;
        const size_t Nx;
        const size_t Ny;
        const size_t Nz;
        const size_t Nc;

        ClusterMatrixCD_t operator()(const double &kTildeX, const double &kTildeY, const double &kTildeZ) // return t(ktilde)^2
        {
            const ClusterMatrixCD_t tmp = h0_(kTildeX, kTildeY, kTildeZ);
            return (FourierDCA::RtoK(tmp * tmp, h0_.RSites(), h0_.KWaveVectors()));
        }
    };
};

} // namespace Models
