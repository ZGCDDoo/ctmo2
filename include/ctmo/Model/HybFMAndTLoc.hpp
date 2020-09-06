#pragma once

#include "ctmo/Model/ABC_H0.hpp"
#include "ctmo/Foundations/Conventions.hpp"
#include "ctmo/Foundations/Fourier_DCA.hpp"

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
        const size_t NKPTS = h0.NKPTS();
        const double kxCenter = M_PI / static_cast<double>(h0.Nx);
        const double kyCenter = M_PI / static_cast<double>(h0.Ny);
        const double kzCenter = M_PI / static_cast<double>(h0.Nz);

        const size_t kxtildepts = (std::abs(h0.txVec().at(0)) < 1e-10) ? 1 : NKPTS;
        const size_t kytildepts = (std::abs(h0.tyVec().at(0)) < 1e-10) ? 1 : NKPTS;
        const size_t kztildepts = (std::abs(h0.tzVec().at(0)) < 1e-10) ? 1 : NKPTS;

        ClusterMatrixCD_t tlocK(h0.n_rows(), h0.n_rows());
        tlocK.zeros();
        ClusterMatrixCD_t hybFMK(h0.n_rows(), h0.n_rows());
        hybFMK.zeros();

        for (size_t KIndex = 0; KIndex < h0.KWaveVectors().size(); KIndex++)
        {
            const double Kx = h0.KWaveVectors().at(KIndex)(0);
            const double Ky = h0.KWaveVectors().at(KIndex)(1);
            const double Kz = h0.KWaveVectors().at(KIndex)(2);

            for (size_t kxindex = 0; kxindex < kxtildepts; kxindex++)
            {
                const double kx = (Kx - kxCenter) + static_cast<double>(kxindex) / static_cast<double>(NKPTS - 1) * 2.0 * kxCenter;
                for (size_t kyindex = 0; kyindex < kytildepts; kyindex++)
                {
                    const double ky = (Ky - kyCenter) + static_cast<double>(kyindex) / static_cast<double>(NKPTS - 1) * 2.0 * kyCenter;
                    for (size_t kzindex = 0; kzindex < kztildepts; kzindex++)
                    {
                        const double kz = (Kz - kzCenter) + static_cast<double>(kzindex) / static_cast<double>(NKPTS - 1) * 2.0 * kzCenter;
                        const double epsk = h0.Eps0k(kx, ky, kz);
                        tlocK(KIndex, KIndex) += epsk;
                        hybFMK(KIndex, KIndex) += epsk * epsk;
                    }
                }
            }
            tlocK(KIndex, KIndex) /= static_cast<double>(kxtildepts * kytildepts * kztildepts);
            hybFMK(KIndex, KIndex) /= static_cast<double>(kxtildepts * kytildepts * kztildepts);
        }
        hybFMK -= tlocK * tlocK;

        tlocK.save(tlocFName, arma::arma_ascii);
        hybFMK.save(hybFMFName, arma::arma_ascii);
        Logging::Debug("End of CalculateHybFMAndTLoc");
    }
}; // namespace Models

} // namespace Models
