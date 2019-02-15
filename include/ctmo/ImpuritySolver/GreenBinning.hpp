#pragma once

#include "ctmo/Foundations/LinAlg.hpp"
#include "ctmo/ImpuritySolver/ISData.hpp"
#include "ctmo/Foundations/UtilitiesRandom.hpp"

namespace Markov
{
namespace Obs
{
const size_t N_BIN_TAU = 100000;
using IOModel_t = IO::Base_IOModel;
using Model_t = Models::ABC_Model_2D;

class GreenBinning
{

  public:
    GreenBinning(const std::shared_ptr<ISDataCT> &dataCT, const Json &jjSim, const FermionSpin_t &spin)
        : dataCT_(dataCT), modelPtr_(dataCT_->modelPtr_), ioModelPtr_(modelPtr_->ioModelPtr()),
          NMat_(0.5 * (jjSim["solver"]["eCutGreen"].get<double>() * dataCT_->beta() / M_PI - 1.0)), spin_(spin),
          NOrb_(jjSim["model"]["nOrb"].get<size_t>())
    {

        const size_t LL = ioModelPtr_->GetNIndepSuperSites(NOrb_);
        M0Bins_.resize(LL);
        M1Bins_.resize(LL);
        M2Bins_.resize(LL);
        M3Bins_.resize(LL);

        for (size_t ii = 0; ii < M0Bins_.size(); ++ii)
        {
            M0Bins_.at(ii).resize(N_BIN_TAU, 0.0);
            M1Bins_.at(ii).resize(N_BIN_TAU, 0.0);
            M2Bins_.at(ii).resize(N_BIN_TAU, 0.0);
            M3Bins_.at(ii).resize(N_BIN_TAU, 0.0);
        }

        Logging::Trace("GreenBinning Created.");
    }

    ClusterCubeCD_t greenCube() const { return greenCube_; };

    void MeasureGreenBinning(const Matrix<double> &Mmat)
    {

        const size_t kkSpin = (spin_ == FermionSpin_t::Up) ? dataCT_->vertices_.NUp() : dataCT_->vertices_.NDown();
        const double DeltaInv = N_BIN_TAU / dataCT_->beta_;
        if (kkSpin)
        {
            for (size_t p1 = 0; p1 < kkSpin; ++p1)
            {
                for (size_t p2 = 0; p2 < kkSpin; ++p2)
                {
                    const SuperSite_t s1 =
                        (spin_ == FermionSpin_t::Up) ? dataCT_->vertices_.atUp(p1).superSite() : dataCT_->vertices_.atDown(p1).superSite();
                    const SuperSite_t s2 =
                        (spin_ == FermionSpin_t::Up) ? dataCT_->vertices_.atUp(p2).superSite() : dataCT_->vertices_.atDown(p2).superSite();
                    const size_t ll = ioModelPtr_->FindIndepSuperSiteIndex(s1, s2, NOrb_);
                    double temp = static_cast<double>(dataCT_->sign_) * Mmat(p1, p2);

                    double tau = (spin_ == FermionSpin_t::Up) ? dataCT_->vertices_.atUp(p1).tau() - dataCT_->vertices_.atUp(p2).tau()
                                                              : dataCT_->vertices_.atDown(p1).tau() - dataCT_->vertices_.atDown(p2).tau();
                    if (tau < 0.0)
                    {
                        temp *= -1.0;
                        tau += dataCT_->beta_;
                    }

                    const int index = DeltaInv * tau;
                    const double dTau = tau - (static_cast<double>(index) + 0.5) / DeltaInv;

                    M0Bins_.at(ll).at(index) += temp;
                    temp *= dTau;
                    M1Bins_[ll][index] += temp;
                    temp *= dTau;
                    M2Bins_[ll][index] += temp;
                    temp *= dTau;
                    M3Bins_[ll][index] += temp;
                }
            }
        }
    }

    ClusterCubeCD_t FinalizeGreenBinning(const double &signMeas, const size_t &NMeas)
    {
        Logging::Debug("Start of GreenBinning.FinalizeGreenBinning()");

        const double dTau = dataCT_->beta_ / N_BIN_TAU;
        SiteVectorCD_t indep_M_matsubara_sampled(ioModelPtr_->GetNIndepSuperSites(NOrb_));
        const ClusterCubeCD_t green0CubeMatsubara =
            spin_ == FermionSpin_t::Up ? modelPtr_->greenCluster0MatUp().data() : modelPtr_->greenCluster0MatDown().data();
        ClusterCubeCD_t greenCube(NOrb_ * ioModelPtr_->Nc, NOrb_ * ioModelPtr_->Nc, NMat_);
        greenCube.zeros();

        for (size_t n = 0; n < NMat_; ++n)
        {
            const double omega_n = M_PI * (2.0 * n + 1.0) / dataCT_->beta_;
            const cd_t iomega_n(0.0, omega_n);
            const cd_t fact = std::exp(iomega_n * dTau);
            const double lambda =
                2.0 * std::sin(omega_n * dTau / 2.0) / (dTau * omega_n * (1.0 - omega_n * omega_n * dTau * dTau / 24.0) * NMeas);

            for (size_t ll = 0; ll < ioModelPtr_->GetNIndepSuperSites(NOrb_); ++ll)
            {
                cd_t temp_matsubara = 0.0;

                const size_t llSite = ll % ioModelPtr_->indepSites().size(); // ll / ioModelPtr_->indepSites().size();
                cd_t exp_factor = std::exp(iomega_n * dTau / 2.0) /
                                  (static_cast<double>(ioModelPtr_->nOfAssociatedSites().at(llSite))); // watch out important factor!
                for (size_t ii = 0; ii < N_BIN_TAU; ii++)
                {
                    cd_t coeff = lambda * exp_factor;

                    temp_matsubara += coeff * M0Bins_.at(ll).at(ii);
                    temp_matsubara += coeff * M1Bins_[ll][ii] * iomega_n;
                    temp_matsubara += coeff * M2Bins_[ll][ii] * iomega_n * iomega_n / 2.0;
                    temp_matsubara += coeff * M3Bins_[ll][ii] * iomega_n * iomega_n * iomega_n / 6.0;

                    exp_factor *= fact;
                }
                indep_M_matsubara_sampled(ll) = temp_matsubara;
            }
            const ClusterMatrixCD_t dummy1 = ioModelPtr_->IndepToFull(indep_M_matsubara_sampled, NOrb_);
            const ClusterMatrixCD_t green0 = green0CubeMatsubara.slice(n);

            greenCube.slice(n) = green0 - green0 * dummy1 * green0 / (dataCT_->beta_ * signMeas);
        }

        greenCube_ = greenCube; // in case it is needed later on

        Logging::Debug("End of GreenBinning.FinalizeGreenBinning()");
        return greenCube; // the  measured interacting green function
    }

  private:
    std::shared_ptr<ISDataCT> dataCT_;
    std::shared_ptr<Model_t> modelPtr_;
    std::shared_ptr<IOModel_t> ioModelPtr_;

    std::vector<std::vector<double>> M0Bins_;
    std::vector<std::vector<double>> M1Bins_;
    std::vector<std::vector<double>> M2Bins_;
    std::vector<std::vector<double>> M3Bins_;

    ClusterCubeCD_t greenCube_;

    const size_t NMat_;
    const FermionSpin_t spin_;
    const size_t NOrb_;
};

} // namespace Obs
} // namespace Markov