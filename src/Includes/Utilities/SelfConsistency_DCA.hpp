#pragma once

#include "Integrator.hpp"
#include "MPITools.hpp"
#include "ABC_SelfConsistency.hpp"
#include "../Models/ABC_Model.hpp"
#include "Fourier_DCA.hpp"

namespace SelfCon
{

using Utilities::GetSpinName;

class SelfConsistencyDCA : public ABC_SelfConsistency
{

    using Model_t = Models::ABC_Model_2D;
    using IOModel_t = IO::Base_IOModel;

  public:
    const size_t hybSavePrecision = 14;

    SelfConsistencyDCA(const Json &jjSim, const Model_t &model, const ClusterCubeCD_t &greenImpurity, const FermionSpin_t &spin) : model_(model),
                                                                                                                                   ioModel_(jjSim),
                                                                                                                                   h0_(model_.h0()),
                                                                                                                                   greenImpurity_(FourierDCA::RtoK(greenImpurity, h0_.RSites(), h0_.KWaveVectors())),
                                                                                                                                   hybridization_(spin == FermionSpin_t::Up ? model_.hybridizationMatUp() : model_.hybridizationMatDown()),
                                                                                                                                   selfEnergy_(),
                                                                                                                                   hybNext_(),
                                                                                                                                   spin_(spin),
                                                                                                                                   weights_(jjSim["selfCon"]["weightsR"].get<double>(), jjSim["selfCon"]["weightsI"].get<double>()),
                                                                                                                                   NOrb_(model.NOrb()),
                                                                                                                                   Nc_(ioModel_.Nc)
    {
        Logging::Debug("Start of SC constructor.");

        const size_t NGreen = greenImpurity_.n_slices;
        const size_t NSelfCon = NGreen;

        selfEnergy_.resize(Nc_, Nc_, NSelfCon);
        selfEnergy_.zeros();
        hybridization_.PatchHF(NSelfCon, model_.beta());

        //0.) Extraire la self jusqu'a NGreen
        for (size_t nn = 0; nn < NGreen; nn++)
        {
            const cd_t zz = cd_t(model_.mu(), (2.0 * nn + 1.0) * M_PI / model_.beta());
            selfEnergy_.slice(nn) = -greenImpurity_.slice(nn).i() + zz * ClusterMatrixCD_t(Nc_, Nc_).eye() - model_.tLoc() - hybridization_.slice(nn);
        }

        if (mpiUt::Tools::Rank() == mpiUt::Tools::master)
        {
            ioModel_.SaveK("self" + GetSpinName(spin_), selfEnergy_, model_.beta(), NOrb_, hybSavePrecision);
            Logging::Debug("In Selfonsistency constructor, after save selfenery. ");
        }

        Logging::Debug("After SC constructor.");
    }

    void DoSCGrid() override
    {
#ifdef HAVEMPI
        DoSCGridParallel();
#else
        DoSCGridSerial();
#endif
    }

    void DoSCGridSerial()
    {
        if (mpiUt::Tools::Rank() == mpiUt::Tools::master)
        {
            Logging::Info("In Selfonsistency DOSC serial.");
            const size_t NSelfCon = selfEnergy_.n_slices;
            const size_t NKPTS = h0_.NKPTS();
            ClusterCubeCD_t gImpUpNext(Nc_, Nc_, NSelfCon);
            assert(Nc_ == h0_.KWaveVectors().size());
            gImpUpNext.zeros();
            hybNext_ = gImpUpNext;

            const double kxCenter = M_PI / static_cast<double>(h0_.Nx);
            const double kyCenter = M_PI / static_cast<double>(h0_.Ny);
            const double kzCenter = M_PI / static_cast<double>(h0_.Nz);

            const size_t kxtildepts = (std::abs(h0_.txVec().at(0)) < 1e-10) ? 1 : NKPTS;
            const size_t kytildepts = (std::abs(h0_.tyVec().at(0)) < 1e-10) ? 1 : NKPTS;
            const size_t kztildepts = (std::abs(h0_.tzVec().at(0)) < 1e-10) ? 1 : NKPTS;

            for (size_t KIndex = 0; KIndex < h0_.KWaveVectors().size(); KIndex++)
            {
                const double Kx = h0_.KWaveVectors().at(KIndex)(0);
                const double Ky = h0_.KWaveVectors().at(KIndex)(1);
                const double Kz = h0_.KWaveVectors().at(KIndex)(2);

                for (size_t nn = 0; nn < NSelfCon; nn++)
                {
                    const cd_t zz = cd_t(model_.mu(), (2.0 * nn + 1.0) * M_PI / model_.beta());
                    for (size_t kxindex = 0; kxindex < kxtildepts; kxindex++)
                    {
                        const double kx = (Kx - kxCenter) + static_cast<double>(kxindex) / static_cast<double>(NKPTS - 1) * 2.0 * kxCenter;
                        for (size_t kyindex = 0; kyindex < kytildepts; kyindex++)
                        {
                            const double ky = (Ky - kyCenter) + static_cast<double>(kyindex) / static_cast<double>(NKPTS - 1) * 2.0 * kyCenter;
                            for (size_t kzindex = 0; kzindex < kztildepts; kzindex++)
                            {
                                const double kz = (Kz - kzCenter) + static_cast<double>(kzindex) / static_cast<double>(NKPTS - 1) * 2.0 * kzCenter;
                                gImpUpNext(KIndex, KIndex, nn) += 1.0 / (zz - h0_.Eps0k(kx, ky, kz, 0) - selfEnergy_(KIndex, KIndex, nn));
                            }
                        }
                    }
                    gImpUpNext(KIndex, KIndex, nn) /= static_cast<double>(kxtildepts * kytildepts * kztildepts);
                }
            }

            ioModel_.SaveK("green" + GetSpinName(spin_), gImpUpNext, model_.beta(), NOrb_, hybSavePrecision);

            for (size_t nn = 0; nn < gImpUpNext.n_slices; nn++)
            {
                const cd_t zz = cd_t(model_.mu(), (2.0 * nn + 1.0) * M_PI / model_.beta());
                hybNext_.slice(nn) = -gImpUpNext.slice(nn).i() - selfEnergy_.slice(nn) + zz * ClusterMatrixCD_t(Nc_, Nc_).eye() - model_.tLoc();
            }

            ioModel_.SaveK("hybNext" + GetSpinName(spin_), hybNext_, model_.beta(), NOrb_, hybSavePrecision);

            Logging::Info("After Selfonsistency DOSC serial.");
        }
    }

#ifdef HAVEMPI
    void DoSCGridParallel()
    {

        mpi::communicator world;

        Logging::Info("In Selfonsistency DOSC Parallel");
        const size_t NSelfCon = selfEnergy_.n_slices;

        if (static_cast<size_t>(mpiUt::Tools::NWorkers()) > NSelfCon)
        {
            DoSCGridSerial();
            return;
        }

        const size_t NSelfConRank = mpiUt::Tools::Rank() == mpiUt::Tools::master ? (NSelfCon / mpiUt::Tools::NWorkers() + NSelfCon % mpiUt::Tools::NWorkers()) : NSelfCon / mpiUt::Tools::NWorkers();

        ClusterCubeCD_t gImpUpNextRank(Nc_, Nc_, NSelfConRank);
        gImpUpNextRank.zeros();
        ClusterCubeCD_t hybNextRank(Nc_, Nc_, NSelfConRank);
        hybNextRank.zeros();

        const size_t NKPTS = h0_.NKPTS();

        const double kxCenter = M_PI / static_cast<double>(h0_.Nx);
        const double kyCenter = M_PI / static_cast<double>(h0_.Ny);
        const double kzCenter = M_PI / static_cast<double>(h0_.Nz);

        const size_t kxtildepts = (std::abs(h0_.txVec().at(0)) < 1e-10) ? 1 : NKPTS;
        const size_t kytildepts = (std::abs(h0_.tyVec().at(0)) < 1e-10) ? 1 : NKPTS;
        const size_t kztildepts = (std::abs(h0_.tzVec().at(0)) < 1e-10) ? 1 : NKPTS;

        const size_t nnStart = mpiUt::Tools::Rank() == mpiUt::Tools::master ? 0 : NSelfCon % mpiUt::Tools::NWorkers() + (NSelfCon / mpiUt::Tools::NWorkers()) * mpiUt::Tools::Rank();
        const size_t nnEnd = nnStart + NSelfConRank;

        Logging::Trace("Just before loops");

        for (size_t KIndex = 0; KIndex < h0_.KWaveVectors().size(); KIndex++)
        {

            const double Kx = h0_.KWaveVectors().at(KIndex)(0);
            const double Ky = h0_.KWaveVectors().at(KIndex)(1);
            const double Kz = h0_.KWaveVectors().at(KIndex)(2);

            for (size_t nn = nnStart; nn < nnEnd; nn++)
            {
                const cd_t zz = cd_t(model_.mu(), (2.0 * nn + 1.0) * M_PI / model_.beta());
                for (size_t kxindex = 0; kxindex < kxtildepts; kxindex++)
                {
                    const double kx = (Kx - kxCenter) + static_cast<double>(kxindex) / static_cast<double>(NKPTS - 1) * 2.0 * kxCenter;
                    for (size_t kyindex = 0; kyindex < kytildepts; kyindex++)
                    {
                        const double ky = (Ky - kyCenter) + static_cast<double>(kyindex) / static_cast<double>(NKPTS - 1) * 2.0 * kyCenter;

                        for (size_t kzindex = 0; kzindex < kztildepts; kzindex++)
                        {
                            const double kz = (Kz - kzCenter) + static_cast<double>(kzindex) / static_cast<double>(NKPTS - 1) * 2.0 * kzCenter;
                            gImpUpNextRank(KIndex, KIndex, nn - nnStart) += 1.0 / (zz - h0_.Eps0k(kx, ky, kz, 0) - selfEnergy_(KIndex, KIndex, nn));
                        }
                    }
                }
                gImpUpNextRank(KIndex, KIndex, nn - nnStart) /= static_cast<double>(kxtildepts * kytildepts * kztildepts);
                hybNextRank(KIndex, KIndex, nn - nnStart) = -1.0 / gImpUpNextRank(KIndex, KIndex, nn - nnStart) - selfEnergy_(KIndex, KIndex, nn) + zz - model_.tLoc()(KIndex, KIndex);
            }
        }

        Logging::Trace("After loops");

        std::vector<std::vector<cd_t>> tmpMemGImpVec;
        std::vector<std::vector<cd_t>> tmpMemHybNextVec;
        std::vector<cd_t> tmpMemGImp = mpiUt::Tools::CubeCDToVecCD(gImpUpNextRank);
        std::vector<cd_t> tmpMemHybNext = mpiUt::Tools::CubeCDToVecCD(hybNextRank);

        if (mpiUt::Tools::Rank() == mpiUt::Tools::master)
        {
            mpi::gather(world, tmpMemGImp, tmpMemGImpVec, mpiUt::Tools::master);
            mpi::gather(world, tmpMemHybNext, tmpMemHybNextVec, mpiUt::Tools::master);
        }
        else
        {
            mpi::gather(world, tmpMemGImp, mpiUt::Tools::master);
            mpi::gather(world, tmpMemHybNext, mpiUt::Tools::master);
        }

        if (mpiUt::Tools::Rank() == mpiUt::Tools::master)
        {
            ClusterCubeCD_t gImpUpNext(Nc_, Nc_, NSelfCon);
            gImpUpNext.zeros();
            hybNext_.resize(Nc_, Nc_, NSelfCon);
            hybNext_.zeros();

            Logging::Trace("Before merge.");

            for (size_t ii = 0; ii < static_cast<size_t>(mpiUt::Tools::NWorkers()); ii++)
            {
                ClusterCubeCD_t tmpGImpNextRank = mpiUt::Tools::VecCDToCubeCD(tmpMemGImpVec.at(ii), Nc_, Nc_, tmpMemGImpVec.at(ii).size() / (Nc_ * Nc_));
                ClusterCubeCD_t tmpHybNextRank = mpiUt::Tools::VecCDToCubeCD(tmpMemHybNextVec.at(ii), Nc_, Nc_, tmpMemHybNextVec.at(ii).size() / (Nc_ * Nc_));

                Logging::Trace("In Loop.");

                const size_t jjStart = ii == 0 ? 0 : NSelfCon % mpiUt::Tools::NWorkers() + (NSelfCon / mpiUt::Tools::NWorkers()) * ii;
                const size_t jjEnd = jjStart + tmpGImpNextRank.n_slices;
                for (size_t jj = jjStart; jj < jjEnd; jj++)
                {
                    gImpUpNext.slice(jj) = tmpGImpNextRank.slice(jj - jjStart);
                    hybNext_.slice(jj) = tmpHybNextRank.slice(jj - jjStart);
                }
            }

            Logging::Trace("After merge.");

            hybNext_ *= (1.0 - weights_);
            hybNext_ += weights_ * hybridization_.data();
            ioModel_.SaveK("green" + GetSpinName(spin_), gImpUpNext, model_.beta(), NOrb_, hybSavePrecision);
            ioModel_.SaveK("hybNext" + GetSpinName(spin_), hybNext_, model_.beta(), NOrb_, hybSavePrecision);

            Logging::Info("After Selfonsistency DOSC Parallel");
        }
    }

#endif

    ClusterCubeCD_t hybNext() const
    {
        return hybNext_;
    };

  private:
    Model_t model_;
    IOModel_t ioModel_;
    Models::ABC_H0 h0_;

    const ClusterCubeCD_t greenImpurity_;
    GreenMat::HybridizationMat hybridization_;
    ClusterCubeCD_t selfEnergy_;
    ClusterCubeCD_t hybNext_;
    const FermionSpin_t spin_;
    const cd_t weights_;
    const size_t NOrb_;
    const size_t Nc_;
};

} // namespace SelfCon
