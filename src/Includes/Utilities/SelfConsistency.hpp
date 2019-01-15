#pragma once

#include "Integrator.hpp"
#include "MPITools.hpp"
#include "ABC_SelfConsistency.hpp"
#include "../Models/ABC_Model.hpp"

namespace SelfCon
{

using Utilities::GetSpinName;

class SelfConsistency : public ABC_SelfConsistency
{

    using Model_t = Models::ABC_Model_2D;
    using IOModel_t = IO::Base_IOModel;

  public:
    SelfConsistency(const Json &jjSim, const Model_t &model, const ClusterCubeCD_t &greenImpurity, const FermionSpin_t &spin) : model_(model),
                                                                                                                                ioModel_(jjSim),
                                                                                                                                greenImpurity_(greenImpurity),
                                                                                                                                hybridization_(spin == FermionSpin_t::Up ? model_.hybridizationMatUp() : model_.hybridizationMatDown()),
                                                                                                                                selfEnergy_(),
                                                                                                                                hybNext_(),
                                                                                                                                spin_(spin),
                                                                                                                                weights_(jjSim["selfCon"]["weightsR"].get<double>(), jjSim["selfCon"]["weightsI"].get<double>()),
                                                                                                                                NOrb_(model.NOrb()),
                                                                                                                                NSS_(NOrb_ * ioModel_.Nc)

    {

        Logging::Debug("Start of SC constructor.");

        const size_t NGreen = greenImpurity_.n_slices;
        size_t NSelfConTmp = std::max<double>(0.5 * (jjSim["selfCon"]["eCutSelfCon"].get<double>() * model_.beta() / M_PI - 1.0),
                                              0.5 * (200.0 * model_.beta() / M_PI - 1.0));
        if (NGreen >= NSelfConTmp)
        {
            NSelfConTmp = factNSelfCon_ * static_cast<double>(NGreen);
        }
        const size_t NSelfCon = NGreen; //NSelfConTmp;
        assert(NSelfCon >= NGreen);
        //Patcher la hyb si necessaire
        hybridization_.PatchHF(NSelfCon, model_.beta());
        const size_t NHyb = hybridization_.n_slices();
        assert(NHyb >= NSelfCon);

        selfEnergy_.resize(NSS_, NSS_, NSelfCon);

        //0.) Extraire la self jusqu'a NGreen
        for (size_t nn = 0; nn < NGreen; ++nn)
        {
            const cd_t zz(model_.mu(), (2.0 * nn + 1.0) * M_PI / model_.beta());
            selfEnergy_.slice(nn) = -greenImpurity_.slice(nn).i() + zz * ClusterMatrixCD_t(NSS_, NSS_).eye() - model_.tLoc() - hybridization_.slice(nn);
        }

        if (mpiUt::Tools::Rank() == mpiUt::Tools::master)
        {
            ioModel_.SaveCube("self" + GetSpinName(spin_), selfEnergy_, model_.beta(), NOrb_, hybSavePrecision_);
            Logging::Info("In Selfonsistency constructor, after save selfenery. ");
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

#ifdef HAVEMPI
    void DoSCGridParallel()
    {

        mpi::communicator world;

        Logging::Info("In Selfonsistency DOSC Parallel.");
        const size_t NSelfCon = selfEnergy_.n_slices;

        if (static_cast<size_t>(mpiUt::Tools::NWorkers()) > NSelfCon)
        {
            DoSCGridSerial();
            return;
        }

        const size_t NSelfConRank = mpiUt::Tools::Rank() == mpiUt::Tools::master ? (NSelfCon / mpiUt::Tools::NWorkers() + NSelfCon % mpiUt::Tools::NWorkers()) : NSelfCon / mpiUt::Tools::NWorkers();

        ClusterCubeCD_t gImpUpNextRank(NSS_, NSS_, NSelfConRank);
        gImpUpNextRank.zeros();
        ClusterCubeCD_t hybNextRank(NSS_, NSS_, NSelfConRank);
        hybNextRank.zeros();

        ClusterCubeCD_t tKTildeGrid;
        assert(tKTildeGrid.load("tktilde.arma"));
        const size_t ktildepts = tKTildeGrid.n_slices;

        const size_t nnStart = mpiUt::Tools::Rank() == mpiUt::Tools::master ? 0 : NSelfCon % mpiUt::Tools::NWorkers() + (NSelfCon / mpiUt::Tools::NWorkers()) * mpiUt::Tools::Rank();
        const size_t nnEnd = nnStart + NSelfConRank;
        for (size_t nn = nnStart; nn < nnEnd; ++nn)
        {
            const cd_t zz = cd_t(model_.mu(), (2.0 * nn + 1.0) * M_PI / model_.beta());
            for (size_t ktildeindex = 0; ktildeindex < ktildepts; ++ktildeindex)
            {
                gImpUpNextRank.slice(nn - nnStart) += (zz * ClusterMatrixCD_t(NSS_, NSS_).eye() - tKTildeGrid.slice(ktildeindex) - selfEnergy_.slice(nn)).i();
            }
            gImpUpNextRank.slice(nn - nnStart) /= static_cast<double>(ktildepts);
            hybNextRank.slice(nn - nnStart) = -gImpUpNextRank.slice(nn - nnStart).i() - selfEnergy_.slice(nn) + zz * ClusterMatrixCD_t(NSS_, NSS_).eye() - model_.tLoc();
        }

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
            ClusterCubeCD_t gImpUpNext(NSS_, NSS_, NSelfCon);
            gImpUpNext.zeros();
            hybNext_.resize(NSS_, NSS_, NSelfCon);
            hybNext_.zeros();

            for (size_t ii = 0; ii < static_cast<size_t>(mpiUt::Tools::NWorkers()); ++ii)
            {
                ClusterCubeCD_t tmpGImpNextRank = mpiUt::Tools::VecCDToCubeCD(tmpMemGImpVec.at(ii), NSS_, NSS_, tmpMemGImpVec.at(ii).size() / (NSS_ * NSS_));
                ClusterCubeCD_t tmpHybNextRank = mpiUt::Tools::VecCDToCubeCD(tmpMemHybNextVec.at(ii), NSS_, NSS_, tmpMemHybNextVec.at(ii).size() / (NSS_ * NSS_));

                const size_t jjStart = ii == 0 ? 0 : NSelfCon % mpiUt::Tools::NWorkers() + (NSelfCon / mpiUt::Tools::NWorkers()) * ii;
                const size_t jjEnd = jjStart + tmpGImpNextRank.n_slices;
                for (size_t jj = jjStart; jj < jjEnd; jj++)
                {
                    gImpUpNext.slice(jj) = tmpGImpNextRank.slice(jj - jjStart);
                    hybNext_.slice(jj) = tmpHybNextRank.slice(jj - jjStart);
                }
            }

            hybNext_ *= (1.0 - weights_);
            hybNext_ += weights_ * hybridization_.data();
            ioModel_.SaveCube("green" + GetSpinName(spin_), gImpUpNext, model_.beta(), NOrb_, hybSavePrecision_);
            ioModel_.SaveCube("hybNext" + GetSpinName(spin_), hybNext_, model_.beta(), NOrb_, hybSavePrecision_);

            Logging::Info("After Selfonsistency DOSC Parallel");
        }
    }

#endif

    void DoSCGridSerial()
    {

        if (mpiUt::Tools::Rank() == mpiUt::Tools::master)
        {
            Logging::Info("In Selfonsistency DOSC serial.");
            const size_t NSelfCon = selfEnergy_.n_slices;
            ClusterCubeCD_t gImpUpNext(NSS_, NSS_, NSelfCon);
            gImpUpNext.zeros();
            hybNext_.resize(NSS_, NSS_, NSelfCon);
            hybNext_.zeros();
            ClusterCubeCD_t tKTildeGrid;
            assert(tKTildeGrid.load("tktilde.arma"));
            size_t ktildepts = tKTildeGrid.n_slices;

            for (size_t nn = 0; nn < NSelfCon; ++nn)
            {
                const cd_t zz = cd_t(model_.mu(), (2.0 * static_cast<double>(nn) + 1.0) * M_PI / model_.beta());
                for (size_t ktildeindex = 0; ktildeindex < ktildepts; ++ktildeindex)
                {
                    gImpUpNext.slice(nn) += (zz * ClusterMatrixCD_t(NSS_, NSS_).eye() - tKTildeGrid.slice(ktildeindex) - selfEnergy_.slice(nn)).i();
                }
                gImpUpNext.slice(nn) /= static_cast<double>(ktildepts);
                hybNext_.slice(nn) = -gImpUpNext.slice(nn).i() - selfEnergy_.slice(nn) + zz * ClusterMatrixCD_t(NSS_, NSS_).eye() - model_.tLoc();
            }

            hybNext_ *= (1.0 - weights_);
            hybNext_ += weights_ * hybridization_.data();
            ioModel_.SaveCube("green" + GetSpinName(spin_), gImpUpNext, model_.beta(), NOrb_, hybSavePrecision_);
            ioModel_.SaveCube("hybNext" + GetSpinName(spin_), hybNext_, model_.beta(), NOrb_, hybSavePrecision_);

            Logging::Info("After Selfonsistency DOSC serial.");
        }
    }

    ClusterCubeCD_t
    hybNext() const
    {
        return hybNext_;
    };

  private:
    Model_t model_;
    IOModel_t ioModel_;

    const ClusterCubeCD_t greenImpurity_;
    GreenMat::HybridizationMat hybridization_;
    ClusterCubeCD_t selfEnergy_;
    ClusterCubeCD_t hybNext_;
    const FermionSpin_t spin_;
    const cd_t weights_;
    const size_t NOrb_;
    const size_t NSS_; //Number of super-sites : (orbital and sites)

    const double factNSelfCon_ = 2;
    const size_t hybSavePrecision_ = 14;
};

} // namespace SelfCon