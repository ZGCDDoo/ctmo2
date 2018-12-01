#pragma once

#include "Integrator.hpp"
#include "Utilities.hpp"
#include "MPIUtilities.hpp"
#include "GreenMat.hpp"
#include "ABC_SelfConsistency.hpp"
#include "Fourier_DCA.hpp"

namespace SelfCon
{

using Utilities::GetSpinName;

template <typename TH0>
struct GreenLattice
{

  public:
    static const size_t Nc;
    static const ClusterMatrixCD_t II;

    GreenLattice(cd_t zz, ClusterMatrixCD_t selfEnergy, TH0 h0) : zz_(zz), selfEnergy_(selfEnergy), h0_(h0){};

    ClusterMatrixCD_t operator()(const double &kx, const double &ky)
    {
        return ((zz_ * II - h0_(kx, ky) - selfEnergy_).i());
    }

  private:
    const cd_t zz_;
    ClusterMatrixCD_t selfEnergy_;
    TH0 h0_;
};
template <typename TH0>
const ClusterMatrixCD_t GreenLattice<TH0>::II = ClusterMatrixCD_t(TH0::Nc, TH0::Nc).eye();

template <typename TH0>
const size_t GreenLattice<TH0>::Nc = TH0::Nc;

template <typename TIOModel, typename TModel, typename TH0>
class SelfConsistency : public ABC_SelfConsistency
{

  public:
    static const size_t Nc;
    static const ClusterMatrixCD_t II;
    static const double factNSelfCon;
    const size_t hybSavePrecision = 14;

    SelfConsistency(const Json &jj, const TModel &model, const ClusterCubeCD_t &greenImpurity, const FermionSpin_t &spin) : model_(model),
                                                                                                                            ioModel_(TIOModel()),
                                                                                                                            h0_(model_.h0()),
                                                                                                                            greenImpurity_(FourierDCA::RtoK(greenImpurity, h0_.RSites(), h0_.KWaveVectors())),
                                                                                                                            hybridization_(spin == FermionSpin_t::Up ? model_.hybridizationMatUp() : model_.hybridizationMatDown()),
                                                                                                                            selfEnergy_(),
                                                                                                                            hybNext_(),
                                                                                                                            spin_(spin),
                                                                                                                            weights_(cd_t(jj["WEIGHTSR"].get<double>(), jj["WEIGHTSI"].get<double>())),
                                                                                                                            NOrb_(jj["NOrb"].get<size_t>())
    {
        mpiUt::Print("Start of SC constructor");

        const size_t NGreen = greenImpurity_.n_slices;
        const size_t NSelfCon = NGreen;

        selfEnergy_.resize(Nc, Nc, NSelfCon);
        selfEnergy_.zeros();

        //0.) Extraire la self jusqu'a NGreen
        for (size_t nn = 0; nn < NGreen; nn++)
        {
            const cd_t zz = cd_t(model_.mu(), (2.0 * nn + 1.0) * M_PI / model_.beta());
            selfEnergy_.slice(nn) = -greenImpurity_.slice(nn).i() + zz * ClusterMatrixCD_t(Nc, Nc).eye() - model_.tLoc() - hybridization_.slice(nn);
        }

        if (mpiUt::Rank() == mpiUt::master)
        {
            ioModel_.SaveK("self" + GetSpinName(spin_), selfEnergy_, model_.beta(), NOrb_, hybSavePrecision);
            std::cout << "In Selfonsistency constructor, after save selfenery " << std::endl;
        }

        mpiUt::Print("After SC constructor");
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
        if (mpiUt::Rank() == mpiUt::master)
        {
            std::cout << "In Selfonsistency DOSC serial" << std::endl;
            const size_t NSelfCon = selfEnergy_.n_slices;
            const size_t NKPTS = h0_.NKPTS();
            ClusterCubeCD_t gImpUpNext(Nc, Nc, NSelfCon);
            assert(Nc == h0_.KWaveVectors().size());
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
                    for (size_t kxindex = 0; kxindex < NKPTS; kxindex++)
                    {
                        const double kx = (Kx - kxCenter) + static_cast<double>(kxindex) / static_cast<double>(NKPTS - 1) * 2.0 * kxCenter;
                        for (size_t kyindex = 0; kyindex < NKPTS; kyindex++)
                        {
                            const double ky = (Ky - kyCenter) + static_cast<double>(kyindex) / static_cast<double>(NKPTS - 1) * 2.0 * kyCenter;
                            for (size_t kzindex = 0; kzindex < NKPTS; kzindex++)
                            {
                                const double kz = (Kz - kzCenter) + static_cast<double>(kzindex) / static_cast<double>(NKPTS - 1) * 2.0 * kzCenter;
                                gImpUpNext(KIndex, KIndex, nn) += 1.0 / (zz - h0_.Eps0k(kx, ky, kz) - selfEnergy_(KIndex, KIndex, nn));
                            }
                        }
                    }
                    gImpUpNext(KIndex, KIndex, nn) /= static_cast<double>(NKPTS * NKPTS);
                }
            }

            ioModel_.SaveK("green" + GetSpinName(spin_), gImpUpNext, model_.beta(), NOrb_, hybSavePrecision);

            for (size_t nn = 0; nn < gImpUpNext.n_slices; nn++)
            {
                const cd_t zz = cd_t(model_.mu(), (2.0 * nn + 1.0) * M_PI / model_.beta());
                hybNext_.slice(nn) = -gImpUpNext.slice(nn).i() - selfEnergy_.slice(nn) + zz * ClusterMatrixCD_t(Nc, Nc).eye() - model_.tLoc();
            }

            ioModel_.SaveK("hybNext" + GetSpinName(spin_), hybNext_, model_.beta(), NOrb_, hybSavePrecision);

            std::cout << "After Selfonsistency DOSC serial" << std::endl;
        }
    }

#ifdef HAVEMPI
    void DoSCGridParallel()
    {

        mpi::communicator world;

        mpiUt::Print("In Selfonsistency DOSC Parallel");
        const size_t NSelfCon = selfEnergy_.n_slices;

        if (static_cast<size_t>(mpiUt::NWorkers()) > NSelfCon)
        {
            DoSCGridSerial();
            return;
        }

        const size_t NSelfConRank = mpiUt::Rank() == mpiUt::master ? (NSelfCon / mpiUt::NWorkers() + NSelfCon % mpiUt::NWorkers()) : NSelfCon / mpiUt::NWorkers();

        ClusterCubeCD_t gImpUpNextRank(Nc, Nc, NSelfConRank);
        gImpUpNextRank.zeros();
        ClusterCubeCD_t hybNextRank(Nc, Nc, NSelfConRank);
        hybNextRank.zeros();

        const double kxCenter = M_PI / static_cast<double>(h0_.Nx);
        const double kyCenter = M_PI / static_cast<double>(h0_.Ny);
        const double kzCenter = M_PI / static_cast<double>(h0_.Nz);

        const size_t kxtildepts = (std::abs(h0_.txVec().at(0)) < 1e-10) ? 1 : NKPTS;
        const size_t kytildepts = (std::abs(h0_.tyVec().at(0)) < 1e-10) ? 1 : NKPTS;
        const size_t kztildepts = (std::abs(h0_.tzVec().at(0)) < 1e-10) ? 1 : NKPTS;

        const size_t nnStart = mpiUt::Rank() == mpiUt::master ? 0 : NSelfCon % mpiUt::NWorkers() + (NSelfCon / mpiUt::NWorkers()) * mpiUt::Rank();
        const size_t nnEnd = nnStart + NSelfConRank;
        for (size_t KIndex = 0; KIndex < h0_.KWaveVectors().size(); KIndex++)
        {

            const double Kx = h0_.KWaveVectors().at(KIndex)(0);
            const double Ky = h0_.KWaveVectors().at(KIndex)(1);
            const double Kz = h0_.KWaveVectors().at(KIndex)(2);

            for (size_t nn = nnStart; nn < nnEnd; nn++)
            {
                const cd_t zz = cd_t(model_.mu(), (2.0 * nn + 1.0) * M_PI / model_.beta());
                for (size_t kxindex = 0; kxindex < NKPTS; kxindex++)
                {
                    const double kx = (Kx - kxCenter) + static_cast<double>(kxindex) / static_cast<double>(NKPTS - 1) * 2.0 * kxCenter;
                    for (size_t kyindex = 0; kyindex < NKPTS; kyindex++)
                    {
                        const double ky = (Ky - kyCenter) + static_cast<double>(kyindex) / static_cast<double>(NKPTS - 1) * 2.0 * kyCenter;
                        gImpUpNextRank(KIndex, KIndex, nn - nnStart) += 1.0 / (zz - h0_.Eps0k(kx, ky) - selfEnergy_(KIndex, KIndex, nn));
                    }
                }
                gImpUpNextRank(KIndex, KIndex, nn - nnStart) /= static_cast<double>(NKPTS * NKPTS);
                hybNextRank(KIndex, KIndex, nn - nnStart) = -1.0 / gImpUpNextRank(KIndex, KIndex, nn - nnStart) - selfEnergy_(KIndex, KIndex, nn) + zz - model_.tLoc()(KIndex, KIndex);
            }
        }

        std::vector<std::vector<cd_t>> tmpMemGImpVec;
        std::vector<std::vector<cd_t>> tmpMemHybNextVec;
        std::vector<cd_t> tmpMemGImp = mpiUt::CubeCDToVecCD(gImpUpNextRank);
        std::vector<cd_t> tmpMemHybNext = mpiUt::CubeCDToVecCD(hybNextRank);

        if (mpiUt::Rank() == mpiUt::master)
        {
            mpi::gather(world, tmpMemGImp, tmpMemGImpVec, mpiUt::master);
            mpi::gather(world, tmpMemHybNext, tmpMemHybNextVec, mpiUt::master);
        }
        else
        {
            mpi::gather(world, tmpMemGImp, mpiUt::master);
            mpi::gather(world, tmpMemHybNext, mpiUt::master);
        }

        if (mpiUt::Rank() == mpiUt::master)
        {
            ClusterCubeCD_t gImpUpNext(Nc, Nc, NSelfCon);
            gImpUpNext.zeros();
            hybNext_.resize(Nc, Nc, NSelfCon);
            hybNext_.zeros();

            for (size_t ii = 0; ii < static_cast<size_t>(mpiUt::NWorkers()); ii++)
            {
                ClusterCubeCD_t tmpGImpNextRank = mpiUt::VecCDToCubeCD(tmpMemGImpVec.at(ii), Nc, Nc, tmpMemGImpVec.at(ii).size() / (Nc * Nc));
                ClusterCubeCD_t tmpHybNextRank = mpiUt::VecCDToCubeCD(tmpMemHybNextVec.at(ii), Nc, Nc, tmpMemHybNextVec.at(ii).size() / (Nc * Nc));

                const size_t jjStart = ii == 0 ? 0 : NSelfCon % mpiUt::NWorkers() + (NSelfCon / mpiUt::NWorkers()) * ii;
                const size_t jjEnd = jjStart + tmpGImpNextRank.n_slices;
                for (size_t jj = jjStart; jj < jjEnd; jj++)
                {
                    gImpUpNext.slice(jj) = tmpGImpNextRank.slice(jj - jjStart);
                    hybNext_.slice(jj) = tmpHybNextRank.slice(jj - jjStart);
                }
            }

            hybNext_ *= (1.0 - weights_);
            hybNext_ += weights_ * hybridization_.data();
            ioModel_.SaveK("green" + GetSpinName(spin_), gImpUpNext, model_.beta(), NOrb_, hybSavePrecision);
            ioModel_.SaveK("hybNext" + GetSpinName(spin_), hybNext_, model_.beta(), NOrb_, hybSavePrecision);

            mpiUt::Print("After Selfonsistency DOSC Parallel");
        }
    }

#endif

    ClusterCubeCD_t hybNext() const
    {
        return hybNext_;
    };

  private:
    TModel model_;
    TIOModel ioModel_;
    TH0 h0_;

    const ClusterCubeCD_t greenImpurity_;
    GreenMat::HybridizationMat hybridization_;
    ClusterCubeCD_t selfEnergy_;
    ClusterCubeCD_t hybNext_;
    const FermionSpin_t spin_;
    const cd_t weights_;
    const size_t NOrb_;
};
template <typename TIOModel, typename TModel, typename TH0>
const ClusterMatrixCD_t SelfConsistency<TIOModel, TModel, TH0>::II = ClusterMatrixCD_t(TH0::Nc, TH0::Nc).eye();

template <typename TIOModel, typename TModel, typename TH0>
const size_t SelfConsistency<TIOModel, TModel, TH0>::Nc = TH0::Nc;

template <typename TIOModel, typename TModel, typename TH0>
const double SelfConsistency<TIOModel, TModel, TH0>::factNSelfCon = 2;

} // namespace SelfCon
