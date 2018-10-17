#pragma once

#define GREEN_STYLE

#include <valarray>

#include "../Utilities/Utilities.hpp"
#include "../Utilities/LinAlg.hpp"
#include "../Utilities/Matrix.hpp"
#include "../Utilities/MPIUtilities.hpp"
#include "../Utilities/Fourier.hpp"
#include "../Utilities/GreenTau.hpp"
#include "Obs/Observables.hpp"
#include "ISData.hpp"

//#define DEBUG_TEST

namespace Markov
{

using Fourier::MatToTau;
using Fourier::MatToTauCluster;
using Vertex = Diagrammatic::Vertex;
using VertexPart = Diagrammatic::VertexPart;
using NFData = Diagrammatic::NFData;

typedef LinAlg::Matrix_t Matrix_t;

struct UpdData
{

    UpdData() : NQUp_(), NQDown_(), newLastRowUp_(), newLastRowDown_(){};
    SiteVector_t NQUp_;
    SiteVector_t NQDown_;
    SiteVector_t newLastRowUp_;
    SiteVector_t newLastRowDown_;
    double sTildeUpI_;
    double sTildeDownI_;
};

template <typename TIOModel, typename TModel>
class ABC_MarkovChain
{

    using GreenTau_t = GreenTau::GreenCluster0Tau<TIOModel>;

  public:
    const size_t Nc = TModel::Nc;
    const double PROBFLIP = 0.25;
    const double PROBINSERT = 0.333;
    const double PROBREMOVE = 1.0 - PROBINSERT;

    ABC_MarkovChain(const Json &jj, const size_t &seed) : modelPtr_(new TModel(jj)),
                                                          rng_(seed),
                                                          urng_(rng_, Utilities::UniformDistribution_t(0.0, 1.0)),
                                                          nfdata_(),
                                                          upddata_(),
                                                          dataCT_(
                                                              new Obs::ISDataCT<TIOModel, TModel>(
                                                                  jj,
                                                                  *modelPtr_)),
                                                          obs_(dataCT_, jj),
                                                          vertexBuilder_(jj, TModel::Nc),
                                                          updsamespin_(0)
    {
        const std::valarray<size_t> zeroPair = {0, 0};
        updStats_["Inserts"] = zeroPair;
        updStats_["Removes"] = zeroPair;
        updStats_["Flips"] = zeroPair;
        updatesProposed_ = 0;

        mpiUt::Print("MarkovChain Created \n");
    }

    virtual ~ABC_MarkovChain() = 0;

    //Getters
    TModel model() const
    {
        return (*modelPtr_);
    };

    Matrix_t Nup() const
    {
        return nfdata_.Nup_;
    };

    Matrix_t Ndown() const
    {
        return nfdata_.Ndown_;
    };

    size_t updatesProposed() const { return updatesProposed_; }

    double beta() const
    {
        return dataCT_->beta_;
    };

    virtual double gammaTrad(const FermionSpin_t &spin, const AuxSpin_t &auxTo, const AuxSpin_t &vauxFrom) = 0;
    virtual double FAux(const FermionSpin_t &spin, const AuxSpin_t &aux) = 0;

    // void ThermalizeFromConfig()
    // {
    //     if (mpiUt::LoadConfig(dataCT_->vertices_))
    //     {
    //         const size_t kk = dataCT_->vertices_.size();
    //         nfdata_.FVup_ = SiteVector_t(kk);
    //         nfdata_.FVdown_ = SiteVector_t(kk);
    //         for (size_t i = 0; i < kk; i++)
    //         {
    //             AuxSpin_t aux = dataCT_->vertices_.at(i).aux();
    //             nfdata_.FVup_(i) = FAux(FermionSpin_t::Up, aux);
    //             nfdata_.FVdown_(i) = FAux(FermionSpin_t::Down, aux);
    //         }

    //         nfdata_.Nup_.Resize(kk, kk);
    //         nfdata_.Ndown_.Resize(kk, kk);
    //         CleanUpdate();
    //     }
    // }

    void DoStep()
    {

        urng_() < PROBINSERT ? InsertVertex() : RemoveVertex();
        updatesProposed_++;
    }

    void AssertSizes()
    {
        const size_t kk = dataCT_->vertices_.size();
        // std::cout << "kk, Nup_.n_rows, Ndown_.n_rows = " << kk << ", " << nfdata_.Nup_.n_rows() << ", " << nfdata_.Ndown_.n_rows() << std::endl;
        assert(nfdata_.Nup_.n_rows() + nfdata_.Ndown_.n_rows() == 2 * kk);
        assert(2 * dataCT_->vertices_.size() == dataCT_->vertices_.NUp() + dataCT_->vertices_.NDown());
        // assert(dataCT_->vertices_.NUp() == dataCT_->vertices_.NDown());
        // std::cout << "dataCT_->vertices_.NUp(), nfdata_.FVup_.n_elem = " << dataCT_->vertices_.NUp() << ", " << nfdata_.FVup_.n_elem << std::endl;
        // std::cout << "dataCT_->vertices_.NDown(), nfdata_.FVdown_.n_elem = " << dataCT_->vertices_.NDown() << ", " << nfdata_.FVdown_.n_elem << std::endl;

        assert(dataCT_->vertices_.NUp() == nfdata_.Nup_.n_rows());

        assert(dataCT_->vertices_.NDown() == nfdata_.Ndown_.n_rows());

        // assert(!nfdata_.Nup_.HasInfOrNan());

        // assert(!nfdata_.Ndown_.HasInfOrNan());
    }

    void InsertVertex()
    {
        // std::cout << "Start InsertVertex " << std::endl;
        AssertSizes();
        updStats_["Inserts"][0]++;
        const Vertex vertex = vertexBuilder_.BuildVertex(urng_);
        const VertexPart x = vertex.vStart();
        const VertexPart y = vertex.vEnd();

        if (x.spin() == y.spin())
        {
            x.spin() == (FermionSpin_t::Up) ? InsertVertexSameSpin(vertex, nfdata_.Nup_, nfdata_.FVup_) : InsertVertexSameSpin(vertex, nfdata_.Ndown_, nfdata_.FVdown_);
        }
        else
        {
            InsertVertexDiffSpin(vertex);
        }

        // std::cout << "End InsertVertex " << std::endl;
    }

    void CalculateUpdDataInsertDiffSpin(const VertexPart &x)
    {
        const double faux = FAux(x.spin(), x.aux());

        if (x.spin() == FermionSpin_t::Up)
        {
            upddata_.sTildeUpI_ = -faux + GetGreenTau0(x, x);

            if (nfdata_.Nup_.n_rows())
            {
                const size_t kkoldUp = dataCT_->vertices_.NUp();
                upddata_.newLastRowUp_.set_size(kkoldUp);
                SiteVector_t newLastColUp_(kkoldUp);
                upddata_.NQUp_.set_size(kkoldUp);

                for (size_t iUp = 0; iUp < dataCT_->vertices_.NUp(); iUp++)
                {
                    upddata_.newLastRowUp_(iUp) = GetGreenTau0(x, dataCT_->vertices_.atUp(iUp));
                    newLastColUp_(iUp) = GetGreenTau0(dataCT_->vertices_.atUp(iUp), x);
                }
                MatrixVectorMult(nfdata_.Nup_, newLastColUp_, 1.0, upddata_.NQUp_);
                upddata_.sTildeUpI_ -= LinAlg::DotVectors(upddata_.newLastRowUp_, upddata_.NQUp_);
            }
        }
        else
        {
            upddata_.sTildeDownI_ = -faux + GetGreenTau0(x, x);

            if (nfdata_.Ndown_.n_rows())
            {
                const size_t kkoldDown = dataCT_->vertices_.NDown();
                upddata_.newLastRowDown_.set_size(kkoldDown);
                SiteVector_t newLastColDown_(kkoldDown);
                upddata_.NQDown_.set_size(kkoldDown);

                for (size_t iDown = 0; iDown < dataCT_->vertices_.NDown(); iDown++)
                {
                    upddata_.newLastRowDown_(iDown) = GetGreenTau0(x, dataCT_->vertices_.atDown(iDown));
                    newLastColDown_(iDown) = GetGreenTau0(dataCT_->vertices_.atDown(iDown), x);
                }

                MatrixVectorMult(nfdata_.Ndown_, newLastColDown_, 1.0, upddata_.NQDown_);
                upddata_.sTildeDownI_ -= LinAlg::DotVectors(upddata_.newLastRowDown_, upddata_.NQDown_);
            }
        }
    }

    void InsertVertexDiffSpin(const Vertex &vertex)
    {
        // std::cout << "Start InsertVertexDiffSpin " << std::endl;
        // const auto x = vertex.vStart();
        // const auto y = vertex.vEnd();
        AssertSizes();

        const VertexPart vertexPartUp = vertex.vStart();
        const VertexPart vertexPartDown = vertex.vEnd();
        const size_t kkoldUp = nfdata_.Nup_.n_rows();
        const size_t kknewUp = kkoldUp + 1;
        const size_t kkoldDown = nfdata_.Ndown_.n_rows();
        const size_t kknewDown = kkoldDown + 1;
        const double fauxup = FAux(vertexPartUp.spin(), vertexPartUp.aux());
        const double fauxdown = FAux(vertexPartDown.spin(), vertexPartDown.aux());

        assert(vertexPartUp.spin() == FermionSpin_t::Up);
        assert(vertexPartDown.spin() == FermionSpin_t::Down);
        assert(std::abs(vertexPartUp.tau() - vertexPartDown.tau()) < 1e-10);

        CalculateUpdDataInsertDiffSpin(vertexPartUp);
        CalculateUpdDataInsertDiffSpin(vertexPartDown);

        const double ratio = upddata_.sTildeUpI_ * upddata_.sTildeDownI_;
        const double ratioAcc = PROBREMOVE / PROBINSERT * vertex.probProb() / static_cast<size_t>(dataCT_->vertices_.size() + 1) * ratio;

        if (urng_() < std::abs(ratioAcc))
        {
            AssertSizes();

            updStats_["Inserts"][1]++;
            if (ratioAcc < .0)
            {
                dataCT_->sign_ *= -1;
            }

            if (nfdata_.Nup_.n_rows())
            {
                LinAlg::BlockRankOneUpgrade(nfdata_.Nup_, upddata_.NQUp_, upddata_.newLastRowUp_, 1.0 / upddata_.sTildeUpI_);
            }
            else
            {
                nfdata_.Nup_ = Matrix_t(1, 1);
                nfdata_.Nup_(0, 0) = 1.0 / upddata_.sTildeUpI_;
            }

            if (nfdata_.Ndown_.n_rows())
            {
                LinAlg::BlockRankOneUpgrade(nfdata_.Ndown_, upddata_.NQDown_, upddata_.newLastRowDown_, 1.0 / upddata_.sTildeDownI_);
            }
            else
            {
                nfdata_.Ndown_ = Matrix_t(1, 1);
                nfdata_.Ndown_(0, 0) = 1.0 / upddata_.sTildeDownI_;
            }

            dataCT_->vertices_.AppendVertex(vertex);
            AssertSizes();

            // if (nfdata_.Nup_.n_rows() > 10)
            // {
            //     dataCT_->vertices_.SwapSpin(2, 7, FermionSpin_t::Up);
            //     nfdata_.Nup_.SwapRowsAndCols(2, 7);
            //     nfdata_.FVup_.swap_rows(2, 7);
            // }
        }

        // std::cout << "End InsertVertexDiffSpin " << std::endl;
    }

    void InsertVertexSameSpin(const Vertex &vertex, Matrix_t &Nspin, SiteVector_t &FVspin)
    {
        // std::cout << "Start InsertVertexSameSpin " << std::endl;
        // return;
        // std::cout << "\n\n";
        // Nspin.Print();
        // std::cout << "\n\n";
        AssertSizes();

        const VertexPart x = vertex.vStart();
        const VertexPart y = vertex.vEnd();
        assert(x.spin() == y.spin());

        const double faux = FAux(x.spin(), vertex.aux());
        const double s00 = -faux + GetGreenTau0(x, x);
        const double s01 = GetGreenTau0(x, y);
        const double s10 = GetGreenTau0(y, x);
        const double s11 = -faux + GetGreenTau0(y, y);

        if (Nspin.n_rows())
        {
            AssertSizes();
            assert(Nspin.n_rows());
            const size_t kkold = dataCT_->vertices_.size();
            const size_t kknew = kkold + 1;
            const size_t kkoldspin = Nspin.n_rows();

            Matrix_t Q_(kkoldspin, 2);
            Matrix_t R_(2, kkoldspin);

            for (size_t i = 0; i < kkoldspin; i++)
            {

                const VertexPart vPartI = (x.spin() == FermionSpin_t::Up) ? dataCT_->vertices_.atUp(i) : dataCT_->vertices_.atDown(i);
                Q_(i, 0) = GetGreenTau0(vPartI, x);
                Q_(i, 1) = GetGreenTau0(vPartI, y);

                assert(vPartI.spin() == x.spin());

                R_(0, i) = GetGreenTau0(x, vPartI);
                R_(1, i) = GetGreenTau0(y, vPartI);
            }
            //     // std::cout << "In INsertvertex After loop " << std::endl;

            //     //Watch out, we are calculating two times the matrix NQ, once here and once in ranktwoupgrade. In a next version, only calculate here, not in ranktwoupgrade.
            //     // Matrix_t NQ(2 * kkold, 2); //NQ = N*Q
            //     // MatrixVectorMult(nfdata_.N_, Q_, 1.0, NQUp);

            //     // Matrix_t RNQ(2, 2); //R*NQ
            //     //     // Matrix_t RNQ(2, 2); //R*NQ

            Matrix_t sTilde = Matrix_t({{s00, s01}, {s10, s11}}) - LinAlg::DotRank2(R_, Nspin, Q_);
            sTilde.Inverse();

            // std::cout << "\n\n";
            // sTilde.Print();
            // std::cout << "\n\n";

            const double ratioAcc = PROBREMOVE / PROBINSERT * vertex.probProb() / kknew * 1.0 / sTilde.Determinant();
            AssertSizes();
            if (urng_() < std::abs(ratioAcc))
            {
                updStats_["Inserts"][1]++;
                if (ratioAcc < .0)
                {
                    dataCT_->sign_ *= -1;
                }

                LinAlg::BlockRankTwoUpgrade(Nspin, Q_, R_, sTilde);
                dataCT_->vertices_.AppendVertex(vertex);
                AssertSizes();
                updsamespin_++;
                // std::cout << "InsertVertexSameSpin accepted " << std::endl;
            }
        }
        else
        {
            // return;
            AssertSizes();
            Matrix_t sTilde = Matrix_t({{s00, s01}, {s10, s11}});
            sTilde.Inverse();
            const double ratioAcc = PROBREMOVE / PROBINSERT * vertex.probProb() * 1.0 / sTilde.Determinant();
            if (urng_() < std::abs(ratioAcc))
            {
                if (ratioAcc < 0.0)
                {
                    dataCT_->sign_ *= -1;
                }

                Nspin = sTilde;

                dataCT_->vertices_.AppendVertex(vertex);
            }
            AssertSizes();
        }
        // std::cout << "End InsertVertexSameSpin " << std::endl;

        // std::cout << "After insertvertex" << std::endl;
    }

    void RemoveVertex()
    {
        // std::cout << "Start RemoveVertex " << std::endl;

        AssertSizes();
        const size_t kk = dataCT_->vertices_.size();
        if (kk)
        {
            updStats_["Removes"][0]++;

            const size_t pp = static_cast<int>(urng_() * dataCT_->vertices_.size());
            const Vertex vertex = dataCT_->vertices_.at(pp);
            const VertexPart x = vertex.vStart();
            const VertexPart y = vertex.vEnd();
            assert(std::abs(x.tau() - y.tau()) < 1e-10);

            if (x.spin() == y.spin())
            {
                // PrepareToRemove(pp);
                x.spin() == (FermionSpin_t::Up) ? RemoveVertexSameSpin(pp, nfdata_.Nup_, nfdata_.FVup_) : RemoveVertexSameSpin(pp, nfdata_.Ndown_, nfdata_.FVdown_);
            }
            else
            {
                // PrepareToRemove(pp);
                RemoveVertexDiffSpin(pp);
            }
        }
        AssertSizes();
        // std::cout << "End RemoveVertex " << std::endl;
    }

    void RemoveVertexDiffSpin(const size_t &pp)
    {
        // std::cout << "Start RemoveVertexDiffSpin " << std::endl;
        // return;
        AssertSizes();

        const Vertex vertex = dataCT_->vertices_.at(pp);
        const size_t ppUp = dataCT_->vertices_.GetIndicesSpins(pp, FermionSpin_t::Up).at(0);
        const size_t ppDown = dataCT_->vertices_.GetIndicesSpins(pp, FermionSpin_t::Down).at(0);
        const auto x = dataCT_->vertices_.atUp(ppUp);
        const auto y = dataCT_->vertices_.atDown(ppDown);
        assert(std::abs(x.tau() - y.tau()) < 1e-10);
        // assert(x.superSite() == y.superSite());

        // std::cout << "here 1" << std::endl;
        // std::cout << "vertices.size(), NUp, NDown, pp , ppUp, ppDown = " << dataCT_->vertices_.size() << ", " << dataCT_->vertices_.NUp() << ", " << dataCT_->vertices_.NDown() << ", " << pp << ", " << ppUp << ", " << ppDown << std::endl;

        //In theory we should find the proper index for each spin
        const double ratioAcc = PROBINSERT / PROBREMOVE * static_cast<double>(dataCT_->vertices_.size()) / vertex.probProb() * nfdata_.Nup_(ppUp, ppUp) * nfdata_.Ndown_(ppDown, ppDown);
        // std::cout << "here 2" << std::endl;
        if (urng_() < std::abs(ratioAcc))
        {
            //AssertSizes();
            updStats_["Removes"][1]++;
            if (ratioAcc < .0)
            {
                dataCT_->sign_ *= -1;
            }

            //The update matrices of size k-1 x k-1 with the pp row and col deleted and the last row and col now at index pp

            const size_t kkUp = dataCT_->vertices_.NUp();
            const size_t kkUpm1 = kkUp - 1;
            const size_t kkDown = dataCT_->vertices_.NDown();
            const size_t kkDownm1 = kkDown - 1;
            dataCT_->vertices_.PrepareToRemove(nfdata_, pp); //vertices are already poped back

            LinAlg::BlockRankOneDowngrade(nfdata_.Nup_, kkUpm1);
            LinAlg::BlockRankOneDowngrade(nfdata_.Ndown_, kkDownm1);
            nfdata_.FVup_.resize(kkUpm1);
            nfdata_.FVdown_.resize(kkDownm1);

            AssertSizes();
        }
        // std::cout << "End RemoveVertexDiffSpin " << std::endl;
    }

    void RemoveVertexSameSpin(const size_t &pp, Matrix_t &Nspin, SiteVector_t &FVspin)
    {
        // std::cout << "Start RemoveVertexSameSpin " << std::endl;
        // return;
        AssertSizes();
        // std::cout << "Here 1 " << std::endl;
        assert(Nspin.n_rows() >= 2);
        // assert(FVspin.n_elem >= 2);

        const Vertex vertex = dataCT_->vertices_.at(pp);
        const VertexPart x = vertex.vStart();
        const VertexPart y = vertex.vEnd();
        assert(x.spin() == y.spin());

        assert(std::abs(x.tau() - y.tau()) < 1e-10);
        assert(x.orbital() != y.orbital());
        assert(x.site() == y.site());

        const auto indicesPP = dataCT_->vertices_.GetIndicesSpins(pp, x.spin());
        const size_t pp1Spin = indicesPP.at(0);
        const size_t pp2Spin = indicesPP.at(1);
        const size_t kk = dataCT_->vertices_.size();
        const size_t kkSpin = (x.spin() == FermionSpin_t::Up) ? dataCT_->vertices_.NUp() : dataCT_->vertices_.NDown();
        const size_t kkSpinm1 = kkSpin - 1;
        // const size_t kkSpinm2 = kkSpin - 2;

        const ClusterMatrix_t STildeInverse = {{Nspin(pp1Spin, pp1Spin), Nspin(pp1Spin, pp2Spin)}, {Nspin(pp2Spin, pp1Spin), Nspin(pp2Spin, pp2Spin)}};
        const double ratioAcc = PROBINSERT / PROBREMOVE * static_cast<double>(kk) / vertex.probProb() * arma::det(STildeInverse);
        // std::cout << "Here 2 " << std::endl;

        if (urng_() < std::abs(ratioAcc))
        {
            AssertSizes();
            updStats_["Removes"][1]++;
            if (ratioAcc < 0.0)
            {
                dataCT_->sign_ *= -1;
            }

            dataCT_->vertices_.PrepareToRemove(nfdata_, pp); //vertices are already poped back

            FVspin.resize(kkSpin - 2);
            LinAlg::BlockRankTwoDowngrade(Nspin);

            assert(Nspin.n_rows() == FVspin.n_elem);
        }

        AssertSizes();
        // std::cout << "End RemoveVertexSameSpin " << std::endl;
    }

    void CleanUpdate()
    {
        // mpiUt::Print("Cleaning, sign, k =  " + std::to_string(dataCT_->sign_) + ",  " + std::to_string(dataCT_->vertices_.size()));
        AssertSizes();
        // std::cout << "updsamespin = " << updsamespin_ << std::endl;
        // for (size_t ii = 0; ii < dataCT_->vertices_.size(); ii++)
        // {
        //     std::cout << dataCT_->vertices_.at(ii).vStart().orbital() << ", " << dataCT_->vertices_.at(ii).vEnd().orbital() << std::endl;
        // }

        const size_t kkup = dataCT_->vertices_.NUp();
        const size_t kkdown = dataCT_->vertices_.NDown();

        if (kkup != 0)
        {

            for (size_t iUp = 0; iUp < kkup; iUp++)
            {
                for (size_t jUp = 0; jUp < kkup; jUp++)
                {

                    nfdata_.Nup_(iUp, jUp) = GetGreenTau0(dataCT_->vertices_.atUp(iUp), dataCT_->vertices_.atUp(jUp));

                    if (iUp == jUp)
                    {
                        nfdata_.Nup_(iUp, iUp) -= FAux(FermionSpin_t::Up, dataCT_->vertices_.atUp(iUp).aux());
                    }
                }
            }
            nfdata_.Nup_.Inverse();
        }

        if (kkdown != 0)
        {
            for (size_t iDown = 0; iDown < kkdown; iDown++)
            {
                for (size_t jDown = 0; jDown < kkdown; jDown++)
                {

                    nfdata_.Ndown_(iDown, jDown) = GetGreenTau0(dataCT_->vertices_.atDown(iDown), dataCT_->vertices_.atDown(jDown));

                    if (iDown == jDown)
                    {
                        nfdata_.Ndown_(iDown, iDown) -= FAux(FermionSpin_t::Down, dataCT_->vertices_.atDown(iDown).aux());
                    }
                }
            }
            nfdata_.Ndown_.Inverse();
        }

        // std::cout << "End cleaning " << std::endl;
    }

    double GetGreenTau0(const VertexPart &x, const VertexPart &y) const
    {
        assert(x.spin() == y.spin());
        if (x.spin() == FermionSpin_t::Up)
        {
            return (dataCT_->green0CachedUp_(x.superSite(), y.superSite(), x.tau() - y.tau()));
        }
        else
        {
#ifdef AFM
            return (dataCT_->green0CachedDown_(x.superSite(), y.superSite(), x.tau() - y.tau()));
#else
            return (dataCT_->green0CachedUp_(x.superSite(), y.superSite(), x.tau() - y.tau()));
#endif
        }
    }

    void Measure()
    {
        AssertSizes();
        *(dataCT_->MupPtr_) = nfdata_.Nup_;
        *(dataCT_->MdownPtr_) = nfdata_.Ndown_;

        obs_.Measure();
    }

    void SaveMeas()
    {

        obs_.Save();
        std::cout << "updsamespin = " << updsamespin_ << std::endl;
        SaveUpd("upd.meas");
        if (mpiUt::Rank() == mpiUt::master)
        {
            dataCT_->vertices_.SaveConfig("Config.dat");
        }
    }

    void SaveTherm()
    {

        SaveUpd("upd.therm");
        for (UpdStats_t::iterator it = updStats_.begin(); it != updStats_.end(); ++it)
        {
            std::string key = it->first;
            updStats_[key] = 0.0;
        }
    }

    void SaveUpd(const std::string fname)
    {
        std::vector<UpdStats_t> updStatsVec;
#ifdef HAVEMPI

        mpi::communicator world;
        if (mpiUt::Rank() == mpiUt::master)
        {
            mpi::gather(world, updStats_, updStatsVec, mpiUt::master);
        }
        else
        {
            mpi::gather(world, updStats_, mpiUt::master);
        }
        if (mpiUt::Rank() == mpiUt::master)
        {
            mpiUt::SaveUpdStats(fname, updStatsVec);
        }

#else
        updStatsVec.push_back(updStats_);
        mpiUt::SaveUpdStats(fname, updStatsVec);
#endif

        mpiUt::Print("Finished Saving MarkovChain.");
    }

  protected:
    //attributes
    std::shared_ptr<TModel> modelPtr_;
    Utilities::EngineTypeMt19937_t rng_;
    Utilities::UniformRngMt19937_t urng_;
    NFData nfdata_;
    UpdData upddata_;
    std::shared_ptr<Obs::ISDataCT<TIOModel, TModel>> dataCT_;
    Obs::Observables<TIOModel, TModel> obs_;
    Diagrammatic::VertexBuilder vertexBuilder_;

    UpdStats_t updStats_; //[0] = number of propsed, [1]=number of accepted

    size_t updatesProposed_;
    size_t updsamespin_;
}; // namespace Markov

template <typename TIOModel, typename TModel>
ABC_MarkovChain<TIOModel, TModel>::~ABC_MarkovChain() {} //destructors must exist

} // namespace Markov
