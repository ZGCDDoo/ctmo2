#pragma once
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
typedef LinAlg::Matrix_t Matrix_t;

struct NFData
{

    NFData() : FVup_(), FVdown_(), Nup_(), Ndown_(), dummy_(){};
    SiteVector_t FVup_;
    SiteVector_t FVdown_;
    Matrix_t Nup_;
    Matrix_t Ndown_;
    Matrix_t dummy_;
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
        assert(nfdata_.Nup_.n_rows() + nfdata_.Ndown_.n_rows() == 2 * kk);
        assert(2 * dataCT_->vertices_.size() == dataCT_->vertices_.NUp() + dataCT_->vertices_.NDown());
        // assert(dataCT_->vertices_.NUp() == dataCT_->vertices_.NDown());
        assert(dataCT_->vertices_.NUp() == nfdata_.FVup_.n_elem);
        assert(dataCT_->vertices_.NUp() == nfdata_.Nup_.n_rows());

        assert(dataCT_->vertices_.NDown() == nfdata_.FVdown_.n_elem);
        assert(dataCT_->vertices_.NDown() == nfdata_.Ndown_.n_rows());
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

    void InsertVertexDiffSpin(const Vertex &vertex)
    {
        // std::cout << "Start InsertVertexDiffSpin " << std::endl;

        const double fauxup = FAux(FermionSpin_t::Up, vertex.aux());
        const double fauxdown = FAux(FermionSpin_t::Down, vertex.aux());
        const double fauxupM1 = fauxup - 1.0;
        const double fauxdownM1 = fauxdown - 1.0;

        const VertexPart vertexPartUp = vertex.vStart();
        const VertexPart vertexPartDown = vertex.vEnd();
        assert(vertexPartUp.spin() == FermionSpin_t::Up);
        assert(vertexPartDown.spin() == FermionSpin_t::Down);

        const double sUp = -fauxup + GetGreenTau0(vertexPartUp, vertexPartUp) * fauxupM1;
        const double sDown = -fauxdown + GetGreenTau0(vertexPartDown, vertexPartDown) * fauxdownM1;

        if (dataCT_->vertices_.size())
        {
            AssertSizes();
            const size_t kkoldUp = dataCT_->vertices_.NUp();
            const size_t kknewUp = kkoldUp + 1;
            const size_t kkoldDown = dataCT_->vertices_.NDown();
            const size_t kknewDown = kkoldDown + 1;

            SiteVector_t newLastColUp_(kkoldUp);
            SiteVector_t newLastRowUp_(kkoldUp);
            SiteVector_t newLastColDown_(kkoldDown);
            SiteVector_t newLastRowDown_(kkoldDown);

            //Probably put this in a method
            for (size_t iUp = 0; iUp < dataCT_->vertices_.NUp(); iUp++)
            {
                newLastRowUp_(iUp) = GetGreenTau0(vertexPartUp, dataCT_->vertices_.atUp(iUp)) * (nfdata_.FVup_(iUp) - 1.0);
                newLastColUp_(iUp) = GetGreenTau0(dataCT_->vertices_.atUp(iUp), vertexPartUp) * fauxupM1;
            }

            for (size_t iDown = 0; iDown < dataCT_->vertices_.NDown(); iDown++)
            {
                newLastRowDown_(iDown) = GetGreenTau0(vertexPartDown, dataCT_->vertices_.atDown(iDown)) * (nfdata_.FVdown_(iDown) - 1.0);
                newLastColDown_(iDown) = GetGreenTau0(dataCT_->vertices_.atDown(iDown), vertexPartDown) * fauxdownM1;
            }

            SiteVector_t NQUp(kkoldUp); //NQ = N*Q
            SiteVector_t NQDown(kkoldDown);
            MatrixVectorMult(nfdata_.Nup_, newLastColUp_, 1.0, NQUp);
            MatrixVectorMult(nfdata_.Ndown_, newLastColDown_, 1.0, NQDown);
            const double sTildeUpI = sUp - LinAlg::DotVectors(newLastRowUp_, NQUp);
            const double sTildeDownI = sDown - LinAlg::DotVectors(newLastRowDown_, NQDown);

            const double ratio = sTildeUpI * sTildeDownI;
            const double ratioAcc = PROBREMOVE / PROBINSERT * vertex.probProb() / static_cast<size_t>(dataCT_->vertices_.size() + 1) * ratio;
            if (urng_() < std::abs(ratioAcc))
            {
                updStats_["Inserts"][1]++;
                if (ratioAcc < .0)
                {
                    dataCT_->sign_ *= -1;
                }

                LinAlg::BlockRankOneUpgrade(nfdata_.Nup_, NQUp, newLastRowUp_, 1.0 / sTildeUpI);
                LinAlg::BlockRankOneUpgrade(nfdata_.Ndown_, NQDown, newLastRowDown_, 1.0 / sTildeDownI);
                nfdata_.FVup_.resize(kknewUp);
                nfdata_.FVdown_.resize(kknewDown);
                nfdata_.FVup_(kkoldUp) = fauxup;
                nfdata_.FVdown_(kkoldDown) = fauxdown;
                dataCT_->vertices_.AppendVertex(vertex);
                AssertSizes();
            }
        }
        else
        {
            AssertSizes();
            const double ratioAcc = PROBREMOVE / PROBINSERT * vertex.probProb() * sUp * sDown;

            if (urng_() < std::abs(ratioAcc))
            {
                if (ratioAcc < 0.0)
                {
                    dataCT_->sign_ *= -1;
                }

                nfdata_.Nup_ = Matrix_t(1, 1);
                nfdata_.Ndown_ = Matrix_t(1, 1);
                nfdata_.Nup_(0, 0) = 1.0 / sUp;
                nfdata_.Ndown_(0, 0) = 1.0 / sDown;

                nfdata_.FVup_ = SiteVector_t(1);
                nfdata_.FVdown_ = SiteVector_t(1);
                nfdata_.FVup_(0) = fauxup;
                nfdata_.FVdown_(0) = fauxdown;

                dataCT_->vertices_.AppendVertex(vertex);
            }
        }
        // std::cout << "End InsertVertexDiffSpin " << std::endl;
    }

    void InsertVertexSameSpin(const Vertex &vertex, Matrix_t &Nspin, SiteVector_t &FVspin)
    {
        // std::cout << "Start InsertVertexSameSpin " << std::endl;
        // return;
        AssertSizes();

        const VertexPart x = vertex.vStart();
        const VertexPart y = vertex.vEnd();
        assert(x.spin() == y.spin());

        const double faux = FAux(x.spin(), vertex.aux());
        const double fauxM1 = faux - 1.0;
        const double s00 = -faux + GetGreenTau0(x, x) * fauxM1;
        const double s01 = GetGreenTau0(x, y);
        const double s10 = GetGreenTau0(y, x);
        const double s11 = -faux + GetGreenTau0(y, y) * fauxM1;

        if (dataCT_->vertices_.size())
        {
            AssertSizes();
            const size_t kkold = dataCT_->vertices_.size();
            const size_t kknew = kkold + 1;
            const size_t kkoldspin = Nspin.n_rows();

            Matrix_t Q_(kkoldspin, 2);
            Matrix_t R_(2, kkoldspin);

            for (size_t i = 0; i < kkoldspin; i++)
            {

                const VertexPart vPartI = (x.spin() == FermionSpin_t::Up) ? dataCT_->vertices_.atUp(i) : dataCT_->vertices_.atDown(i);
                const double fauxIm1 = FVspin(i) - 1.0; // Faux_i - 1.0
                Q_(i, 0) = GetGreenTau0(vPartI, x) * fauxM1;
                Q_(i, 1) = GetGreenTau0(vPartI, y) * fauxM1;

                assert(vPartI.spin() == x.spin());

                R_(0, i) = GetGreenTau0(x, vPartI) * fauxIm1;
                R_(1, i) = GetGreenTau0(y, vPartI) * fauxIm1;
            }
            //     // std::cout << "In INsertvertex After loop " << std::endl;

            //     //Watch out, we are calculating two times the matrix NQ, once here and once in ranktwoupgrade. In a next version, only calculate here, not in ranktwoupgrade.
            //     // Matrix_t NQ(2 * kkold, 2); //NQ = N*Q
            //     // MatrixVectorMult(nfdata_.N_, Q_, 1.0, NQUp);

            //     // Matrix_t RNQ(2, 2); //R*NQ
            //     //     // Matrix_t RNQ(2, 2); //R*NQ
            Matrix_t sTilde = Matrix_t({{s00, s01}, {s10, s11}}) - LinAlg::DotRank2(R_, Nspin, Q_);
            sTilde.Inverse();
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
                FVspin.resize(kkoldspin + 2);
                FVspin(kkoldspin) = faux;
                FVspin(kkoldspin + 1) = faux;
                dataCT_->vertices_.AppendVertex(vertex);
                AssertSizes();
                updsamespin_++;
                // std::cout << "InsertVertexSameSpin accepted " << std::endl;
            }
        }
        else
        {
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

                FVspin = SiteVector_t(2);
                FVspin(0) = faux;
                FVspin(1) = faux;

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
        updStats_["Removes"][0]++;
        if (dataCT_->vertices_.size())
        {
            const size_t pp = static_cast<int>(urng_() * dataCT_->vertices_.size());
            const Vertex vertex = dataCT_->vertices_.at(pp);
            const VertexPart x = vertex.vStart();
            const VertexPart y = vertex.vEnd();

            if (x.spin() == y.spin())
            {
                x.spin() == (FermionSpin_t::Up) ? RemoveVertexSameSpin(pp, nfdata_.Nup_, nfdata_.FVup_) : RemoveVertexSameSpin(pp, nfdata_.Ndown_, nfdata_.FVdown_);
            }
            else
            {
                RemoveVertexDiffSpin(pp);
            }
        }
        // std::cout << "End RemoveVertex " << std::endl;
    }

    void RemoveVertexDiffSpin(const size_t &pp)
    {
        // std::cout << "Start RemoveVertexDiffSpin " << std::endl;
        AssertSizes();

        const Vertex vertex = dataCT_->vertices_.at(pp);
        const size_t ppUp = dataCT_->vertices_.GetIndicesSpins(pp, FermionSpin_t::Up);
        const size_t ppDown = dataCT_->vertices_.GetIndicesSpins(pp, FermionSpin_t::Down);
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
            // std::cout << "here 3" << std::endl;

            LinAlg::BlockRankOneDowngrade(nfdata_.Nup_, ppUp);
            LinAlg::BlockRankOneDowngrade(nfdata_.Ndown_, ppDown);

            // std::cout << "here 4" << std::endl;

            nfdata_.FVup_.swap_rows(ppUp, kkUpm1);
            nfdata_.FVdown_.swap_rows(ppDown, kkDownm1);
            nfdata_.FVup_.resize(kkUpm1);
            nfdata_.FVdown_.resize(kkDownm1);

            // std::cout << "here 5" << std::endl;

            dataCT_->vertices_.RemoveVertex(pp);
            AssertSizes();
            // std::cout << "here 6" << std::endl;
        }
        // std::cout << "End RemoveVertexDiffSpin " << std::endl;
    }

    void RemoveVertexSameSpin(const size_t &pp, Matrix_t &Nspin, SiteVector_t &FVspin)
    {
        // std::cout << "Start RemoveVertexSameSpin " << std::endl;
        return;
        AssertSizes();
        assert(Nspin.n_rows() >= 2);
        assert(FVspin.n_elem >= 2);

        const Vertex vertex = dataCT_->vertices_.at(pp);
        const VertexPart x = vertex.vStart();
        const VertexPart y = vertex.vEnd();
        assert(x.spin() == y.spin());

        const size_t ppSpin = dataCT_->vertices_.GetIndicesSpins(pp, x.spin());

        const size_t kk = dataCT_->vertices_.size();
        const size_t kkSpin = (x.spin() == FermionSpin_t::Up) ? dataCT_->vertices_.NUp() : dataCT_->vertices_.NDown();
        const size_t kkSpinm2 = kkSpin - 2;

        const ClusterMatrix_t STildeInverse = {{Nspin(ppSpin, ppSpin), Nspin(ppSpin, ppSpin + 1)}, {Nspin(ppSpin + 1, ppSpin), Nspin(ppSpin + 1, ppSpin + 1)}};
        const double ratioAcc = PROBINSERT / PROBREMOVE * static_cast<double>(kk) / vertex.probProb() * arma::det(STildeInverse);

        if (urng_() < std::abs(ratioAcc))
        {
            AssertSizes();
            updStats_["Removes"][1]++;
            if (ratioAcc < 0.0)
            {
                dataCT_->sign_ *= -1;
            }

            LinAlg::BlockDowngrade(Nspin, ppSpin, 2);

            FVspin.swap_rows(ppSpin, kkSpin - 2);
            FVspin.swap_rows(ppSpin + 1, kkSpin - 1);
            FVspin.resize(kkSpinm2);

            dataCT_->vertices_.RemoveVertex(pp);
            // std::cout << "RemoveVertexSameSpin accepted " << std::endl;
        }

        AssertSizes();
        // std::cout << "End RemoveVertexSameSpin " << std::endl;
    }

    void CleanUpdate()
    {
        mpiUt::Print("Cleaning, sign, k =  " + std::to_string(dataCT_->sign_) + ",  " + std::to_string(dataCT_->vertices_.size()));
        AssertSizes();
        std::cout << "updsamespin = " << updsamespin_ << std::endl;

        const size_t kkup = dataCT_->vertices_.NUp();
        const size_t kkdown = dataCT_->vertices_.NDown();

        if (kkup != 0)
        {

            for (size_t iUp = 0; iUp < kkup; iUp++)
            {
                for (size_t jUp = 0; jUp < kkup; jUp++)
                {

                    nfdata_.Nup_(iUp, jUp) = GetGreenTau0(dataCT_->vertices_.atUp(iUp), dataCT_->vertices_.atUp(jUp)) * (nfdata_.FVup_(jUp) - 1.0);

                    if (iUp == jUp)
                    {
                        nfdata_.Nup_(iUp, iUp) -= nfdata_.FVup_(iUp);
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

                    nfdata_.Ndown_(iDown, jDown) = GetGreenTau0(dataCT_->vertices_.atDown(iDown), dataCT_->vertices_.atDown(jDown)) * (nfdata_.FVdown_(jDown) - 1.0);

                    if (iDown == jDown)
                    {
                        nfdata_.Ndown_(iDown, iDown) -= nfdata_.FVdown_(iDown);
                    }
                }
            }
            nfdata_.Ndown_.Inverse();
        }

        std::cout << "End cleaning " << std::endl;
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
        const SiteVector_t FVupM1 = (nfdata_.FVup_ - 1.0);
        const SiteVector_t FVdownM1 = (nfdata_.FVdown_ - 1.0);
        DDMGMM(FVupM1, nfdata_.Nup_, *(dataCT_->MupPtr_));
        DDMGMM(FVdownM1, nfdata_.Ndown_, *(dataCT_->MdownPtr_));
        obs_.Measure();
    }

    void SaveMeas()
    {

        obs_.Save();
        // mpiUt::SaveConfig(dataCT_->vertices_);
        SaveUpd("upd.meas");
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
