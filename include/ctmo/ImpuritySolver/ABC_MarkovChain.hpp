#pragma once

#include "ctmo/Foundations/UtilitiesRandom.hpp"
#include "ctmo/Foundations/Matrix.hpp"
#include "ctmo/Foundations/GreenTau.hpp"

#ifdef SLMC
#include "ctmo/ImpuritySolver/ConfigParser.hpp"
#endif

#include "ctmo/ImpuritySolver/Observables.hpp"

namespace Markov
{

using Fourier::MatToTau;
using Fourier::MatToTauCluster;
using Vertex = Diagrammatic::Vertex;
using VertexPart = Diagrammatic::VertexPart;

using Matrix_t = LinAlg::Matrix_t;

struct UpdData
{

    UpdData() = default;
    SiteVector_t NQUp_;
    SiteVector_t NQDown_;
    SiteVector_t newLastRowUp_;
    SiteVector_t newLastRowDown_;
    double sTildeUpI_{0.0};
    double sTildeDownI_{0.0};
};

struct NFData
{

    NFData() = default;
    SiteVector_t FVup_;
    SiteVector_t FVdown_;
    Matrix_t Nup_;
    Matrix_t Ndown_;
    Matrix_t dummy_;
};

class ABC_MarkovChain
{

    using GreenTau_t = GreenTau::GreenCluster0Tau;
    using Model_t = Models::ABC_Model_2D;

  public:
    const double PROBINSERT = 0.3333333333;
    const double PROBREMOVE = 1.0 - PROBINSERT;

    ABC_MarkovChain(const Json &jj, const size_t &seed)
        : modelPtr_(new Model_t(jj)), rng_(seed), urng_(rng_, Utilities::UniformDistribution_t(0.0, 1.0)),
          dataCT_(new Obs::ISDataCT(jj, modelPtr_)), obs_(dataCT_, jj), vertexBuilder_(jj, modelPtr_->Nc()),
#ifdef SLMC
          configParser_(), logDeterminant_(0.0),
#endif
          updsamespin_(0), isOneOrbitalOptimized_(jj["solver"]["isOneOrbitalOptimized"].get<bool>())
    {
        const std::valarray<size_t> zeroPair = {0, 0};
        updStats_["Inserts"] = zeroPair;
        updStats_["Removes"] = zeroPair;
        updStats_["Flips"] = zeroPair;
        updatesProposed_ = 0;

#ifdef SLMC
        if (!isOneOrbitalOptimized_)
        {
            throw std::runtime_error("Should be one Orbital optimized for SLMC.");
        }
#endif
        if (isOneOrbitalOptimized_)
        {
            Logging::Trace("Optimized for one orbital and will crash if not One-band Hubbard Model.");
        }
        Logging::Debug("MarkovChain Created.");
    }

    ABC_MarkovChain(const ABC_MarkovChain &abc_markovchain) = default;
    ABC_MarkovChain &operator=(const ABC_MarkovChain &abc_markovchain) = delete;
    ABC_MarkovChain(ABC_MarkovChain &&abc_markovChain) = default;
    ABC_MarkovChain &operator=(ABC_MarkovChain &&abc_markovChain) = delete;

    virtual ~ABC_MarkovChain() = default;

    // Getters
    Model_t model() const { return (*modelPtr_); }

    Matrix_t Nup() const { return nfdata_.Nup_; }

    Matrix_t Ndown() const { return nfdata_.Ndown_; }

    size_t updatesProposed() const { return updatesProposed_; }

    double beta() const { return dataCT_->beta_; }

#ifdef SLMC
    double logDeterminant() const { return logDeterminant_; }
#endif

    // End Getters
    virtual double FAux(const VertexPart &vPart) const = 0;

    virtual double FAuxBar(const VertexPart &vPart) const = 0;

    virtual double gamma(const VertexPart &vpI, const VertexPart &vpJ) const = 0; // little gamma

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

        assert(dataCT_->vertices_.NUp() == nfdata_.FVup_.n_elem);
        assert(dataCT_->vertices_.NUp() == nfdata_.Nup_.n_rows());

        assert(dataCT_->vertices_.NDown() == nfdata_.FVdown_.n_elem);
        assert(dataCT_->vertices_.NDown() == nfdata_.Ndown_.n_rows());
    }

    void InsertVertex()
    {
        AssertSizes();
        updStats_["Inserts"][0]++;
        const Vertex vertex = vertexBuilder_.BuildVertex(urng_);
        const VertexPart x = vertex.vStart();
        const VertexPart y = vertex.vEnd();

        if (x.spin() == y.spin())
        {
            x.spin() == (FermionSpin_t::Up) ? InsertVertexSameSpin(vertex, nfdata_.Nup_, nfdata_.FVup_)
                                            : InsertVertexSameSpin(vertex, nfdata_.Ndown_, nfdata_.FVdown_);
        }
        else
        {
            InsertVertexDiffSpin(vertex);
        }
    }

    void CalculateUpdDataInsertDiffSpin(const VertexPart &x)
    {
        const double faux = FAux(x);
        const double fauxM1 = faux - 1.0;

        if (x.spin() == FermionSpin_t::Up)
        {
            upddata_.sTildeUpI_ = -faux + GetGreenTau0(x, x) * fauxM1;

            if (static_cast<bool>(nfdata_.Nup_.n_rows()))
            {
                const size_t kkoldUp = dataCT_->vertices_.NUp();
                upddata_.newLastRowUp_.set_size(kkoldUp);
                SiteVector_t newLastColUp_(kkoldUp);
                upddata_.NQUp_.set_size(kkoldUp);

                for (size_t iUp = 0; iUp < dataCT_->vertices_.NUp(); iUp++)
                {
                    upddata_.newLastRowUp_(iUp) = GetGreenTau0(x, dataCT_->vertices_.atUp(iUp)) * (nfdata_.FVup_(iUp) - 1.0);
                    newLastColUp_(iUp) = GetGreenTau0(dataCT_->vertices_.atUp(iUp), x) * fauxM1;
                }
                MatrixVectorMult(nfdata_.Nup_, newLastColUp_, 1.0, upddata_.NQUp_);
                upddata_.sTildeUpI_ -= LinAlg::DotVectors(upddata_.newLastRowUp_, upddata_.NQUp_);
            }
        }
        else
        {
            upddata_.sTildeDownI_ = -faux + GetGreenTau0(x, x) * fauxM1;

            if (static_cast<bool>(nfdata_.Ndown_.n_rows()))
            {
                const size_t kkoldDown = dataCT_->vertices_.NDown();
                upddata_.newLastRowDown_.set_size(kkoldDown);
                SiteVector_t newLastColDown_(kkoldDown);
                upddata_.NQDown_.set_size(kkoldDown);

                for (size_t iDown = 0; iDown < dataCT_->vertices_.NDown(); iDown++)
                {
                    upddata_.newLastRowDown_(iDown) = GetGreenTau0(x, dataCT_->vertices_.atDown(iDown)) * (nfdata_.FVdown_(iDown) - 1.0);
                    newLastColDown_(iDown) = GetGreenTau0(dataCT_->vertices_.atDown(iDown), x) * fauxM1;
                }

                MatrixVectorMult(nfdata_.Ndown_, newLastColDown_, 1.0, upddata_.NQDown_);
                upddata_.sTildeDownI_ -= LinAlg::DotVectors(upddata_.newLastRowDown_, upddata_.NQDown_);
            }
        }
    }

    void InsertVertexDiffSpin(const Vertex &vertex)
    {

        const VertexPart vertexPartUp = vertex.vStart();
        const VertexPart vertexPartDown = vertex.vEnd();
        const size_t kkoldUp = nfdata_.Nup_.n_rows();
        const size_t kknewUp = kkoldUp + 1;
        const size_t kkoldDown = nfdata_.Ndown_.n_rows();
        const size_t kknewDown = kkoldDown + 1;
        const double fauxup = FAux(vertexPartUp);
        const double fauxdown = FAux(vertexPartDown);

        CalculateUpdDataInsertDiffSpin(vertexPartUp);
        CalculateUpdDataInsertDiffSpin(vertexPartDown);

        const double ratio = upddata_.sTildeUpI_ * upddata_.sTildeDownI_;
        const double ratioAcc = PROBREMOVE / PROBINSERT * vertex.probProb() / static_cast<size_t>(dataCT_->vertices_.size() + 1) * ratio;

        if (urng_() < std::abs(ratioAcc))
        {
            AssertSizes();
#ifdef SLMC
            logDeterminant_ += std::log(std::abs(ratio));

#endif
            updStats_["Inserts"][1]++;
            if (ratioAcc < 0.0)
            {
                dataCT_->sign_ *= -1;
            }

            if (static_cast<bool>(nfdata_.Nup_.n_rows()))
            {
                LinAlg::BlockRankOneUpgrade(nfdata_.Nup_, upddata_.NQUp_, upddata_.newLastRowUp_, 1.0 / upddata_.sTildeUpI_);
                nfdata_.FVup_.resize(kknewUp);
                nfdata_.FVup_(kkoldUp) = fauxup;
            }
            else
            {
                nfdata_.Nup_ = Matrix_t(1, 1);
                nfdata_.Nup_(0, 0) = 1.0 / upddata_.sTildeUpI_;
                nfdata_.FVup_ = SiteVector_t(1);
                nfdata_.FVup_(0) = fauxup;
            }

            if (static_cast<bool>(nfdata_.Ndown_.n_rows()))
            {
                LinAlg::BlockRankOneUpgrade(nfdata_.Ndown_, upddata_.NQDown_, upddata_.newLastRowDown_, 1.0 / upddata_.sTildeDownI_);
                nfdata_.FVdown_.resize(kknewDown);
                nfdata_.FVdown_(kkoldDown) = fauxdown;
            }
            else
            {
                nfdata_.Ndown_ = Matrix_t(1, 1);
                nfdata_.Ndown_(0, 0) = 1.0 / upddata_.sTildeDownI_;

                nfdata_.FVdown_ = SiteVector_t(1);
                nfdata_.FVdown_(0) = fauxdown;
            }

            dataCT_->vertices_.AppendVertex(vertex);
        }
    }

    void InsertVertexSameSpin(const Vertex &vertex, Matrix_t &Nspin, SiteVector_t &FVspin)
    {
        const VertexPart x = vertex.vStart();
        const VertexPart y = vertex.vEnd();

        if ((std::abs(x.tau() - y.tau()) < 1e-14) && (x.vtype() != Diagrammatic::VertexType::Phonon))
        {
            assert(x.orbital() != y.orbital());
        }

        assert(x.spin() == y.spin());

        const double faux = FAux(x);
        const double fauxM1 = faux - 1.0;
        const double fauxBar = FAuxBar(x);
        const double fauxM1Bar = fauxBar - 1.0;

        const double s00 = -faux + GetGreenTau0(x, x) * fauxM1;
        const double s01 = GetGreenTau0(x, y) * fauxM1Bar;
        const double s10 = GetGreenTau0(y, x) * fauxM1;
        const double s11 = -fauxBar + GetGreenTau0(y, y) * fauxM1Bar;

        if (static_cast<bool>(Nspin.n_rows()))
        {
            assert(Nspin.n_rows());
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
                Q_(i, 1) = GetGreenTau0(vPartI, y) * fauxM1Bar;

                assert(vPartI.spin() == x.spin());

                R_(0, i) = GetGreenTau0(x, vPartI) * fauxIm1;
                R_(1, i) = GetGreenTau0(y, vPartI) * fauxIm1;
            }

            Matrix_t NQ_(kkoldspin, 2); // NQ = N*Q
            NQ_.Zeros();
            DGEMM(1.0, 0.0, Nspin, Q_, NQ_);

            Matrix_t sTilde(2, 2);
            DGEMM(-1.0, 0.0, R_, NQ_, sTilde);
            sTilde += Matrix_t({{s00, s01}, {s10, s11}});
            sTilde.Inverse();

            const double ratioAcc = PROBREMOVE / PROBINSERT * vertex.probProb() / kknew * 1.0 / sTilde.Determinant();
            if (urng_() < std::abs(ratioAcc))
            {
                AssertSizes();

                updStats_["Inserts"][1]++;
                if (ratioAcc < .0)
                {
                    dataCT_->sign_ *= -1;
                }

                LinAlg::BlockRankTwoUpgrade(Nspin, NQ_, R_, sTilde);
                FVspin.resize(kkoldspin + 2);
                FVspin(kkoldspin) = faux;
                FVspin(kkoldspin + 1) = fauxBar;
                dataCT_->vertices_.AppendVertex(vertex);
                updsamespin_++;
            }
        }
        else
        {
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
                FVspin(1) = fauxBar;

                dataCT_->vertices_.AppendVertex(vertex);
            }
        }
    }

    void RemoveVertex()
    {
        AssertSizes();
        const size_t kk = dataCT_->vertices_.size();
        if (static_cast<bool>(kk))
        {
            updStats_["Removes"][0]++;

            const auto pp = static_cast<size_t>(urng_() * dataCT_->vertices_.size());
            const Vertex vertex = dataCT_->vertices_.at(pp);
            const VertexPart x = vertex.vStart();
            const VertexPart y = vertex.vEnd();

            if (x.spin() == y.spin())
            {
                x.spin() == (FermionSpin_t::Up) ? RemoveVertexSameSpin(pp, nfdata_.Nup_, nfdata_.FVup_)
                                                : RemoveVertexSameSpin(pp, nfdata_.Ndown_, nfdata_.FVdown_);
            }
            else
            {
                RemoveVertexDiffSpin(pp);
            }
        }
    }

    void RemoveVertexDiffSpin(const size_t &pp)
    {

        const Vertex vertex = dataCT_->vertices_.at(pp);
        size_t ppUp = pp;
        size_t ppDown = pp;

        if (!isOneOrbitalOptimized_)
        {
            const UInt64_t vertexKey = dataCT_->vertices_.GetKey(pp);
            ppUp = dataCT_->vertices_.GetKeyIndex(vertexKey, FermionSpin_t::Up);
            ppDown = dataCT_->vertices_.GetKeyIndex(vertexKey, FermionSpin_t::Down);
        }

        const auto x = dataCT_->vertices_.atUp(ppUp);

        const auto y = dataCT_->vertices_.atDown(ppDown);
        assert(x.spin() == FermionSpin_t::Up);
        assert(y.spin() == FermionSpin_t::Down);
        if (x.vtype() != Diagrammatic::VertexType::Phonon)
        {
            assert(x.vtype() == y.vtype());
            assert(std::abs(x.tau() - y.tau()) < 1e-14);
        }

        const double ratio = nfdata_.Nup_(ppUp, ppUp) * nfdata_.Ndown_(ppDown, ppDown);
        const double ratioAcc = PROBINSERT / PROBREMOVE * static_cast<double>(dataCT_->vertices_.size()) / vertex.probProb() * ratio;

        if (urng_() < std::abs(ratioAcc))
        {

#ifdef SLMC
            logDeterminant_ += std::log(std::abs(ratio));

#endif
            updStats_["Removes"][1]++;
            if (ratioAcc < .0)
            {
                dataCT_->sign_ *= -1;
            }

            // The update matrices of size k-1 x k-1 with the pp row and col deleted and the last row and col now at index pp

            const size_t kkUp = dataCT_->vertices_.NUp();
            const size_t kkUpm1 = kkUp - 1;
            const size_t kkDown = dataCT_->vertices_.NDown();
            const size_t kkDownm1 = kkDown - 1;
            dataCT_->vertices_.SwapVertexPart(ppUp, kkUpm1, x.spin());
            dataCT_->vertices_.SwapVertexPart(ppDown, kkDownm1, y.spin());

            LinAlg::BlockRankOneDowngrade(nfdata_.Nup_, ppUp);
            LinAlg::BlockRankOneDowngrade(nfdata_.Ndown_, ppDown);
            nfdata_.FVup_.swap_rows(ppUp, kkUpm1);
            nfdata_.FVdown_.swap_rows(ppDown, kkDownm1);
            nfdata_.FVup_.resize(kkUpm1);
            nfdata_.FVdown_.resize(kkDownm1);

            dataCT_->vertices_.RemoveVertex(pp);
            dataCT_->vertices_.PopBackVertexPart(x.spin());
            dataCT_->vertices_.PopBackVertexPart(y.spin());

            AssertSizes();
        }
    }

    void RemoveVertexSameSpin(const size_t &pp, Matrix_t &Nspin, SiteVector_t &FVspin)
    {
        assert(Nspin.n_rows() >= 2);
        assert(FVspin.n_elem >= 2);

        const Vertex vertex = dataCT_->vertices_.at(pp);
        const UInt64_t vertexKey = dataCT_->vertices_.GetKey(pp);
        const VertexPart x = vertex.vStart();
        const VertexPart y = vertex.vEnd();
        assert(x.spin() == y.spin());

        if (x.vtype() != Diagrammatic::VertexType::Phonon)
        {
            assert(x.vtype() == y.vtype());
            assert(std::abs(x.tau() - y.tau()) < 1e-14);
            assert(x.orbital() != y.orbital());
        }

        assert(x.site() == y.site());

        const size_t pp1Spin = dataCT_->vertices_.GetKeyIndex(vertexKey, x.spin());
        const size_t pp2Spin = dataCT_->vertices_.GetKeyIndex(vertexKey + 1, y.spin());
        const size_t kk = dataCT_->vertices_.size();
        const size_t kkSpin = (x.spin() == FermionSpin_t::Up) ? dataCT_->vertices_.NUp() : dataCT_->vertices_.NDown();
        const size_t kkSpinm1 = kkSpin - 1;
        const size_t kkSpinm2 = kkSpin - 2;

        const ClusterMatrix_t STildeInverse = {{Nspin(pp1Spin, pp1Spin), Nspin(pp1Spin, pp2Spin)},
                                               {Nspin(pp2Spin, pp1Spin), Nspin(pp2Spin, pp2Spin)}};
        const double ratioAcc = PROBINSERT / PROBREMOVE * static_cast<double>(kk) / vertex.probProb() * arma::det(STildeInverse);

        if (urng_() < std::abs(ratioAcc))
        {
            AssertSizes();
            updStats_["Removes"][1]++;
            if (ratioAcc < 0.0)
            {
                dataCT_->sign_ *= -1;
            }

            dataCT_->vertices_.SwapVertexPart(pp2Spin, kkSpinm1, x.spin());
            FVspin.swap_rows(pp2Spin, kkSpinm1);
            Nspin.SwapRowsAndCols(pp2Spin, kkSpinm1);

            const size_t pp1SpinNew = dataCT_->vertices_.GetKeyIndex(vertexKey, y.spin());
            dataCT_->vertices_.SwapVertexPart(pp1SpinNew, kkSpinm2, y.spin());
            FVspin.swap_rows(pp1SpinNew, kkSpinm2);
            Nspin.SwapRowsAndCols(pp1SpinNew, kkSpinm2);

            LinAlg::BlockRankTwoDowngrade(Nspin);

            FVspin.resize(kkSpinm2);

            dataCT_->vertices_.RemoveVertex(pp);
            dataCT_->vertices_.PopBackVertexPart(x.spin());
            dataCT_->vertices_.PopBackVertexPart(y.spin());

            assert(Nspin.n_rows() == FVspin.n_elem);
        }
    }

    void CleanUpdate()
    {

        const size_t kkup = dataCT_->vertices_.NUp();
        const size_t kkdown = dataCT_->vertices_.NDown();

        if (kkup != 0)
        {

            for (size_t iUp = 0; iUp < kkup; iUp++)
            {
                for (size_t jUp = 0; jUp < kkup; jUp++)
                {

                    nfdata_.Nup_(iUp, jUp) =
                        GetGreenTau0(dataCT_->vertices_.atUp(iUp), dataCT_->vertices_.atUp(jUp)) * (nfdata_.FVup_(jUp) - 1.0);

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

                    nfdata_.Ndown_(iDown, jDown) =
                        GetGreenTau0(dataCT_->vertices_.atDown(iDown), dataCT_->vertices_.atDown(jDown)) * (nfdata_.FVdown_(jDown) - 1.0);

                    if (iDown == jDown)
                    {
                        nfdata_.Ndown_(iDown, iDown) -= nfdata_.FVdown_(iDown);
                    }
                }
            }
            nfdata_.Ndown_.Inverse();
        }
    }

    double GetGreenTau0(const VertexPart &x, const VertexPart &y) const
    {
        assert(x.spin() == y.spin());
#ifndef AFM
        return (dataCT_->green0CachedUp_(x.superSite(), y.superSite(), x.tau() - y.tau()));
#else
        if (x.spin() == FermionSpin_t::Up)
        {
            return (dataCT_->green0CachedUp_(x.superSite(), y.superSite(), x.tau() - y.tau()));
        }
        else
        {
            return (dataCT_->green0CachedDown_(x.superSite(), y.superSite(), x.tau() - y.tau()));
        }
#endif
    }

    void Measure()
    {
        AssertSizes();
#ifdef SLMC
        configParser_.SaveConfig(dataCT_->vertices_, logDeterminant_, dataCT_->sign_);
#else
        const SiteVector_t FVupM1 = (nfdata_.FVup_ - 1.0);
        const SiteVector_t FVdownM1 = (nfdata_.FVdown_ - 1.0);
        DDMGMM(FVupM1, nfdata_.Nup_, *(dataCT_->MupPtr_));
        DDMGMM(FVdownM1, nfdata_.Ndown_, *(dataCT_->MdownPtr_));
        obs_.Measure();
#endif
    }

    void SaveMeas()
    {

        obs_.Save();
        Logging::Trace("updsamespin = " + std::to_string(updsamespin_));
        SaveUpd("Measurements");
        if (mpiUt::Tools::Rank() == mpiUt::Tools::master)
        {
            dataCT_->vertices_.SaveConfig("Config.dat");
        }
        Logging::Info("Finished Saving MarkovChain.");
    }

    void SaveTherm()
    {

        SaveUpd("Thermalization");
        for (auto &updStat : updStats_)
        {
            updStat.second = 0;
        }
    }

    void SaveUpd(const std::string &updType)
    {
        std::vector<UpdStats_t> updStatsVec;
#ifdef HAVEMPI

        mpi::communicator world;
        if (mpiUt::Tools::Rank() == mpiUt::Tools::master)
        {
            mpi::gather(world, updStats_, updStatsVec, mpiUt::Tools::master);
        }
        else
        {
            mpi::gather(world, updStats_, mpiUt::Tools::master);
        }
        if (mpiUt::Tools::Rank() == mpiUt::Tools::master)
        {
            Logging::Info("\n\n Statistics Updates of " + updType + ":\n" + mpiUt::Tools::SaveUpdStats(updStatsVec) + "\n\n");
        }

#else
        updStatsVec.push_back(updStats_);
        Logging::Info("\n\n Statistics Updates of " + updType + ":\n" + mpiUt::Tools::SaveUpdStats(updStatsVec) + "\n\n");

#endif
    }

  protected:
    // attributes
    std::shared_ptr<Model_t> modelPtr_;
    Utilities::EngineTypeMt19937_t rng_;
    Utilities::UniformRngMt19937_t urng_;
    NFData nfdata_;
    UpdData upddata_;
    std::shared_ptr<Obs::ISDataCT> dataCT_;
    Obs::Observables obs_;
    Diagrammatic::VertexBuilder vertexBuilder_;

#ifdef SLMC
    Diagrammatic::ConfigParser configParser_;
    double logDeterminant_;

#endif
    UpdStats_t updStats_; //[0] = number of propsed, [1]=number of accepted

    size_t updatesProposed_;
    size_t updsamespin_;

    const bool isOneOrbitalOptimized_;
}; // namespace Markov

} // namespace Markov
