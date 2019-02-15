#pragma once

#include "ctmo/Foundations/Conventions.hpp"
#include "ctmo/Foundations/GreenMat.hpp"
#include "ctmo/Foundations/IO.hpp"
#include "ctmo/Model/ABC_H0.hpp"
#include "ctmo/Model/UTensorSimple.hpp"

#ifdef DCA
#include "ctmo/Model/HybFMAndTLoc.hpp"
#endif

namespace Models
{

class ABC_Model_2D
{

  public:
    explicit ABC_Model_2D(const Json &jjSim)
        : ioModelPtr_(new IO::Base_IOModel(jjSim)), h0_(jjSim), hybFM_(), tLoc_(), beta_(jjSim["model"]["beta"].get<double>()),
          mu_(jjSim["model"]["mu"].get<double>()), NOrb_(jjSim["model"]["nOrb"].get<size_t>()), Nc_(h0_.Nc)
    {
        Logging::Debug("Start ABC_Model Constructor. ");

        if (mpiUt::Tools::Rank() == mpiUt::Tools::master)
        {

#ifdef DCA
            HybFMAndTLoc::CalculateHybFMAndTLoc(h0_);
#else
            h0_.SaveTKTildeAndHybFM();
#endif
        }

        // tLoc and hybFM should have been calculated by now.

        Conventions::MapSS_t mapNames = Conventions::BuildFileNameConventions();

        assert(tLoc_.load(mapNames.at("tlocFile")));
        assert(hybFM_.load(mapNames.at("hybFMFile")));

        FinishConstructor(jjSim);
        Logging::Debug(" End of ABC_Model Constructor. ");
    };

    void FinishConstructor(const Json &jjSim)
    {
        std::string hybNameUp = jjSim["model"]["hybUpFile"].get<std::string>();
#ifdef DCA
        ClusterCubeCD_t hybtmpUp = ioModelPtr_->ReadGreenKDat(hybNameUp, NOrb_);
#else
        ClusterCubeCD_t hybtmpUp = ioModelPtr_->ReadGreenDat(hybNameUp, NOrb_);
#endif

#ifdef AFM
        std::string hybNameDown = jjSim["model"]["hybDownFile"].get<std::string>();
        ClusterCubeCD_t hybtmpDown = ioModelPtr_->ReadGreenDat(hybNameDown, NOrb_);
#endif

        const size_t NHyb = hybtmpUp.n_slices;
        const double factNHyb = 3.0;
        const size_t NHyb_HF = std::max<double>(factNHyb * static_cast<double>(NHyb), 0.5 * (MIN_EHYB_ * beta_ / M_PI - 1.0));
        hybtmpUp.resize(Nc_ * NOrb_, Nc_ * NOrb_, NHyb_HF);
#ifdef AFM
        hybtmpDown.resize(Nc_ * NOrb_, Nc_ * NOrb_, NHyb_HF);
#endif

        for (size_t nn = NHyb; nn < NHyb_HF; nn++)
        {
            const cd_t iwn(0.0, (2.0 * nn + 1.0) * M_PI / beta_);
            hybtmpUp.slice(nn) = hybFM_ / iwn;
#ifdef AFM
            hybtmpDown.slice(nn) = hybtmpUp.slice(nn);
#endif
        }

        this->hybridizationMatUp_ = GreenMat::HybridizationMat(hybtmpUp, this->hybFM_);
#ifdef AFM
        this->hybridizationMatDown_ = GreenMat::HybridizationMat(hybtmpDown, this->hybFM_);
#endif

        // this is in fact greencluster tilde.
        const UTensor ut(jjSim);
        this->greenCluster0MatUp_ = GreenMat::GreenCluster0Mat(this->hybridizationMatUp_, this->tLoc_, ut.auxMu(), this->beta_);
#ifdef DCA
        greenCluster0MatUp_.FourierTransform(h0_.RSites(), h0_.KWaveVectors());
#endif
        // save green0mat
        if (mpiUt::Tools::Rank() == mpiUt::Tools::master)
        {
            ioModelPtr_->SaveCube("giwn", this->greenCluster0MatUp_.data(), this->beta_, NOrb_);
        }
#ifndef AFM
        this->greenCluster0MatDown_ = greenCluster0MatUp_;
#endif
#ifdef AFM
        this->greenCluster0MatDown_ = GreenMat::GreenCluster0Mat(this->hybridizationMatDown_, this->tLoc_, ut.auxMu(), this->beta_);
#endif
    }

    ~ABC_Model_2D() = default;

    // Getters
    double mu() const { return mu_; }
    double beta() const { return beta_; }
    size_t NOrb() const { return NOrb_; }
    ClusterMatrixCD_t tLoc() const { return tLoc_; }
    GreenMat::GreenCluster0Mat const greenCluster0MatUp() const { return greenCluster0MatUp_; }
    GreenMat::GreenCluster0Mat const greenCluster0MatDown() const { return greenCluster0MatDown_; }
    GreenMat::HybridizationMat const hybridizationMatUp() const { return hybridizationMatUp_; }
    GreenMat::HybridizationMat const hybridizationMatDown() const { return hybridizationMatDown_; }
    Models::ABC_H0 const h0() const { return h0_; }
    const std::shared_ptr<IO::Base_IOModel> ioModelPtr() const { return ioModelPtr_; }
    size_t Nc() const { return Nc_; }

  protected:
    std::shared_ptr<IO::Base_IOModel> ioModelPtr_;
    GreenMat::HybridizationMat hybridizationMatUp_;
    GreenMat::HybridizationMat hybridizationMatDown_;
    GreenMat::GreenCluster0Mat greenCluster0MatUp_;
    GreenMat::GreenCluster0Mat greenCluster0MatDown_;
    Models::ABC_H0 h0_;

    ClusterMatrixCD_t hybFM_;
    ClusterMatrixCD_t tLoc_;

    const double beta_;
    const double mu_;
    const size_t NOrb_;
    const double MIN_EHYB_ = 300;
    const size_t Nc_;
};

} // namespace Models
