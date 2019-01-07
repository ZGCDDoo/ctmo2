#pragma once

#include "../Utilities/Utilities.hpp"
#include "../Utilities/MPIUtilities.hpp"
#include "../Utilities/Integrator.hpp"
#include "../Utilities/GreenMat.hpp"
#include "../Utilities/IO.hpp"
#include "ABC_H0.hpp"
#include "UTensorSimple.hpp"
#include "HybFMAndTLoc.hpp"

namespace Models
{

class ABC_Model_2D
{

      public:
        ABC_Model_2D(const Json &jjSim) : ioModel_(jjSim),
                                          h0_(jjSim),
                                          hybFM_(),
                                          tLoc_(),
                                          U_(jjSim["U"].get<double>()),
                                          beta_(jjSim["beta"].get<double>()),
                                          mu_(jjSim["mu"].get<double>()),
                                          NOrb_(jjSim["NOrb"].get<size_t>()),
                                          Nc_(h0_.Nc)
        {
                mpiUt::Print("start abc_model constructor ");

                if (mpiUt::Rank() == mpiUt::master)
                {

#ifdef DCA
                        HybFMAndTLoc<TH0>::CalculateHybFMAndTLoc(h0_);
#else
                        h0_.SaveTKTildeAndHybFM();
#endif
                }

                //tLoc and hybFM should have been calculated by now.

                Conventions::MapSS_t mapNames = Conventions::BuildFileNameConventions();

                assert(tLoc_.load(mapNames["tlocFile"]));
                assert(hybFM_.load(mapNames["hybFMFile"]));

                FinishConstructor(jjSim);
                mpiUt::Print(" End of ABC_Model Constructor ");
        };

        void FinishConstructor(const Json &jjSim)
        {
                std::string hybNameUp = jjSim["HybFileUp"].get<std::string>();
#ifdef DCA
                ClusterCubeCD_t hybtmpUp = ioModel_.ReadGreenKDat(hybNameUp + ".dat", NOrb_);
#else
                ClusterCubeCD_t hybtmpUp = ioModel_.ReadGreenDat(hybNameUp + ".dat", NOrb_);
#endif

#ifdef AFM
                std::string hybNameDown = jjSim["HybFileDown"].get<std::string>();
                ClusterCubeCD_t hybtmpDown = ioModel_.ReadGreenDat(hybNameDown + ".dat", NOrb_);
#endif

                const size_t NHyb = hybtmpUp.n_slices;
                const double factNHyb = 3.0;
                const size_t NHyb_HF = std::max<double>(factNHyb * static_cast<double>(NHyb),
                                                        0.5 * (MIN_EHYB_ * beta_ / M_PI - 1.0));
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

                //this is in fact greencluster tilde.
                const UTensor ut(jjSim);
                this->greenCluster0MatUp_ = GreenMat::GreenCluster0Mat(this->hybridizationMatUp_, this->tLoc_, ut.auxMu(), this->beta_);
#ifdef DCA
                greenCluster0MatUp_.FourierTransform(h0_.RSites(), h0_.KWaveVectors());
#endif
                //save green0mat
                if (mpiUt::Rank() == mpiUt::master)
                {
                        ioModel_.SaveCube("giwn", this->greenCluster0MatUp_.data(), this->beta_, NOrb_);
                }
#ifndef AFM
                this->greenCluster0MatDown_ = greenCluster0MatUp_;
#endif
#ifdef AFM
                this->greenCluster0MatDown_ = GreenMat::GreenCluster0Mat(this->hybridizationMatDown_, this->tLoc_, ut.auxMu(), this->beta_);
#endif
        }

        ~ABC_Model_2D() = default;

        //Getters
        double mu() const { return mu_; }
        double U() const { return U_; }
        double beta() const { return beta_; }
        size_t NOrb() const { return NOrb_; }
        ClusterMatrixCD_t tLoc() const { return tLoc_; }
        GreenMat::GreenCluster0Mat const greenCluster0MatUp() const { return greenCluster0MatUp_; }
        GreenMat::GreenCluster0Mat const greenCluster0MatDown() const { return greenCluster0MatDown_; }
        GreenMat::HybridizationMat const hybridizationMatUp() const { return hybridizationMatUp_; }
        GreenMat::HybridizationMat const hybridizationMatDown() const { return hybridizationMatDown_; }
        Models::ABC_H0 const h0() const { return h0_; }
        IO::Base_IOModel const ioModel() const { return ioModel_; }
        size_t Nc() const { return Nc_; }

        double auxU() const { return U_ / 2.0; }

      protected:
        IO::Base_IOModel ioModel_;
        GreenMat::HybridizationMat hybridizationMatUp_;
        GreenMat::HybridizationMat hybridizationMatDown_;
        GreenMat::GreenCluster0Mat greenCluster0MatUp_;
        GreenMat::GreenCluster0Mat greenCluster0MatDown_;
        Models::ABC_H0 h0_;

        ClusterMatrixCD_t hybFM_;
        ClusterMatrixCD_t tLoc_;

        const double U_;
        const double beta_;
        const double mu_;
        const size_t NOrb_;
        const double MIN_EHYB_ = 300;
        const size_t Nc_;
};

} // namespace Models
