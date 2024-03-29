#pragma once

#include "ctmo/ImpuritySolver/MPIResult.hpp"
#include "ctmo/ImpuritySolver/GreenBinning.hpp"
#include "ctmo/ImpuritySolver/FillingAndDocc.hpp"
#include "ctmo/ImpuritySolver/KineticEnergy.hpp"
#include "ctmo/ImpuritySolver/ISResult.hpp"

namespace Markov
{
namespace Obs
{

using Matrix_t = LinAlg::Matrix<double>;

class Observables
{

  public:
    // Observables(){};
    Observables(std::shared_ptr<ISDataCT> dataCT, const Json &jjSim)
        : dataCT_(std::move(dataCT)), modelPtr_(dataCT_->modelPtr_), ioModelPtr_(modelPtr_->ioModelPtr()),
          rng_(jjSim["monteCarlo"]["seed"].get<size_t>() + mpiUt::Tools::Rank() * mpiUt::Tools::Rank()),
          urngPtr_(new Utilities::UniformRngFibonacci3217_t(rng_, Utilities::UniformDistribution_t(0.0, 1.0))),
          greenBinningUp_(dataCT_, jjSim, FermionSpin_t::Up), greenBinningDown_(dataCT_, jjSim, FermionSpin_t::Down),
          fillingAndDocc_(dataCT_, urngPtr_, jjSim["solver"]["n_tau_sampling"].get<size_t>()), signMeas_(0.0), expOrder_(0.0), NMeas_(0),
          NOrb_(jjSim["model"]["nOrb"].get<size_t>()), averageOrbitals_(jjSim["solver"]["averageOrbitals"].get<bool>())
    {

        Logging::Debug("In Obs constructor ");

        Logging::Debug("After Obs  constructor ");
    }

    // Getters
    double signMeas() const { return signMeas_; };
    double expOrder() const { return expOrder_; };

    void Measure()
    {

        ++NMeas_;
        signMeas_ += static_cast<double>(dataCT_->sign_);
        expOrder_ += static_cast<double>(dataCT_->vertices_.size()) * static_cast<double>(dataCT_->sign_);

        fillingAndDocc_.MeasureFillingAndDocc();

        greenBinningUp_.MeasureGreenBinning(*dataCT_->MupPtr_);
        greenBinningDown_.MeasureGreenBinning(*dataCT_->MdownPtr_);
    }

    void Save()
    {
        Logging::Info("Start of Observables.Save()");
        signMeas_ /= NMeas_;

        fillingAndDocc_.Finalize(signMeas_, NMeas_);
        std::map<std::string, double> obsScal;

        obsScal = fillingAndDocc_.GetObs();

        obsScal["sign"] = signMeas_;
        obsScal["NMeas"] = NMeas_;

        // dont forget that the following obs have not been finalized (multiplied by following factor)
        const double fact = 1.0 / (NMeas_ * signMeas_);
        obsScal["k"] = fact * expOrder_;

        ClusterCubeCD_t greenCubeMatUp = greenBinningUp_.FinalizeGreenBinning(signMeas_, NMeas_);
        ClusterCubeCD_t greenCubeMatDown = greenBinningDown_.FinalizeGreenBinning(signMeas_, NMeas_);

        // Average the green Functions if orbitals have the same parameters
        if (averageOrbitals_)
        {
            greenCubeMatUp = ioModelPtr_->AverageOrbitals(greenCubeMatUp);
            greenCubeMatDown = ioModelPtr_->AverageOrbitals(greenCubeMatDown);
        }

        ClusterMatrixCD_t greenMatsubaraUp = ioModelPtr_->FullCubeToIndep(greenCubeMatUp);
        ClusterMatrixCD_t greenMatsubaraDown = ioModelPtr_->FullCubeToIndep(greenCubeMatDown);

// Gather and stats of all the results for all cores
#ifndef AFM
        greenMatsubaraUp = 0.5 * (greenMatsubaraUp + greenMatsubaraDown);
        greenMatsubaraDown = greenMatsubaraUp;
#endif

        Result::ISResult isResult(obsScal, greenMatsubaraUp, greenMatsubaraDown, fillingAndDocc_.fillingUp(),
                                  fillingAndDocc_.fillingDown());
        std::vector<Result::ISResult> isResultVec;
#ifdef HAVEMPI

        mpi::communicator world;
        if (mpiUt::Tools::Rank() == mpiUt::Tools::master)
        {
            mpi::gather(world, isResult, isResultVec, mpiUt::Tools::master);
        }
        else
        {
            mpi::gather(world, isResult, mpiUt::Tools::master);
        }
        if (mpiUt::Tools::Rank() == mpiUt::Tools::master)
        {
            mpiUt::IOResult::SaveISResults(isResultVec, *ioModelPtr_, dataCT_->beta_);
        }

#else
        isResultVec.push_back(isResult);
        mpiUt::IOResult::SaveISResults(isResultVec, *ioModelPtr_, dataCT_->beta_);
#endif

        // Start: This should be in PostProcess.cpp ?
        // Start of observables that are easier and ok to do once all has been saved (for exemples, depends only on final green function)
        // Get KinecticEnergy
        // #ifndef DCA
        //                 if (mpiUt::Tools::Rank() == mpiUt::Tools::master)
        //                 {
        //                         std::ifstream fin("Obs.json");
        //                         Json results;
        //                         fin >> results;
        //                         fin.close();

        //                         std::cout << "Start Calculating Kinetic Energy " << std::endl;
        //                         KineticEnergy<TModel, TIOModel> kEnergy(modelPtr_, ioModelPtr_->ReadGreenDat("greenUp.dat", NOrb_));
        //                         results["KEnergy"] = {kEnergy.GetKineticEnergy(), 0.0};
        //                         std::cout << "End Calculating Kinetic Energy " << std::endl;

        //                         std::ofstream fout("Obs.json");
        //                         fout << std::setw(4) << results << std::endl;
        //                         fout.close();
        //                 }

        //                 //End: This should be in PostProcess.cpp ?
        // #endif
        // ioModelPtr_->SaveCube("greenUp.dat", modelPtr_->greenCluster0MatUp().data(), modelPtr_->beta());
        Logging::Info("End of Observables.Save()");
    }

  private:
    std::shared_ptr<ISDataCT> dataCT_;
    std::shared_ptr<Model_t> modelPtr_;
    std::shared_ptr<IOModel_t> ioModelPtr_;
    Utilities::EngineTypeFibonacci3217_t rng_;
    std::shared_ptr<Utilities::UniformRngFibonacci3217_t> urngPtr_;

    GreenBinning greenBinningUp_;
    GreenBinning greenBinningDown_;
    FillingAndDocc fillingAndDocc_;

    Matrix_t Maveraged_;

    //=======Measured quantities
    double signMeas_;
    double expOrder_;

    size_t NMeas_;

    const size_t NOrb_;
    const bool averageOrbitals_;
};

} // namespace Obs
} // namespace Markov
