#pragma once

#include "../Utilities/Logging.hpp"
#include "ABC_MonteCarlo.hpp"
#include <chrono>
#include <ctime>

namespace MC
{

struct Timer
{
    Timer() = default;
    void Start(double duration)
    {
        duration_ = duration;
        start_ = std::chrono::steady_clock::now();
    };

    static void PrintTime()
    {
        if (mpiUt::Tools::Rank() == mpiUt::Tools::master)
        {
            auto timeChrono = std::chrono::system_clock::now();
            std::time_t timeNow = std::chrono::system_clock::to_time_t(timeChrono);
            std::cout << "\t " << std::ctime(&timeNow) << std::endl;
        }

        return;
    }

    bool End() { return std::chrono::duration_cast<std::chrono::seconds>(std::chrono::steady_clock::now() - start_).count() > duration_; };

  private:
    double duration_;
    std::chrono::steady_clock::time_point start_;
};

template <typename TMarkovChain_t> class MonteCarlo : public ABC_MonteCarlo
{
  public:
    MonteCarlo(const std::shared_ptr<TMarkovChain_t> &markovchainPtr, const Json &jj)
        : markovchainPtr_(markovchainPtr), thermalizationTime_(jj["monteCarlo"]["thermalizationTime"].get<double>()),
          measurementTime_(jj["monteCarlo"]["measurementTime"].get<double>()), updatesMeas_(jj["solver"]["updatesMeas"].get<size_t>()),
          cleanUpdate_(jj["solver"]["cleanUpdate"].get<size_t>()), NMeas_(0), NCleanUpdates_(0),
          thermFromConfig_(jj["monteCarlo"]["thermFromConfig"].get<bool>())
    {
    }

    ~MonteCarlo(){};

    void RunMonteCarlo()
    {
        Timer timer;

        // if (thermFromConfig_)
        // {
        //     // markovchainPtr_->ThermalizeFromConfig();
        // }
        // else
        // {

        Logging::Info("Start Thermalization. ");

        timer.Start(60.0 * thermalizationTime_);
        while (true)
        {
            markovchainPtr_->DoStep();
            if (markovchainPtr_->updatesProposed() % updatesMeas_ == 0)
            {
                if (timer.End())
                {
                    break;
                }
                ++NMeas_;
            }

            if (markovchainPtr_->updatesProposed() % cleanUpdate_ == 0)
            {
                markovchainPtr_->CleanUpdate();
            }
        }

        markovchainPtr_->SaveTherm();
        Logging::Info("End Thermalization.: ");
        // }

        NMeas_ = 0;
        timer.Start(60.0 * measurementTime_);
        Logging::Info("Start Measurements. ");

        while (true)
        {
            markovchainPtr_->DoStep(); // One simple sweep

            if (markovchainPtr_->updatesProposed() % updatesMeas_ == 0)
            {
                if (timer.End())
                {
                    break;
                }
                markovchainPtr_->Measure();
                NMeas_++;
            }

            if (markovchainPtr_->updatesProposed() % cleanUpdate_ == 0)
            {
                markovchainPtr_->CleanUpdate();
                ++NCleanUpdates_;
            }
        }

        Logging::Debug("NCleanUpdates = " + std::to_string(NCleanUpdates_));
        Logging::Info("End Measurements.");
        markovchainPtr_->SaveMeas();
    }

    // Getters
    size_t NMeas() const { return NMeas_; }
    size_t NCleanUpdates() const { return NCleanUpdates_; }
    size_t updatesProposed() const { return markovchainPtr_->updatesProposed(); }

  private:
    // attributes
    const std::shared_ptr<TMarkovChain_t> markovchainPtr_;
    const double thermalizationTime_;
    const double measurementTime_;
    const size_t updatesMeas_;
    const size_t cleanUpdate_;

    size_t NMeas_;
    size_t NCleanUpdates_;
    bool thermFromConfig_;
}; // namespace MC
} // namespace MC
