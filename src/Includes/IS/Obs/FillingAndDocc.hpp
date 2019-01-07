#pragma once

#include "../../Utilities/Utilities.hpp"
#include "../../Utilities/LinAlg.hpp"
#include "../ISData.hpp"

namespace Markov
{
namespace Obs
{

using IOModel_t = IO::Base_IOModel;
using Model_t = Models::ABC_Model_2D;

class FillingAndDocc
{

  public:
    FillingAndDocc(const std::shared_ptr<ISDataCT> &dataCT, const std::shared_ptr<IOModel_t> &ioModelPtr,
                   std::shared_ptr<Utilities::UniformRngFibonacci3217_t> urngPtr, const size_t &N_T_INV) : dataCT_(dataCT),
                                                                                                           ioModelPtr_(ioModelPtr),
                                                                                                           urngPtr_(urngPtr),
                                                                                                           N_T_INV_(N_T_INV)
    {
        const size_t LL = dataCT_->NOrb_ * ioModelPtr_->fillingSites().size();

        fillingUpCurrent_.resize(LL, 0.0);
        fillingDownCurrent_.resize(LL, 0.0);
        doccCurrent_.resize(LL, 0.0);
        SzCurrent_.resize(LL, 0.0);

        fillingUp_.resize(LL, 0.0);
        fillingDown_.resize(LL, 0.0);
        docc_.resize(LL, 0.0);
        Sz_.resize(LL, 0.0);
    }

    std::vector<double> fillingUp() const { return fillingUp_; }
    std::vector<double> fillingDown() const { return fillingDown_; }
    std::valarray<double> fillingDownCurrent() const { return fillingDownCurrent_; }
    std::valarray<double> fillingUpCurrent() const { return fillingUpCurrent_; }
    std::map<std::string, double> GetObs() const { return obsmap_; }

    double fillingUpTotalCurrent()
    {
        double fillingUpTotalCurrent = 0.0;

        size_t index = 0;
        for (size_t oo = 0; oo < dataCT_->NOrb_; oo++)
        {
            for (size_t kk = 0; kk < ioModelPtr_->fillingSites().size(); kk++)
            {

                const double factFilling = static_cast<double>(ioModelPtr_->nOfAssociatedSites().at(ioModelPtr_->fillingSitesIndex().at(kk)));
                fillingUpTotalCurrent += factFilling * fillingUpCurrent_[index];
                index++;
            }
        }

        return fillingUpTotalCurrent;
    }

    double fillingDownTotalCurrent()
    {
        double fillingDownTotalCurrent = 0.0;
        size_t index = 0;

        for (size_t oo = 0; oo < dataCT_->NOrb_; oo++)
        {
            for (size_t kk = 0; kk < ioModelPtr_->fillingSites().size(); kk++)
            {

                const double factFilling = static_cast<double>(ioModelPtr_->nOfAssociatedSites().at(ioModelPtr_->fillingSitesIndex().at(kk)));
                fillingDownTotalCurrent += factFilling * fillingDownCurrent_[index];
                index++;
            }
        }
        return fillingDownTotalCurrent;
    }

    double doccTotalCurrent()
    {
        double doccTotalCurrent = 0.0;
        size_t index = 0;
        for (size_t oo = 0; oo < dataCT_->NOrb_; oo++)
            for (size_t kk = 0; kk < ioModelPtr_->fillingSites().size(); kk++)
            {

                const double factFilling = static_cast<double>(ioModelPtr_->nOfAssociatedSites().at(ioModelPtr_->fillingSitesIndex().at(kk)));
                doccTotalCurrent += factFilling * doccCurrent_[index];
                index++;
            }
        return doccTotalCurrent;
    };

    void ResetCurrent()
    {
        doccCurrent_ = 0.0;
        SzCurrent_ = 0.0;
        fillingUpCurrent_ = 0.0;
        fillingDownCurrent_ = 0.0;
    }

    void MeasureFillingAndDocc()
    {
        // mpiUt::Print("Start of MeasureFillingAndDocc ");
        ResetCurrent();

        const size_t KK = dataCT_->vertices_.size();
        const size_t KKUp = dataCT_->vertices_.NUp();
        const size_t KKDown = dataCT_->vertices_.NDown();

        const double eps = 1e-12;

        // assert(KKUp == KKDown);
        assert(2 * KK == KKUp + KKDown);
        SiteVector_t vec1Up(KKUp);
        SiteVector_t vec2Up(KKUp);
        SiteVector_t vec1Down(KKDown);
        SiteVector_t vec2Down(KKDown);

        const double sign = static_cast<double>(dataCT_->sign_);
        const size_t fillingSize = ioModelPtr_->fillingSites().size();

        for (size_t oIndex = 0; oIndex < dataCT_->NOrb_; oIndex++)
        {
            for (Site_t ii = 0; ii < fillingSize; ii++)
            {
                const size_t index = oIndex * fillingSize + ii;

                const Site_t s1 = ioModelPtr_->fillingSites()[ii];
                const SuperSite_t superSite1{s1, oIndex};

                for (size_t nsamples = 0; nsamples < N_T_INV_; nsamples++)
                {
                    const double tauRng = (*urngPtr_)() * dataCT_->beta_;
                    const Site_t siteRng = ioModelPtr_->FindSitesRng(s1, s1, (*urngPtr_)()).first;
                    const SuperSite_t superSiteRng{siteRng, oIndex};

                    for (size_t iUp = 0; iUp < KKUp; iUp++)
                    {
                        const SuperSite_t superSite = dataCT_->vertices_.atUp(iUp).superSite();
                        const Tau_t tt = dataCT_->vertices_.atUp(iUp).tau();

                        vec1Up(iUp) = dataCT_->green0CachedUp_(superSiteRng, superSite, tauRng - tt);
                        vec2Up(iUp) = dataCT_->green0CachedUp_(superSite, superSiteRng, tt - tauRng);
                    }

                    for (size_t iDown = 0; iDown < KKDown; iDown++)
                    {
                        const SuperSite_t superSite = dataCT_->vertices_.atDown(iDown).superSite();
                        const Tau_t tt = dataCT_->vertices_.atDown(iDown).tau();
#ifdef AFM
                        vec1Down(iDown) = dataCT_->green0CachedDown_(superSiteRng, superSite, tauRng - tt);
                        vec2Down(iDown) = dataCT_->green0CachedDown_(superSite, superSiteRng, tt - tauRng);
#else
                        vec1Down(iDown) = dataCT_->green0CachedUp_(superSiteRng, superSite, tauRng - tt);
                        vec2Down(iDown) = dataCT_->green0CachedUp_(superSite, superSiteRng, tt - tauRng);
#endif
                    }

                    double dotup = 0.0;
                    double dotdown = 0.0;

                    if (KK)
                    {
                        dotup = LinAlg::Dot(vec1Up, *(dataCT_->MupPtr_), vec2Up);
                        dotdown = LinAlg::Dot(vec1Down, *(dataCT_->MdownPtr_), vec2Down);
                    }

                    const double green00Up = dataCT_->green0CachedUp_(superSite1, superSite1, -eps);

#ifdef AFM
                    const double green00Down = dataCT_->green0CachedDown_(superSite1, superSite1, -eps);
#else
                    const double green00Down = dataCT_->green0CachedUp_(superSite1, superSite1, -eps);
#endif
                    const double nUptmp = green00Up - dotup;
                    const double nDowntmp = green00Down - dotdown;

                    fillingUpCurrent_[index] += sign * nUptmp;
                    fillingDownCurrent_[index] += sign * nDowntmp;
                    doccCurrent_[index] += sign * (nUptmp * nDowntmp);

                    const double ndiff = nUptmp - nDowntmp;
                    SzCurrent_[index] += sign * ndiff;
                }

                fillingUpCurrent_[index] /= static_cast<double>(N_T_INV_);
                fillingDownCurrent_[index] /= static_cast<double>(N_T_INV_);
                doccCurrent_[index] /= static_cast<double>(N_T_INV_);
                SzCurrent_[index] /= static_cast<double>(N_T_INV_);

                fillingUp_.at(index) += fillingUpCurrent_[index];
                fillingDown_[index] += fillingDownCurrent_[index];
                docc_[index] += doccCurrent_[index];
                Sz_[index] += SzCurrent_[index];
            }
        }

        // mpiUt::Print("End of MeasureFillingAndDocc");
    }

    void Finalize(const double &signMeas, const double &NMeas)
    {

        const double fact = NMeas * signMeas;
        double fillingTotal = 0.0;
        double doccTotal = 0.0;
        double SzTotal = 0.0;
        std::string nName;

        for (size_t oIndex = 0; oIndex < dataCT_->NOrb_; oIndex++)
        {
            for (size_t kk = 0; kk < ioModelPtr_->fillingSites().size(); kk++)
            {
                const size_t index = oIndex * ioModelPtr_->fillingSites().size() + kk;

                fillingUp_.at(index) /= fact;
                fillingDown_.at(index) /= fact;

                nName = "n" + std::to_string(ioModelPtr_->fillingSites().at(kk)) + "_" + std::to_string(oIndex);
                obsmap_[nName + "Up"] = fillingUp_.at(index);
                obsmap_[nName + "Down"] = fillingDown_.at(index);

                docc_.at(index) /= fact;
                nName = "docc" + std::to_string(ioModelPtr_->fillingSites().at(kk)) + "_" + std::to_string(oIndex);
                obsmap_[nName] = docc_.at(index);

                Sz_.at(index) /= fact;
                nName = "Sz" + std::to_string(ioModelPtr_->fillingSites().at(kk)) + "_" + std::to_string(oIndex);
                obsmap_[nName] = Sz_.at(index);

                const double factFilling = static_cast<double>(ioModelPtr_->nOfAssociatedSites().at(ioModelPtr_->fillingSitesIndex().at(kk))) / static_cast<double>(ioModelPtr_->Nc);
                fillingTotal += factFilling * (fillingUp_.at(index) + fillingDown_.at(index));
                doccTotal += factFilling * docc_.at(index);
                SzTotal += factFilling * Sz_.at(index);
            }
        }

        obsmap_["n"] = fillingTotal;
        obsmap_["docc"] = doccTotal;
        obsmap_["Sz"] = SzTotal;
    }

  private:
    std::shared_ptr<ISDataCT> dataCT_;
    std::shared_ptr<IOModel_t> ioModelPtr_;
    std::shared_ptr<Utilities::UniformRngFibonacci3217_t> urngPtr_;

    std::map<std::string, double> obsmap_;

    //current values for the current config
    std::valarray<double> fillingUpCurrent_;
    std::valarray<double> fillingDownCurrent_;
    std::valarray<double> doccCurrent_;
    std::valarray<double> SzCurrent_;

    //values that are accumulated
    std::vector<double> fillingUp_;
    std::vector<double> fillingDown_;
    std::vector<double> docc_;
    std::vector<double> Sz_;

    const size_t N_T_INV_;
}; // class FillingAndDocc

} // namespace Obs
} // namespace Markov