#pragma once

#include "GreenMat.hpp"
#include "Fourier.hpp"
#include "Utilities.hpp"
#include "MPIUtilities.hpp"
#include "IO.hpp"

namespace GreenTau
{

using namespace GreenMat;
using Vector_t = std::vector<double>;
using Data_t = std::vector<Vector_t>;

class GreenCluster0Tau
{
    //definit par la fct hyb, tloc, mu et beta et un Nombre de slice de temps NTau

  public:
    const double EPS = 1e-13;
    const double deltaTau = 0.008;

    // GreenCluster0Tau() : gfMatCluster_(), beta_(), NTau_(){};

    GreenCluster0Tau(const GreenCluster0Mat &gfMatCluster, const Json &jj) : ioModel_(jj),
                                                                             gfMatCluster_(gfMatCluster),
                                                                             beta_(gfMatCluster.beta()),
                                                                             NTau_(std::max<double>(jj["NTAU"].get<double>(), beta_ / deltaTau)),
                                                                             NOrb_(gfMatCluster_.n_rows() / ioModel_.Nc)
    {
        mpiUt::Print("Creating gtau ");
        assert(NOrb_ >= 1);
        data_.resize(ioModel_.GetNIndepSuperSites(NOrb_));

#ifdef HAVEMPI
        mpi::communicator world;
        BuildParallel();
#else
        BuildSerial();
#endif

        if (mpiUt::Rank() == mpiUt::master)
        {
            Save("gtau.dat");
        }

        gfMatCluster_.clear();
        mpiUt::Print("gtau Created");
    };

    Vector_t BuildOneGTau(const size_t &indepSuperSiteIndex) //return g_i(tau)
    {
        Vector_t result(NTau_ + 1);
        const std::pair<size_t, size_t> indices = ioModel_.GetIndices(indepSuperSiteIndex, NOrb_);
        const size_t s1 = indices.first;
        const size_t s2 = indices.second;

        for (size_t tt = 0; tt < NTau_ + 1; tt++)
        {
            Tau_t tau = gfMatCluster_.beta() * (static_cast<double>(tt)) / static_cast<double>(NTau_);
            if (tt == 0)
            {
                tau += EPS;
            }
            if (tt == NTau_)
            {
                tau -= EPS;
            }

            const SiteVectorCD_t greenMat = gfMatCluster_.data().tube(s1, s2);
            const double fm = gfMatCluster_.fm()(s1, s2).real();
            const double sm = gfMatCluster_.sm()(s1, s2).real();
            const double tm = gfMatCluster_.tm()(s1, s2).real();
            result.at(tt) = Fourier::MatToTauAnalytic(greenMat, tau, beta_, fm, sm, tm);
        }

        return result;
    }

    void BuildSerial()
    {
        for (size_t i = 0; i < ioModel_.GetNIndepSuperSites(NOrb_); i++)
        {
            data_.at(i) = BuildOneGTau(i);
        }
    }

#ifdef HAVEMPI
    void BuildParallel()
    {
        mpi::communicator world;

        std::vector<Data_t> dataVec;
        size_t ii = 0;
        while (ii * mpiUt::NWorkers() < ioModel_.GetNIndepSuperSites(NOrb_))
        {
            const size_t indepSuperSiteIndex = mpiUt::Rank() + ii * mpiUt::NWorkers();
            Vector_t g0Tau;

            if (indepSuperSiteIndex < ioModel_.GetNIndepSuperSites(NOrb_))
            {
                g0Tau = BuildOneGTau(indepSuperSiteIndex);
            }

            Data_t dataResult;
            mpi::all_gather(world, g0Tau, dataResult);
            dataVec.push_back(dataResult);
            ii++;
        }

        //There will be empty vectors, important not to count them here.
        size_t jj = 0;
        for (size_t ll = 0; ll < dataVec.size() && jj < ioModel_.GetNIndepSuperSites(NOrb_); ll++)
        {

            for (size_t kk = 0; kk < dataVec.at(ll).size() && jj < ioModel_.GetNIndepSuperSites(NOrb_); kk++)
            {
                data_.at(jj) = dataVec.at(ll).at(kk);
                jj++;
            }
        }
    }
#endif

    ~GreenCluster0Tau() = default;

    GreenCluster0Mat gfMatCluster() const { return gfMatCluster_; };
    size_t NTau() const { return NTau_; };

    void Clear()
    {
        data_.clear();
        gfMatCluster_.clear();
    }

    double operator()(const SuperSite_t &s1, const SuperSite_t &s2, const Tau_t &tauIn)
    {
        assert(NOrb_ >= 1);
        double tau = tauIn - EPS;

        double aps = 1.0;

        if (tau < 0.0)
        {
            tau += beta_;
            aps = -1.0;
        }

        const size_t ll = ioModel_.FindIndepSuperSiteIndex(s1, s2, NOrb_);
        const double nt = std::abs(tau) / beta_ * static_cast<double>(NTau_);
        const size_t n0 = static_cast<size_t>(nt);
        const double greentau0 = aps * ((1.0 - (nt - n0)) * data_.at(ll).at(n0) + (nt - n0) * data_[ll].at(n0 + 1));
        return greentau0;
    }

    const GreenCluster0Tau &operator=(const GreenCluster0Tau &gf)
    {
        if (this == &gf)
            return *this; //évite les boucles infinies
        gfMatCluster_ = gf.gfMatCluster_;
        NTau_ = gf.NTau_;
        beta_ = gf.beta_;
        NOrb_ = gf.NOrb_;
        data_ = gf.data_;
        return *this;
    }

    void Save(std::string fileName)
    {

        std::ofstream fout(fileName);

        for (size_t tt = 0; tt < NTau_ + 1; tt++)
        {
            fout << beta_ * double(tt) / (static_cast<double>(NTau_)) << " ";
            for (size_t ii = 0; ii < data_.size(); ii++)
            {
                fout << data_.at(ii).at(tt) << " ";
            }
            fout << "\n";
        }
        fout.close();
    }

  private:
    IO::Base_IOModel ioModel_;
    GreenCluster0Mat gfMatCluster_;
    Data_t data_;

    double beta_;
    size_t NTau_;
    size_t NOrb_;
};
} // namespace GreenTau
