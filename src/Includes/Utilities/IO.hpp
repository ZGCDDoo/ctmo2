#pragma once

#include "Utilities.hpp"
// #include "MPIUtilities.hpp"
#include "IOConstruct.hpp"

namespace IO
{

class Base_IOModel
{

  public:
    const size_t INVALID = 999;
    const size_t Nx;
    const size_t Ny;
    const size_t Nz;
    const size_t Nc;

    //get the number of indep orbitals
    static size_t GetNOrbIndep(const size_t &NOrb)
    {
        // size_t NOrbIndep = 0;
        // for (Orbital_t o1 = 0; o1 < NOrb; ++o1)
        // {
        //     for (Orbital_t o2 = o1; o2 < NOrb; ++o2)
        //     {
        //         NOrbIndep++;
        //     }
        // }

        //The number of indep orbitals should be equal to the number of elements in a triangular matrix, i.e:
        // assert(NOrbIndep == NOrb * (NOrb + 1) / 2);
        // return NOrbIndep;
        return (NOrb * (NOrb + 1) / 2);
    }

    Base_IOModel(const Json &jj) : Nx(jj["Nx"].get<size_t>()),
                                   Ny(jj["Ny"].get<size_t>()),
                                   Nz(jj["Nz"].get<size_t>()),
                                   Nc(Nx * Ny * Nz)
    {
        std::cout << "Start Base_IOModel construction " << std::endl;
        GreenSites_ = BuildGreenSites(jj["ModelFile"].get<std::string>());
        indepSites_ = BuildIndepSites(GreenSites_);
        FinishConstructor();
        std::cout << "End Base_IOModel construction " << std::endl;
    };

    void FinishConstructor()
    {
        ConstructfullSiteToIndepSite();
        SetnOfAssociatedSites();
        ConstructFillingSites();
        AssertSanity();
    }

    void ConstructFillingSites()
    {
        const size_t KK = indepSites_.size();
        for (size_t ii = 0; ii < KK; ii++)
        {
            Site_t s1 = indepSites_.at(ii).first;
            Site_t s2 = indepSites_.at(ii).second;
            if (s1 == s2)
            {
                fillingSites_.push_back(s1);
                fillingSitesIndex_.push_back(ii);
            }
        }
    }

    void ConstructfullSiteToIndepSite()
    {
        //construct equivalentSites_ also.
        equivalentSites_.clear();
        equivalentSites_.resize(indepSites_.size());

        for (size_t ii = 0; ii < indepSites_.size(); ii++)
        {
            std::pair<size_t, size_t> pairii = indepSites_.at(ii);
            equivalentSites_.at(ii).push_back(pairii);
        }

        for (Site_t s1 = 0; s1 < Nc; ++s1)
        {
            for (Site_t s2 = 0; s2 < Nc; ++s2)
            {

                const std::pair<Site_t, Site_t> pairSites = GreenSites_.at(s1).at(s2);
                std::vector<std::pair<size_t, size_t>>::iterator llit = std::find(indepSites_.begin(), indepSites_.end(), pairSites);
                if (llit == indepSites_.end())
                {
                    //Try the pair with first and second exchanged:
                    llit = std::find(indepSites_.begin(), indepSites_.end(), std::make_pair(pairSites.second, pairSites.first));
                    if (llit == indepSites_.end())
                    {
                        throw std::runtime_error("Bad index in FindIndepSiteIndex!");
                    }
                }
                const size_t llDistance = std::distance(indepSites_.begin(), llit);
                fullSiteToIndepSite_.push_back(llDistance);
                equivalentSites_.at(llDistance).push_back(std::make_pair(s1, s2));

                if (s1 == s2)
                {
                    assert(indepSites_.at(llDistance).first == indepSites_.at(llDistance).second);
                }
            }
        }
    }

    size_t FindIndepSiteIndex(const Site_t &s1, const Site_t &s2) const
    {
        return fullSiteToIndepSite_.at(s1 * Nc + s2);
    }

    size_t FindIndepSuperSiteIndex(const SuperSite_t &s1, const SuperSite_t &s2, const size_t &NOrb) const
    {

        const size_t NSitesIndep = indepSites_.size();
        const size_t siteIndex = fullSiteToIndepSite_.at(s1.first * Nc + s2.first); //arranged by row major ordering here, one of the only places where this happens
        const size_t orbitalIndex = Utilities::GetIndepOrbitalIndex(s1.second, s2.second, NOrb);
        return (siteIndex + orbitalIndex * NSitesIndep);
    }

    void SetnOfAssociatedSites()
    {
        nOfAssociatedSites_.resize(indepSites_.size());
        // nOfFillingSites_ = 0;

        for (size_t ii = 0; ii < Nc; ++ii)
        {
            for (size_t jj = 0; jj < Nc; ++jj)
            {
                const size_t ll = FindIndepSiteIndex(ii, jj);
                nOfAssociatedSites_.at(ll) = nOfAssociatedSites_.at(ll) + 1;
            }
        }
    }

#ifdef DCA
    //read a green in .dat format.
    ClusterCubeCD_t ReadGreenKDat(const std::string &filename, const size_t &NOrb) const
    {
        // mpiUt::Print("In IOModel ReadGreenKDat ");

        const size_t shutUpWarning = NOrb;
        std::cout << shutUpWarning << "WARNING, Norb not implemented in ReadGreenKDat" << std::endl;

        ClusterMatrix_t fileMat;
        ClusterMatrixCD_t tmp(Nc, Nc);
        fileMat.load(filename);
        assert(!fileMat.has_nan());
        assert(!fileMat.has_inf());
        assert(fileMat.n_cols == 2 * Nc + 1);
        fileMat.shed_col(0); // we dont want the matsubara frequencies

        ClusterCubeCD_t cubetmp(Nc, Nc, fileMat.n_rows);
        cubetmp.zeros();

        for (size_t n = 0; n < cubetmp.n_slices; ++n)
        {
            tmp.zeros();
            for (size_t KIndex = 0; KIndex < Nc; ++KIndex)
            {
                tmp(KIndex, KIndex) = cd_t(fileMat(n, 2 * KIndex), fileMat(n, 2 * KIndex + 1));
            }

            cubetmp.slice(n) = tmp;
        }

        return cubetmp;
    }
#endif

    //Read green in .dat format.
    ClusterCubeCD_t ReadGreenDat(const std::string &filename, const size_t &NOrb) const
    {
        // mpiUt::Print("In IOModel ReadGreenNDat ");

        const size_t NN = Nc * NOrb;
        const size_t NOrbIndep = GetNOrbIndep(NOrb);
        const size_t NSitesIndep = indepSites_.size();

        ClusterMatrix_t fileMat;
        ClusterMatrixCD_t tmp(NN, NN);
        fileMat.load(filename);
        assert(!fileMat.has_nan());
        assert(!fileMat.has_inf());

        assert(fileMat.n_cols == NOrbIndep * 2 * NSitesIndep + 1);
        fileMat.shed_col(0); // we dont want the matsubara frequencies

        ClusterCubeCD_t cubetmp(NN, NN, fileMat.n_rows);

        for (size_t n = 0; n < cubetmp.n_slices; ++n)
        {
            for (Orbital_t o1 = 0; o1 < NOrb; ++o1)
            {
                for (Orbital_t o2 = o1; o2 < NOrb; ++o2)
                {
                    for (size_t ii = 0; ii < Nc; ++ii)
                    {
                        for (size_t jj = 0; jj < Nc; ++jj)
                        {
                            const size_t indexIndepSuperSite = FindIndepSuperSiteIndex(std::make_pair(ii, o1), std::make_pair(jj, o2), NOrb);
                            tmp(ii + o1 * Nc, jj + o2 * Nc) = cd_t(fileMat(n, 2 * indexIndepSuperSite), fileMat(n, 2 * indexIndepSuperSite + 1));
                            tmp(jj + o2 * Nc, ii + o1 * Nc) = tmp(ii + o1 * Nc, jj + o2 * Nc); //symmetrize
                        }
                    }
                }
            }

            cubetmp.slice(n) = tmp;
        }

        return cubetmp;
    }

    void SaveCube(const std::string &fname, const ClusterCubeCD_t &green, const double &beta,
                  const size_t &NOrb, const size_t &precision = 14, const bool &saveArma = false) const
    {
        assert(!green.has_nan());
        assert(!green.has_inf());
        const size_t NMat = green.n_slices;
        const size_t NOrbIndep = GetNOrbIndep(NOrb);
        const size_t NSitesIndep = this->indepSites_.size();

        assert(green.n_rows == green.n_cols);
        assert(green.n_rows == Nc * NOrb);
        ClusterMatrixCD_t greenOut(NMat, NOrbIndep * NSitesIndep);

        std::ofstream fout;
        fout.open(fname + std::string(".dat"), std::ios::out);
        for (size_t nn = 0; nn < green.n_slices; nn++)
        {
            const double iwn = (2.0 * nn + 1.0) * M_PI / beta;
            fout << std::setprecision(precision) << iwn << " ";

            for (Orbital_t o1 = 0; o1 < NOrb; ++o1)
            {
                for (Orbital_t o2 = o1; o2 < NOrb; ++o2)
                {
                    for (Site_t ii = 0; ii < NSitesIndep; ++ii)
                    {
                        const Site_t r1 = this->indepSites_.at(ii).first;
                        const Site_t r2 = this->indepSites_.at(ii).second;

                        const cd_t value = green(r1 + o1 * Nc, r2 + o2 * Nc, nn);
                        const size_t indexIndepSuperSite = FindIndepSuperSiteIndex(std::make_pair(r1, o1), std::make_pair(r2, o2), NOrb);
                        greenOut(nn, indexIndepSuperSite) = value;
                        fout << std::setprecision(precision) << value.real()
                             << " "
                             << std::setprecision(precision) << value.imag()
                             << " ";
                    }
                }
            }
            fout << "\n";
        }

        fout.close();

        if (saveArma)
        {
            greenOut.save(fname + std::string(".arma"), arma::arma_ascii);
        }
    }

#ifdef DCA
    void SaveK(const std::string &fname, const ClusterCubeCD_t &green, const double &beta, const size_t &NOrb, const size_t &precision = 14) const
    {
        const size_t shutUpWarning = NOrb;
        std::cout << shutUpWarning << "WARNING, Norb not implemented in SaveK" << std::endl;

        std::ofstream fout;
        fout.open(fname + std::string(".dat"), std::ios::out);
        for (size_t nn = 0; nn < green.n_slices; ++nn)
        {
            const double iwn = (2.0 * nn + 1.0) * M_PI / beta;
            fout << std::setprecision(precision) << iwn << " ";

            for (Site_t ii = 0; ii < Nc; ++ii)
            {

                fout << std::setprecision(precision) << green(ii, ii, nn).real()
                     << " "
                     << std::setprecision(precision) << green(ii, ii, nn).imag()
                     << " ";
            }
            fout << "\n";
        }

        fout.close();
    }

#endif

    std::pair<size_t, size_t> GetIndices(const size_t &indepSuperSiteIndex, const size_t &NOrb) const
    {
        const size_t LL = indepSites_.size();

        const size_t indepSiteIndex = indepSuperSiteIndex % LL;
        const size_t indepOrbitalIndex = indepSuperSiteIndex / LL;

        const auto orbitalPair = GetIndicesOrbital(indepOrbitalIndex, NOrb);
        const size_t o1 = orbitalPair.first;
        const size_t o2 = orbitalPair.second;

        const std::pair<Site_t, Site_t> sites = indepSites_.at(indepSiteIndex);
        return {sites.first + o1 * Nc, sites.second + o2 * Nc};
    }

    std::pair<size_t, size_t> GetIndicesOrbital(const size_t &indepOrbitalIndex, const size_t &NOrb) const
    {
        size_t tmp = 0;
        for (size_t o1 = 0; o1 < NOrb; ++o1)
        {
            for (size_t o2 = o1; o2 < NOrb; ++o2)
            {
                if (tmp == indepOrbitalIndex)
                {
                    return {o1, o2};
                }
                ++tmp;
            }
        }

        throw std::runtime_error("Shit man");
    }

    size_t GetNIndepSuperSites(const size_t &NOrb) const
    {
        return (indepSites_.size() * GetNOrbIndep(NOrb));
    }

    template <typename T1_t, typename T2_t = ClusterMatrixCD_t>
    T2_t IndepToFull(const T1_t &indepElements, const size_t &NOrb) const //in practice T1_t will be a Sitevector_t or SitevectorCD_t
    {
        assert(indepElements.n_elem == GetNOrbIndep(NOrb) * indepSites_.size());
        T2_t fullMatrix(NOrb * Nc, NOrb * Nc);
        fullMatrix.zeros();

        for (Orbital_t o1 = 0; o1 < NOrb; ++o1)
        {
            for (Orbital_t o2 = o1; o2 < NOrb; ++o2)
            {
                for (size_t ii = 0; ii < Nc; ++ii)
                {
                    for (size_t jj = 0; jj < Nc; ++jj)
                    {

                        const size_t indexIndepSuperSite = FindIndepSuperSiteIndex(std::make_pair(ii, o1), std::make_pair(jj, o2), NOrb);
                        fullMatrix(ii + o1 * Nc, jj + o2 * Nc) = indepElements(indexIndepSuperSite);
                        fullMatrix(jj + o2 * Nc, ii + o1 * Nc) = fullMatrix(ii + o1 * Nc, jj + o2 * Nc); //symmetrize
                    }
                }
            }
        }

        return fullMatrix;
    }

    //from th full cube return the independant in tabular form
    template <typename T1_t>
    ClusterMatrixCD_t FullCubeToIndep(const T1_t &greenCube) const //T1_t = {ClusterCube_t or ClustercubeCD_t}
    {
        const size_t NOrb = greenCube.n_rows / Nc;
        const size_t NIndepSuperSites = GetNIndepSuperSites(NOrb);
        ClusterMatrixCD_t indepTabular(greenCube.n_slices, NIndepSuperSites);
        indepTabular.zeros();

        for (size_t i = 0; i < NIndepSuperSites; ++i)
        {
            const auto indices = GetIndices(i, NOrb);
            const size_t s1 = indices.first;
            const size_t s2 = indices.second;

            for (size_t n = 0; n < greenCube.n_slices; ++n)
            {
                indepTabular(n, i) = greenCube(s1, s2, n);
            }
        }

        return indepTabular;
    }

    void AssertSanity() const
    {
        size_t sum = 0.0;
        assert(fillingSites_.size() == fillingSitesIndex_.size());

        for (size_t ii = 0; ii < fillingSitesIndex_.size(); ++ii)
        {
            sum += nOfAssociatedSites_.at(fillingSitesIndex_.at(ii));
        }
        assert(sum == Nc);

        //make sure fullSitetoIndepSite is ok
        assert(fullSiteToIndepSite_.size() == Nc * Nc);

        //make sure equivalentSites_ is ok
        assert(equivalentSites_.size() == indepSites_.size());
        for (size_t ii = 0; ii < equivalentSites_.size(); ++ii)
        {
            std::pair<size_t, size_t> pairii = equivalentSites_.at(ii).at(0);
            for (size_t jj = 0; jj < equivalentSites_.at(ii).size(); ++jj)
            {
                const size_t s1 = equivalentSites_.at(ii).at(jj).first;
                const size_t s2 = equivalentSites_.at(ii).at(jj).second;
                assert(pairii == GreenSites_.at(s1).at(s2));
            }
        }
    }

    ClusterCubeCD_t AverageOrbitals(const ClusterCubeCD_t green) const
    {
        //For now only averages the diagonal blocks.
        const size_t n_rows = green.n_rows;
        const size_t n_cols = green.n_cols;
        const size_t n_slices = green.n_slices;
        assert(n_rows == n_cols);
        const size_t NOrb = green.n_rows / Nc;

        ClusterCubeCD_t result(n_rows, n_cols, n_slices);
        result.zeros();

        using arma::span;
        const size_t Ncm1 = Nc - 1;

        for (size_t oo = 0; oo < NOrb; ++oo)
        {
            result.subcube(span(0, Ncm1), span(0, Ncm1), span(0, n_slices - 1)) += green.subcube(span(Nc * oo, Ncm1 + Nc * oo), span(Nc * oo, Ncm1 + Nc * oo), span(0, n_slices - 1));
        }
        result /= static_cast<double>(NOrb);

        for (size_t oo = 1; oo < NOrb; ++oo)
        {
            result.subcube(span(Nc * oo, Ncm1 + Nc * oo), span(Nc * oo, Ncm1 + Nc * oo), span(0, n_slices - 1)) = result.subcube(span(0, Ncm1), span(0, Ncm1), span(0, n_slices - 1));
        }

        return result;
    }

    std::pair<size_t, size_t> FindSitesRng(const size_t &s1, const size_t &s2, const double &rngDouble) const
    {
        const size_t indepSiteIndex = FindIndepSiteIndex(s1, s2);
        const size_t equivalentSize = equivalentSites_.at(indepSiteIndex).size();
        const size_t intRng = rngDouble * equivalentSize;
        return equivalentSites_.at(indepSiteIndex).at(intRng);
    }

    //Getters
    std::vector<std::pair<size_t, size_t>> const indepSites() const
    {
        return indepSites_;
    };
    std::vector<std::vector<std::pair<size_t, size_t>>> const GreenSites() const
    {
        return GreenSites_;
    };
    std::vector<std::vector<std::pair<size_t, size_t>>> const equivalentSites() const
    {
        return equivalentSites_;
    };
    std::vector<size_t> const nOfAssociatedSites() const
    {
        return nOfAssociatedSites_;
    };
    std::vector<size_t> const fillingSites() const
    {
        return fillingSites_;
    };
    std::vector<size_t> const fillingSitesIndex() const
    {
        return fillingSitesIndex_;
    };
    std::vector<size_t> const downEquivalentSites() const
    {
        return downEquivalentSites_;
    };

  protected:
    std::vector<std::pair<size_t, size_t>> indepSites_;
    std::vector<std::vector<std::pair<size_t, size_t>>> GreenSites_;
    std::vector<std::vector<std::pair<size_t, size_t>>> equivalentSites_; // for example for square lattice equivalentsites.at(0) = {{0.0}, {1,1} , {2,2}, {3,3}}
    std::vector<size_t> nOfAssociatedSites_;
    std::vector<size_t> fullSiteToIndepSite_;
    std::vector<size_t> fillingSites_;
    std::vector<size_t> fillingSitesIndex_; //the indexes of the fillingsites in the indepSites_
    std::vector<size_t> downEquivalentSites_;
};

} // namespace IO
