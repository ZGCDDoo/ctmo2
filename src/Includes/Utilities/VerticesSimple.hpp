#pragma once

#include "Utilities.hpp"
#include "../Utilities/LinAlg.hpp"

#include "../Models/UTensorSimple.hpp"
#include "../Utilities/Logging.hpp"
#include <boost/math/special_functions/binomial.hpp>

namespace Diagrammatic
{
typedef LinAlg::Matrix_t Matrix_t;

enum class VertexType
{
    HubbardIntra,     //Hubbard intraorbital
    HubbardInter,     // Hubbard interorbital, different spins (U')
    HubbardInterSpin, // Hubbard interorbital same spin (U'-J_H)
    Phonon,
    Invalid
};

const size_t N_VERTEX_TYPES = 3; //for now, we dont do Phonon
const size_t INVALID = 999;

class VertexPart
{
  public:
    VertexPart(const VertexType vtype, const Tau_t &tau, const Site_t &site, const FermionSpin_t &spin,
               const size_t &orbital, const AuxSpin_t &aux) : vtype_(vtype),
                                                              tau_(tau),
                                                              site_(site),
                                                              spin_(spin),
                                                              orbital_(orbital),
                                                              superSite_{site, orbital},
                                                              aux_(aux) {}

    VertexPart &operator=(const VertexPart &vpart) = default;

    //Getters
    VertexType vtype() const { return vtype_; };
    Tau_t tau() const { return tau_; };
    Site_t site() const { return site_; };
    FermionSpin_t spin() const { return spin_; };
    Orbital_t orbital() const { return orbital_; };
    SuperSite_t superSite() const { return superSite_; };
    AuxSpin_t aux() const { return aux_; };

    bool operator==(const VertexPart &x) const
    {
        return ((x.tau() == tau_) && (x.site() == site_) && (x.spin() == spin_) && (x.orbital() == orbital_) && (x.superSite() == superSite_) && (x.aux() == aux_));
    }

  private:
    VertexType vtype_;
    Tau_t tau_;
    Site_t site_;
    FermionSpin_t spin_;
    Orbital_t orbital_;
    SuperSite_t superSite_;
    AuxSpin_t aux_;
};

class Vertex
{

  public:
    Vertex(const VertexType &vtype, const VertexPart &vStart, const VertexPart &vEnd,
           const double &probProb) : vtype_(vtype),
                                     vStart_(vStart),
                                     vEnd_(vEnd),
                                     probProb_(probProb)

    {
    }

    Vertex &operator=(const Vertex &vertex) = default;

    // Getters
    VertexType vtype() const { return vtype_; };
    double probProb() const { return probProb_; };
    VertexPart vStart() const { return vStart_; };
    VertexPart vEnd() const { return vEnd_; };

    //Setters

  private:
    VertexType vtype_;
    VertexPart vStart_;
    VertexPart vEnd_;
    double probProb_;
};

class Vertices
{

  public:
    Vertices() : key_(0){};

    void AppendVertex(const Vertex &vertex)
    {
        AssertSizes();

        data_.push_back(vertex);
        verticesKeysVec_.push_back(key_);
        AppendVertexPart(vertex.vStart());

        //VertexParts differ by one for their id if spins are the same
        if (vertex.vStart().spin() == vertex.vEnd().spin())
        {
            key_++;
            // assert(false);
        }

        AppendVertexPart(vertex.vEnd());

        //Update the id number once all the vertices parts have been inserted
        key_ += 3;
        AssertSizes();
    }

    void AssertSizes() const
    {

        assert(2 * data_.size() == (vPartUpVec_.size() + vPartDownVec_.size()));
        assert(indexPartUpVec_.size() == vPartUpVec_.size());
        assert(indexPartDownVec_.size() == vPartDownVec_.size());
        assert(data_.size() == verticesKeysVec_.size());
    }

    void Print() const
    {

        std::cout << "Start Print " << std::endl;

        for (size_t ii = 0; ii < indexPartUpVec_.size(); ii++)
        {
            std::cout << "indexPartUpVec_.at(ii)  = " << indexPartUpVec_.at(ii) << std::endl;
        }

        for (size_t ii = 0; ii < indexPartDownVec_.size(); ii++)
        {
            std::cout << "indexPartDownVec__.at(ii)  = " << indexPartDownVec_.at(ii) << std::endl;
        }
        std::cout << "End Print " << std::endl;
    }

    void AppendVertexPart(const VertexPart &vPart)
    {

        if (vPart.spin() == FermionSpin_t::Up)
        {
            vPartUpVec_.push_back(vPart);
            indexPartUpVec_.push_back(key_);
        }
        else
        {
            vPartDownVec_.push_back(vPart);
            indexPartDownVec_.push_back(key_);
        }
    }

    size_t GetKeyIndex(const UInt64_t &key, const FermionSpin_t &spin) const
    {
        std::vector<size_t> indices;
        if (spin == FermionSpin_t::Up)
        {
            //Find in order the vertexParts corresponding to the same vertex
            for (size_t ii = 0; ii < indexPartUpVec_.size(); ii++)
            {
                if (key == indexPartUpVec_.at(ii))
                {
                    indices.push_back(ii);
                }
            }
        }
        else
        {
            for (size_t ii = 0; ii < indexPartDownVec_.size(); ii++)
            {
                if (key == indexPartDownVec_.at(ii))
                {
                    indices.push_back(ii);
                }
            }
        }

        assert(indices.size() == 1);
        return (indices.at(0));
    }

    UInt64_t GetKey(const size_t &pp) const
    {
        AssertSizes();

        return verticesKeysVec_.at(pp);
    }

    std::vector<size_t> GetIndicesSpins(const size_t &pp, const FermionSpin_t &spin) const
    {
        const UInt64_t vertexKey = verticesKeysVec_.at(pp);
        const VertexPart x = data_.at(pp).vStart();
        const VertexPart y = data_.at(pp).vEnd();
        std::vector<size_t> indices;

        indices.push_back(GetKeyIndex(vertexKey, spin));

        if (x.spin() == y.spin())
        {
            indices.push_back(GetKeyIndex(vertexKey + 1, spin));
        }

        return indices;
    }

    void SaveConfig(const std::string &fname)
    {
        std::ofstream fout(fname);
        for (size_t ii = 0; ii < size(); ii++)
        {
            const auto x = data_.at(ii).vStart();
            const auto y = data_.at(ii).vEnd();
            fout << static_cast<int>(x.aux()) << " " << x.site() << " " << x.tau() << " " << static_cast<int>(x.spin()) << " " << static_cast<int>(y.spin()) << " " << x.orbital() << " " << y.orbital() << " " << std::endl;
        }
    }

    void RemoveVertex(const size_t &pp)
    {
        const size_t kkm1 = size() - 1;
        std::iter_swap(data_.begin() + pp, data_.begin() + kkm1);                       //swap the last vertex and the vertex pp in vertices.
        std::iter_swap(verticesKeysVec_.begin() + pp, verticesKeysVec_.begin() + kkm1); //swap the last vertex and the vertex pp in vertices.
        data_.pop_back();
        verticesKeysVec_.pop_back();
    }

    void PopBackVertexPart(const FermionSpin_t &spin)
    {
        if (spin == FermionSpin_t::Up)
        {
            vPartUpVec_.pop_back();
            indexPartUpVec_.pop_back();
        }
        else
        {
            vPartDownVec_.pop_back();
            indexPartDownVec_.pop_back();
        }
    }

    void SwapVertexPart(const size_t &pp1, const size_t &pp2, const FermionSpin_t &spin)
    {

        if (spin == FermionSpin_t::Up)
        {
            std::iter_swap(vPartUpVec_.begin() + pp1, vPartUpVec_.begin() + pp2);
            std::iter_swap(indexPartUpVec_.begin() + pp1, indexPartUpVec_.begin() + pp2);
        }
        else
        {

            std::iter_swap(vPartDownVec_.begin() + pp1, vPartDownVec_.begin() + pp2);
            std::iter_swap(indexPartDownVec_.begin() + pp1, indexPartDownVec_.begin() + pp2);
        }
    }

    //Getters
    size_t size() const
    {
        return data_.size();
    };
    size_t NUp() const { return vPartUpVec_.size(); };
    size_t NDown() const { return vPartDownVec_.size(); };
    std::vector<UInt64_t> verticesKeysVec() const { return verticesKeysVec_; }; // Each vertex has a unqique key identifying it

    Vertex at(const size_t &i) const { return data_.at(i); };

    VertexPart atUp(const size_t &i) { return vPartUpVec_.at(i); };
    VertexPart atDown(const size_t &i) { return vPartDownVec_.at(i); };

    void Clear()
    {
        data_.clear();
        vPartUpVec_.clear();
        vPartDownVec_.clear();
        indexPartUpVec_.clear();
        indexPartDownVec_.clear();
        verticesKeysVec_.clear();
        key_ = 0;
    }

  private:
    std::vector<Vertex> data_;
    std::vector<VertexPart> vPartUpVec_;
    std::vector<VertexPart> vPartDownVec_;
    std::vector<size_t> indexPartUpVec_;
    std::vector<size_t> indexPartDownVec_;  // Ex: The row and col 0 of Ndown_ is associtated to the vertexPart of the vertex given by the id of indexPartDownVec_.at(0)
    std::vector<UInt64_t> verticesKeysVec_; // Each vertex has a unqique key identifying it
    UInt64_t key_;

}; // namespace Diagrammatic

class AuxHelper
{
  public:
    explicit AuxHelper(const double &delta) : delta_(delta){}; //in futur, alpha tensor in constructor

    double auxValue(const FermionSpin_t &spin, const AuxSpin_t &aux) const
    {
        if (spin == FermionSpin_t::Up)
        {
            return ((aux == AuxSpin_t::Up) ? 1.0 + delta_ : -delta_);
        }
        else
        {
            return ((aux == AuxSpin_t::Down) ? 1.0 + delta_ : -delta_);
        }
    }

    double auxValueBar(const FermionSpin_t &spin, const AuxSpin_t &aux) const
    {
        if (spin == FermionSpin_t::Up)
        {
            return ((aux == AuxSpin_t::Up) ? -delta_ : 1.0 + delta_);
        }
        else
        {
            return ((aux == AuxSpin_t::Down) ? -delta_ : 1.0 + delta_);
        }
    }

    double auxPh(const AuxSpin_t &aux) const
    {
        return ((aux == AuxSpin_t::Up) ? 1.0 + delta_ : -delta_);
    }

    double FAux(const VertexPart &vp) const
    {

        if (vp.vtype() == Diagrammatic::VertexType::Phonon)
        {
            return (auxPh(vp.aux()) / (auxPh(vp.aux()) - 1.0));
        }
        else
        {
            return (auxValue(vp.spin(), vp.aux()) / (auxValue(vp.spin(), vp.aux()) - 1.0));
        }
    }

    double FAuxBar(const VertexPart &vp) const
    {
        //return FAux_sigma(-s);
        if (vp.vtype() == VertexType::Phonon)
        {
            return FAux(vp);
        }
        else
        {
            const AuxSpin_t sBar = (vp.aux() == AuxSpin_t::Up) ? AuxSpin_t::Down : AuxSpin_t::Up;
            return (auxValue(vp.spin(), sBar) / (auxValue(vp.spin(), sBar) - 1.0));
        }
    }

    double delta() const { return delta_; };

  private:
    const double delta_;
};

class VertexBuilder
{
  public:
    //must hold the alphas, the values of the U, U' and (U-J_H)
    VertexBuilder(const Json &jj, const size_t &Nc) : Utensor(jj),
                                                      auxHelper_(jj["model"]["delta"].get<double>()),
                                                      delta_(jj["model"]["delta"].get<double>()),
                                                      beta_(jj["model"]["beta"].get<double>()),
                                                      Nc_(Nc),
                                                      NOrb_(jj["model"]["nOrb"].get<size_t>()),
                                                      probU_(0.5),
                                                      factXi_(
                                                          1.0 / probU_ * (NOrb_ * NOrb_ * 2 * 2 / 2 - NOrb_) // Pauli principale and dont double count pairs
                                                          ),
                                                      isOrbitalDiagonal_(jj["solver"]["isOrbitalDiagonal"].get<bool>())

    {
    }

    Vertex BuildVertex(Utilities::UniformRngMt19937_t &urng)
    {
        const Tau_t tau = urng() * beta_;
        const Site_t site = urng() * Nc_;
        const AuxSpin_t aux = urng() < 0.5 ? AuxSpin_t::Up : AuxSpin_t::Down;

        Orbital_t o1 = urng() * NOrb_;
        Orbital_t o2 = urng() * NOrb_;
        FermionSpin_t spin1 = (urng() < 0.5) ? FermionSpin_t::Up : FermionSpin_t::Down;
        FermionSpin_t spin2 = (urng() < 0.5) ? FermionSpin_t::Up : FermionSpin_t::Down;

        VertexType vertextype = VertexType::Invalid;

        if (urng() < probU_) //Then build Electron-Eletron vertex
        {
            while ((o1 == o2) && (spin1 == spin2))
            {
                o1 = urng() * NOrb_;
                o2 = urng() * NOrb_;
                spin1 = (urng() < 0.5) ? FermionSpin_t::Up : FermionSpin_t::Down;
                spin2 = (urng() < 0.5) ? FermionSpin_t::Up : FermionSpin_t::Down;
            }
            if ((o1 == o2) && (spin1 != spin2))
            {
                vertextype = VertexType::HubbardIntra;
                const VertexPart vStart(vertextype, tau, site, FermionSpin_t::Up, o1, aux);
                const VertexPart vEnd(vertextype, tau, site, FermionSpin_t::Down, o2, aux);

                return Vertex(vertextype, vStart, vEnd, GetProbProb(vStart, vEnd));
            }
            else if ((o1 != o2) && (spin1 != spin2))
            {
                vertextype = VertexType::HubbardInter;
                const VertexPart vStart(vertextype, tau, site, FermionSpin_t::Up, o1, aux);
                const VertexPart vEnd(vertextype, tau, site, FermionSpin_t::Down, o2, aux);
                return Vertex(vertextype, vStart, vEnd, GetProbProb(vStart, vEnd));
            }
            else if ((o1 != o2) && (spin1 == spin2))
            {
                vertextype = VertexType::HubbardInterSpin;
                const VertexPart vStart(vertextype, tau, site, spin1, o1, aux);
                const VertexPart vEnd(vertextype, tau, site, spin2, o2, aux);
                return Vertex(vertextype, vStart, vEnd, GetProbProb(vStart, vEnd));
            }
            else
            {
                throw std::runtime_error("Miseria, Error in Vertices. Stupido!");
            }
        }
        else //Then build Electron-Phonon vertex
        {
            vertextype = VertexType::Phonon;
            const double deltaTau = GetDeltaTauPhonon(urng());
            double tau1 = urng() * beta_;
            double tauPrime = tau1 - deltaTau;

            if (tauPrime < 0.0)
            {
                tauPrime += beta_;
            }
            else if (tauPrime > beta_)
            {
                tauPrime -= beta_;
            }

//Use completely random insertion for GREEN_STYLE, testing purpose: ctmo and ctmo_green should give the same results.
#ifdef GREEN_STYLE
            tau1 = urng() * beta_;
            tauPrime = urng() * beta_;
#endif
            const VertexPart vStart(vertextype, tau1, site, spin1, o1, aux);
            const VertexPart vEnd(vertextype, tauPrime, site, spin2, o2, aux);

            return Vertex(vertextype, vStart, vEnd, GetProbProb(vStart, vEnd));
        }

        throw std::runtime_error("Miseria, Error in Vertices. Stupido!");
    }

    double GetProbProb(const VertexPart &x, const VertexPart &y)
    {
        assert(x.aux() == y.aux());
        assert(x.vtype() == y.vtype());

        const VertexType vtype = x.vtype();

        double U_xio1o2 = INVALID;

        if (vtype == VertexType::HubbardIntra)
        {
            U_xio1o2 = Utensor.U() / 2.0;
        }
        else if (vtype == VertexType::HubbardInter)
        {
            U_xio1o2 = isOrbitalDiagonal_ ? 0.0 : Utensor.UPrime() / 2.0;
        }
        else if (vtype == VertexType::HubbardInterSpin)
        {
            if (isOrbitalDiagonal_)
            {
                U_xio1o2 = 0.0;
            }
            else if (std::abs(Utensor.JH()) < 1e-10)
            {
                U_xio1o2 = 0.0;
            }
            else
            {
                U_xio1o2 = (Utensor.UPrime() - Utensor.JH()) / 2.0;
            }
        }
        else if (vtype == VertexType::Phonon)
        {
            const double w0 = Utensor.w0Phonon();
            if ((std::abs(w0) < 1e-10) || (std::abs(Utensor.gPhonon()) < 1e-10))
            {
                return 0.0;
            }

            //Big M is = 1.0
            U_xio1o2 = Utensor.gPhonon() * Utensor.gPhonon() / (4.0 * w0 * w0);
            const double factPh = 1.0 / (1.0 - probU_) * NOrb_ * NOrb_ * 2.0 * 2.0;

#ifdef GREEN_STYLE
            const double gtauPH = PhononPropagator(x.tau() - y.tau());
            return (-static_cast<double>(Nc_) * beta_ * beta_ * factPh * U_xio1o2 * 2.0 * gtauPH); //For testing purposes, green style is defined by sampling the two times uniformaly
#else
            const double fact = 1.0 / (auxHelper_.FAux(x) - 1.0);
            return (static_cast<double>(Nc_) * beta_ * factPh * U_xio1o2 * fact * fact * 2.0);
#endif
        }
        else
        {
            throw std::runtime_error("Ayaya, Miseria, vertextype problem. Stupido !");
        }

#ifdef GREEN_STYLE
        return (-U_xio1o2 * beta_ * static_cast<double>(Nc_) * factXi_ * 2.0);

#else
        //factor of 2 for Ising Spin
        return (-U_xio1o2 * beta_ * static_cast<double>(Nc_) * factXi_ * 2.0 / (((1.0 + delta_) / delta_ - 1.0) * (delta_ / (1.0 + delta_) - 1.0)));
#endif
    }

    double PhononPropagator(const double &tauIn)
    {
        double tau = tauIn;
        if (tau < 0.0)
        {
            tau += beta_;
        }

        const double w0 = Utensor.w0Phonon();

        return (-0.5 * w0 * (std::cosh((tau - beta_ / 2.0) * w0) / std::sinh(beta_ * w0 / 2.0)));
    }

    double GetDeltaTauPhonon(const double &u)
    {

        const double w0 = Utensor.w0Phonon();
        if ((std::abs(w0) < 1e-10) || (std::abs(Utensor.gPhonon()) < 1e-10))
        {
            return 1e-14;
        }
        else if (std::abs(w0 * beta_) > 25.0)
        {
            return 1e-14;
        }
        else
        {
            const double coth = 1.0 / std::tanh(w0 * beta_ / 2.0);
            const double cothfact = 4.0 * coth * coth;

            const double deltaTau = 1.0 / w0 * std::log((-2.0 + 4.0 * u + std::sqrt(16.0 * u * u + cothfact - 16.0 * u)) / (2.0 * coth - 2.0));

            if (deltaTau < 0.0)
            {
                Logging::Critical("GetDeltaTauPhonon: deltaTau < 0.0 ? ");
                return 0.0 + 1e-14;
            }
            else if (deltaTau > beta_)
            {
                Logging::Critical("GetDeltaTauPhonon: deltaTau > beta  ? ");
                return beta_ - 1e-14;
            }

            return deltaTau;
        }
    }

  private:
    const Models::UTensor Utensor; // the interaction tensor in xi=vertextype and orbital1 and orbital2 indices
    const AuxHelper auxHelper_;
    const double delta_;
    const double beta_;
    const size_t Nc_;
    const size_t NOrb_;
    const double probU_;
    const double factXi_;
    const bool isOrbitalDiagonal_;

    std::vector<double> deltaTauVec_;
};

} // namespace Diagrammatic
