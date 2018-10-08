#pragma once

#include "Utilities.hpp"
#include "../Models/UTensorSimple.hpp"
#include <boost/math/special_functions/binomial.hpp>

namespace Diagrammatic
{

enum class VertexType
{
    HubbardIntra,     //Hubbard intraorbital
    HubbardInter,     // Hubbard interorbital, different spins (U')
    HubbardInterSpin, // Hubbard interorbital same spin (U'-J_H)
    //Phonon,
    Invalid
};

const size_t N_VERTEX_TYPES = 3; //for now, we dont do Phonon
const size_t INVALID = 999;

class VertexPart
{
  public:
    VertexPart(const Tau_t &tau, const Site_t &site, const FermionSpin_t spin,
               const size_t orbital) : tau_(tau),
                                       site_(site),
                                       spin_(spin),
                                       orbital_(orbital),
                                       superSite_{site, orbital} {}

    VertexPart &operator=(const VertexPart &vpart) = default;

    //Getters
    Tau_t tau() const { return tau_; };
    Site_t site() const { return site_; };
    FermionSpin_t spin() const { return spin_; };
    Orbital_t orbital() const { return orbital_; };
    SuperSite_t superSite() const { return superSite_; };

  private:
    Tau_t tau_;
    Site_t site_;
    FermionSpin_t spin_;
    Orbital_t orbital_;
    SuperSite_t superSite_;
};

class Vertex
{

  public:
    Vertex(const VertexType &vtype, const VertexPart &vStart, const VertexPart &vEnd,
           const AuxSpin_t &aux, const double &probProb) : vtype_(vtype),
                                                           vStart_(vStart),
                                                           vEnd_(vEnd),
                                                           aux_(aux),
                                                           probProb_(probProb)

    {
    }

    const Vertex &operator=(const Vertex &vertex)
    {
        if (this == &vertex)
            return *this; //évite les boucles infinies
        vStart_ = vertex.vStart_;
        vEnd_ = vertex.vEnd_;
        aux_ = vertex.aux_;

        return *this;
    }

    // Getters
    AuxSpin_t aux() const { return aux_; };
    VertexType vtype() const { return vtype_; };
    double probProb() const { return probProb_; };
    VertexPart vStart() const { return vStart_; };
    VertexPart vEnd() const { return vEnd_; };

    //Setters
    void SetAux(AuxSpin_t aux)
    {
        aux_ = aux;
    }

    void FlipAux() { aux_ == AuxSpin_t::Up ? aux_ = AuxSpin_t::Down : aux_ = AuxSpin_t::Up; };

    double Ising()
    {
        if (aux_ == AuxSpin_t::Zero)
        {
            return 0.0;
        }
        return (aux_ == AuxSpin_t::Up ? 1.0 : -1.0);
    }

  private:
    VertexType vtype_;
    VertexPart vStart_;
    VertexPart vEnd_;
    AuxSpin_t aux_;
    double probProb_;
};

class Vertices
{

  public:
    Vertices(){};

    void AppendVertex(const Vertex &vertex)
    {
        data_.push_back(vertex);
        AppendVertexPart(vertex.vStart());
        AppendVertexPart(vertex.vEnd());
        assert(2 * data_.size() == (vPartUpVec_.size() + vPartDownVec_.size()));
        // assert(data_.size() == vPartDownVec_.size());
        // assert(data_.size() == vPartUpVec_.size());
    }

    void FlipAux(const size_t &p)
    {
        data_.at(p).FlipAux();
    }

    void AppendVertexPart(const VertexPart &vPart)
    {
        (vPart.spin() == FermionSpin_t::Up) ? vPartUpVec_.push_back(vPart) : vPartDownVec_.push_back(vPart);
    }

    void RemoveVertex(const size_t &pp)
    {
        const VertexPart x = data_.at(pp).vStart();
        const VertexPart y = data_.at(pp).vEnd();
        const size_t xIndex = GetIndicesSpins(pp, x.spin());
        const size_t yIndex = GetIndicesSpins(pp, y.spin());

        RemoveVertexPart(xIndex, x.spin());
        RemoveVertexPart(yIndex, y.spin());

        data_.erase(data_.begin() + pp);
    }

    void RemoveVertexPart(const size_t &ppSpin, const FermionSpin_t &spin)
    {
        if (spin == FermionSpin_t::Up)
        {
            vPartUpVec_.erase(vPartUpVec_.begin() + ppSpin);
        }
        else
        {
            vPartDownVec_.erase(vPartDownVec_.begin() + ppSpin);
        }
    }

    /*
    Get the index for vPartUpVec and vPartDownVec corresponding to a given vertex Index (will be the same for only different spin interactions) 
    */
    size_t GetIndicesSpins(const size_t &pp, const FermionSpin_t &spin) const
    {
        size_t indexVertexPartUp = 0;
        size_t indexVertexPartDown = 0;

        for (size_t ii = 0; ii < pp; ii++) //Not sure here ...
        {
            const VertexPart x = data_.at(ii).vStart();
            const VertexPart y = data_.at(ii).vEnd();
            x.spin() == FermionSpin_t::Up ? indexVertexPartUp++ : indexVertexPartDown++;
            y.spin() == FermionSpin_t::Up ? indexVertexPartUp++ : indexVertexPartDown++;
        }

        assert(indexVertexPartUp + indexVertexPartDown == 2 * (pp));

        return ((spin == FermionSpin_t::Up) ? indexVertexPartUp : indexVertexPartDown);
    }

    void Swap(const size_t &ii, const size_t &jj)
    {

        std::iter_swap(data_.begin() + ii, data_.begin() + jj);
    }

    void SwapSpin(const size_t &ii, const size_t &jj, const FermionSpin_t &spin)
    {
        if (spin == FermionSpin_t::Up)
        {
            std::iter_swap(vPartUpVec_.begin() + ii, vPartUpVec_.begin() + jj);
        }
        else
        {
            std::iter_swap(vPartDownVec_.begin() + ii, vPartDownVec_.begin() + jj);
        }
    }

    // vois Save(const std::string& filename)
    // {
    //     for(size_t ii=0;ii<data_.size(); ii++)
    //     {
    //         std::cout << "o1, o2, spin1, spin2, tau1, tau2 = " <<
    //     }

    // }

    //Getters
    size_t size() const { return data_.size(); };
    size_t NUp() const { return vPartUpVec_.size(); };
    size_t NDown() const { return vPartDownVec_.size(); };

    Vertex at(const size_t &i) const { return data_.at(i); };

    VertexPart atUp(const size_t &i) { return vPartUpVec_.at(i); };
    VertexPart atDown(const size_t &i) { return vPartDownVec_.at(i); };

  private:
    std::vector<Vertex> data_;
    std::vector<VertexPart> vPartUpVec_;
    std::vector<VertexPart> vPartDownVec_;
}; // namespace Diagrammatic

class VertexBuilder
{
  public:
    //must hold the alphas, the values of the U, U' and (U-J_H)
    VertexBuilder(const Json &jj, const size_t &Nc) : Utensor(jj),
                                                      delta_(jj["delta"].get<double>()),
                                                      beta_(jj["beta"].get<double>()),
                                                      Nc_(Nc),
                                                      NOrb_(jj["NOrb"].get<size_t>()),
                                                      factXi_(
                                                          NOrb_ //* NOrb_ * 2 * 2 / 2 - NOrb_ // Pauli principale and dont double count pairs
                                                      )

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

        while ((o1 == o2) && (spin1 == spin2))
        {
            o1 = urng() * NOrb_;
            o2 = urng() * NOrb_;
            spin1 = (urng() < 0.5) ? FermionSpin_t::Up : FermionSpin_t::Down;
            spin2 = (urng() < 0.5) ? FermionSpin_t::Up : FermionSpin_t::Down;
        }

        VertexType vertextype = VertexType::Invalid;

        if ((o1 == o2) && spin1 != spin2)
        {
            vertextype = VertexType::HubbardIntra;
            const VertexPart vStart(tau, site, FermionSpin_t::Up, o1);
            const VertexPart vEnd(tau, site, FermionSpin_t::Down, o2);
            // std::cout << "GetKxio1o2(HubbardIntra) = " << GetKxio1o2(vertextype) << std::endl;

            return Vertex(vertextype, vStart, vEnd, aux, GetKxio1o2(vertextype));
        }
        else if ((o1 != o2) && (spin1 != spin2))
        {
            // std::cout << "Here !" << std::endl;
            vertextype = VertexType::HubbardInter;
            const VertexPart vStart(tau, site, FermionSpin_t::Up, o1);
            const VertexPart vEnd(tau, site, FermionSpin_t::Down, o2);
            // std::cout << "GetKxio1o2(HubbardInter) = " << GetKxio1o2(vertextype) << std::endl;
            return Vertex(vertextype, vStart, vEnd, aux, GetKxio1o2(vertextype));
        }
        else if ((o1 != o2) && (spin1 == spin2))
        {
            vertextype = VertexType::HubbardInterSpin;
            const VertexPart vStart(tau, site, spin1, o1);
            const VertexPart vEnd(tau, site, spin2, o2);
            return Vertex(vertextype, vStart, vEnd, aux, GetKxio1o2(vertextype));
        }
        else
        {
            throw std::runtime_error("Miseria, Error in Vertices. Stupido!");
        }
    }

    double GetKxio1o2(const VertexType &vtype)
    {

        double U_xio1o2 = INVALID;

        if (vtype == VertexType::HubbardIntra)
        {
            U_xio1o2 = Utensor.U();
        }
        else if (vtype == VertexType::HubbardInter)
        {
            U_xio1o2 = Utensor.UPrime();
        }
        else if (vtype == VertexType::HubbardInterSpin)
        {
            U_xio1o2 = Utensor.JH();
        }
        else
        {
            throw std::runtime_error("Ayaya, Miseria, vertextype problem. Stupido !");
        }

#ifdef GREEN_STYLE
        return (-U_xio1o2 * beta_ * Nc_ * factXi_);

#else
        return (-U_xio1o2 * beta_ * Nc_ * factXi_ / (((1.0 + delta_) / delta_ - 1.0) * (delta_ / (1.0 + delta_) - 1.0)));
#endif
    }

  private:
    const Models::UTensor Utensor; // the interaction tensor in xi=vertextype and orbital1 and orbital2 indices
    const double delta_;
    const double beta_;
    const size_t Nc_;
    const size_t NOrb_;
    const double factXi_;
};

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

    double FAux(const FermionSpin_t &spin, const AuxSpin_t &aux) const
    {
        if (aux == AuxSpin_t::Zero)
        {
            return 1.0;
        }
        return (auxValue(spin, aux) / (auxValue(spin, aux) - 1.0));
    };

    double gamma(const FermionSpin_t &spin, const AuxSpin_t &auxI, const AuxSpin_t &auxJ) const //little gamma
    {
        const double fsJ = FAux(spin, auxJ);
        return ((FAux(spin, auxI) - fsJ) / fsJ);
    }

    double delta() const { return delta_; };

  private:
    const double delta_;
};

} // namespace Diagrammatic
