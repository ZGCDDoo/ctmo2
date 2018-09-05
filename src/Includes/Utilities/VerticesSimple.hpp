#pragma once

#include "Utilities.hpp"
#include "../Models/UTensorSimple.hpp"

namespace Diagrammatic
{

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
    VertexPart(const Tau_t &tau, const Site_t &site, const FermionSpin_t spin,
               const size_t orbital) : tau_(tau),
                                       site_(site),
                                       spin_(spin),
                                       orbital_(orbital) {}

    VertexPart &operator=(const VertexPart &vpart) = default;

    //Getters
    Tau_t tau() const { return tau_; };
    Site_t site() const { return site_; };
    FermionSpin_t spin() const { return spin_; };
    Orbital_t orbital() const { return orbital_; };

  private:
    Tau_t tau_;
    Site_t site_;
    FermionSpin_t spin_;
    Orbital_t orbital_;
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
        assert(data_.size() == vPartDownVec_.size());
        assert(data_.size() == vPartUpVec_.size());
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
        const size_t kkm1 = data_.size() - 1;
        std::iter_swap(data_.begin() + pp, data_.begin() + kkm1); //swap the last vertex and the vertex pp in vertices.
                                                                  //to be consistent with the updated Mup and dataCT_->Mdown_
        data_.pop_back();

        //Will need to change this when testing for different spin interactions
        std::iter_swap(vPartUpVec_.begin() + pp, vPartUpVec_.begin() + kkm1);     //swap the last vertex and the vertex pp in vertices.
        std::iter_swap(vPartDownVec_.begin() + pp, vPartDownVec_.begin() + kkm1); //swap the last vertex and the vertex pp in vertices.
        vPartUpVec_.pop_back();
        vPartDownVec_.pop_back();
    }

    /*
    Get the index for vPartUpVec corresponding to a given vertex Index (will be the same for only different spin interactions) 
    */
    // size_t GetIndexSpin()
    // {
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
};

class VertexBuilder
{
  public:
    //must hold the alphas, the values of the U, U' and (U-J_H)
    VertexBuilder(const Json &jj, const size_t &Nc) : Utensor(jj),
                                                      delta_(jj["delta"].get<double>()),
                                                      beta_(jj["beta"].get<double>()),
                                                      Nc_(Nc),
                                                      NOrb_(jj["NOrb"].get<size_t>())

    {
    }

    Vertex BuildVertex(Utilities::UniformRngMt19937_t &urng)

    {

        VertexType vertextype = (NOrb_ == 1) ? VertexType::HubbardIntra : static_cast<VertexType>(static_cast<size_t>(urng() * N_VERTEX_TYPES));

        if (vertextype == VertexType::HubbardIntra)
        {
            return BuildHubbardIntra(urng);
        }
        else if (vertextype == VertexType::HubbardInter)
        {
            return BuildHubbardInter(urng);
        }
        else if (vertextype == VertexType::HubbardInterSpin)
        {
            return BuildHubbardInterSpin(urng);
        }
        else
        {
            throw std::runtime_error("Miseria, Error in Vertices. Stupido!");
        }
    }

    Vertex BuildHubbardIntra(Utilities::UniformRngMt19937_t &urng)
    {

        VertexType vtype = VertexType::HubbardIntra;
        const Tau_t tau = urng() * beta_;
        const Site_t site = urng() * Nc_;
        const Orbital_t orbital = urng() * NOrb_;
        VertexPart vStart(tau, site, FermionSpin_t::Up, orbital);
        VertexPart vEnd(tau, site, FermionSpin_t::Down, orbital);
        AuxSpin_t aux = urng() < 0.5 ? AuxSpin_t::Up : AuxSpin_t::Down;
        return Vertex(vtype, vStart, vEnd, aux, GetKxio1o2(vtype));
    }

    Vertex BuildHubbardInter(Utilities::UniformRngMt19937_t &urng)
    {
        VertexType vtype = VertexType::HubbardInter;
        const Tau_t tau = urng() * beta_;
        const Site_t site = urng() * Nc_;
        const Orbital_t o1 = urng() * NOrb_;
        Orbital_t tmp = urng() * NOrb_;
        while (tmp == o1)
        {
            tmp = urng() * NOrb_;
        }

        const Orbital_t o2 = tmp;
        VertexPart vStart(tau, site, FermionSpin_t::Up, o1);
        VertexPart vEnd(tau, site, FermionSpin_t::Down, o2);
        AuxSpin_t aux = urng() < 0.5 ? AuxSpin_t::Up : AuxSpin_t::Down;
        return Vertex(vtype, vStart, vEnd, aux, GetKxio1o2(vtype));
    }

    Vertex BuildHubbardInterSpin(Utilities::UniformRngMt19937_t &urng)
    {
        VertexType vtype = VertexType::HubbardInter;
        const Tau_t tau = urng() * beta_;
        const Site_t site = urng() * Nc_;
        const Orbital_t o1 = urng() * NOrb_;
        Orbital_t tmp = urng() * NOrb_;
        while (tmp == o1)
        {
            tmp = urng() * NOrb_;
        }

        const Orbital_t o2 = tmp;
        const FermionSpin_t spin = urng() < 0.5 ? FermionSpin_t::Up : FermionSpin_t::Down;
        VertexPart vStart(tau, site, spin, o1);
        VertexPart vEnd(tau, site, spin, o2);
        AuxSpin_t aux = urng() < 0.5 ? AuxSpin_t::Up : AuxSpin_t::Down;
        return Vertex(vtype, vStart, vEnd, aux, GetKxio1o2(vtype));
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

        return (-U_xio1o2 * beta_ * Nc_ / (((1.0 + delta_) / delta_ - 1.0) * (delta_ / (1.0 + delta_) - 1.0)));
    }

  private:
    const Models::UTensor Utensor; // the interaction tensor in xi=vertextype and orbital1 and orbital2 indices
    const double delta_;
    const double beta_;
    const size_t Nc_;
    const size_t NOrb_;
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
