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
    VertexPart(const Tau_t &tau, const Site_t &site, const FermionSpin_t &spin,
               const size_t &orbital, const AuxSpin_t &aux) : tau_(tau),
                                                              site_(site),
                                                              spin_(spin),
                                                              orbital_(orbital),
                                                              superSite_{site, orbital},
                                                              aux_(aux) {}

    VertexPart &operator=(const VertexPart &vpart) = default;

    //Getters
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
           const AuxSpin_t &aux, const double &probProb) : vtype_(vtype),
                                                           vStart_(vStart),
                                                           vEnd_(vEnd),
                                                           aux_(aux),
                                                           probProb_(probProb)

    {
    }

    Vertex &operator=(const Vertex &vertex) = default;

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
    Vertices() : key_(0){};

    void AppendVertex(const Vertex &vertex)
    {
        AssertSizes();

        data_.push_back(vertex);
        verticesKeysVec_.push_back(key_);
        AppendVertexPart(vertex.vStart());

        //VertexParts differ by one for their id
        key_++;
        AppendVertexPart(vertex.vEnd());

        //Update the id number once all the vertices parts have been inserted
        key_ += 3;
        AssertSizes();
    }

    void AssertSizes()
    {
        // std::cout << "data_.size(), vPartUpVec_.size(), vPartDownVec_.size() = "
        //   << data_.size() << ", " << vPartUpVec_.size() << ", " << vPartDownVec_.size() << std::endl;
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

    void FlipAux(const size_t &p)
    {
        data_.at(p).FlipAux();
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

    void RemoveVertex(const size_t &pp)
    {
        AssertSizes();

        const VertexPart x = data_.at(pp).vStart();
        const VertexPart y = data_.at(pp).vEnd();
        const std::vector<size_t> xIndex = GetIndicesSpins(pp, x.spin());
        const std::vector<size_t> yIndex = GetIndicesSpins(pp, y.spin());

        // std::cout << "xIndex, yIndex = " << xIndex.at(0) << ", " << yIndex.at(0) << std::endl;
        if (x.spin() == FermionSpin_t::Up)
        {
            assert(x == vPartUpVec_.at(xIndex.at(0)));
            if (y.spin() == FermionSpin_t::Up)
            {
                assert(y == vPartUpVec_.at(xIndex.at(1)));
                assert(xIndex == yIndex);
            }
            else
            {
                assert(y == vPartDownVec_.at(yIndex.at(0)));
            }
        }
        else
        {
            assert(x == vPartDownVec_.at(xIndex.at(0)));
            if (y.spin() == FermionSpin_t::Down)
            {
                assert(y == vPartDownVec_.at(yIndex.at(1)));
                assert(xIndex == yIndex);
            }
            else
            {
                assert(y == vPartUpVec_.at(yIndex.at(0)));
            }
        }

        if (x.spin() != y.spin())
        {
            assert(xIndex.size() == 1);
            assert(yIndex.size() == 1);

            RemoveVertexPart(xIndex.at(0), x.spin());
            RemoveVertexPart(yIndex.at(0), y.spin());
        }
        else
        {
            assert(xIndex == yIndex);
            assert(xIndex.size() == 2);
            RemoveTwoVertexParts(xIndex, x.spin());
        }

        // CorrectIndices(pp);

        const size_t kkm1 = data_.size() - 1;
        std::iter_swap(data_.begin() + pp, data_.begin() + kkm1); //swap the last vertex and the vertex pp in vertices.
        data_.pop_back();
        std::iter_swap(verticesKeysVec_.begin() + pp, verticesKeysVec_.begin() + kkm1); //swap the last vertex and the vertex pp in vertices.
        verticesKeysVec_.pop_back();
        AssertSizes();
    }

    // void CorrectIndices(const size_t &pp)
    // {
    //     std::cout << "In Correct INdices " << std::endl;
    //     std::cout << "pp = " << pp << std::endl;
    //     Print();
    //     for (size_t ii = 0; ii < indexPartUpVec_.size(); ii++)
    //     {
    //         assert(indexPartUpVec_.at(ii) != pp); //this index should have been removed previously

    //         if (indexPartUpVec_.at(ii) > pp)
    //         {
    //             assert(indexPartUpVec_.at(ii)); //should not be zero !
    //             indexPartUpVec_.at(ii) = indexPartUpVec_.at(ii) - 1;
    //         }
    //     }

    //     for (size_t ii = 0; ii < indexPartDownVec_.size(); ii++)
    //     {
    //         assert(indexPartDownVec_.at(ii) != pp); //this index should have been removed previously

    //         if (indexPartDownVec_.at(ii) > pp)
    //         {
    //             assert(indexPartDownVec_.at(ii)); //should not be zero !
    //             indexPartDownVec_.at(ii) = indexPartDownVec_.at(ii) - 1;
    //         }
    //     }
    //     std::cout << "End Correct INdices " << std::endl;
    // }

    void RemoveVertexPart(const size_t &ppSpin, const FermionSpin_t &spin)
    {

        const size_t kkUpm1 = vPartUpVec_.size() - 1;
        const size_t kkDownm1 = vPartDownVec_.size() - 1;

        if (spin == FermionSpin_t::Up)
        {
            std::iter_swap(vPartUpVec_.begin() + ppSpin, vPartUpVec_.begin() + kkUpm1);
            std::iter_swap(indexPartUpVec_.begin() + ppSpin, indexPartUpVec_.begin() + kkUpm1);
            vPartUpVec_.pop_back();
            indexPartUpVec_.pop_back();
        }
        else
        {
            std::iter_swap(vPartDownVec_.begin() + ppSpin, vPartDownVec_.begin() + kkDownm1);
            std::iter_swap(indexPartDownVec_.begin() + ppSpin, indexPartDownVec_.begin() + kkDownm1);
            vPartDownVec_.pop_back();
            indexPartDownVec_.pop_back();
        }
    }

    void RemoveTwoVertexParts(const std::vector<size_t> &indicesToRemove, const FermionSpin_t &spin)
    {
        AssertSizes();
        // assert(false);

        // assert(false);
        // std::cout << "In removeTwoVertexParts " << std::endl;
        // std::cout << "indicesToRemove = " << indicesToRemove.at(0) << ", " << indicesToRemove.at(1) << std::endl;
        assert(indicesToRemove.size() == 2);
        const size_t kkUpm1 = vPartUpVec_.size() - 1;
        const size_t kkDownm1 = vPartDownVec_.size() - 1;

        // std::cout << "before remove " << std::endl;
        //Print();

        if (spin == FermionSpin_t::Up)
        {

            assert(indexPartUpVec_.at(indicesToRemove.at(0)) + 1 == indexPartUpVec_.at(indicesToRemove.at(1)));

            if (indicesToRemove.at(0) == kkUpm1)
            {

                std::iter_swap(vPartUpVec_.begin() + indicesToRemove.at(1), vPartUpVec_.begin() + kkUpm1 - 1);
                std::iter_swap(indexPartUpVec_.begin() + indicesToRemove.at(1), indexPartUpVec_.begin() + kkUpm1 - 1);
            }
            else if (indicesToRemove.at(0) == kkUpm1 - 1)
            {

                std::iter_swap(vPartUpVec_.begin() + indicesToRemove.at(1), vPartUpVec_.begin() + kkUpm1);
                std::iter_swap(indexPartUpVec_.begin() + indicesToRemove.at(1), indexPartUpVec_.begin() + kkUpm1);
            }
            else if (indicesToRemove.at(1) == kkUpm1)
            {
                std::iter_swap(vPartUpVec_.begin() + indicesToRemove.at(0), vPartUpVec_.begin() + kkUpm1 - 1);
                std::iter_swap(indexPartUpVec_.begin() + indicesToRemove.at(0), indexPartUpVec_.begin() + kkUpm1 - 1);
            }
            else if (indicesToRemove.at(1) == kkUpm1 - 1)
            {
                std::iter_swap(vPartUpVec_.begin() + indicesToRemove.at(0), vPartUpVec_.begin() + kkUpm1);
                std::iter_swap(indexPartUpVec_.begin() + indicesToRemove.at(0), indexPartUpVec_.begin() + kkUpm1);
            }
            else
            {
                std::iter_swap(vPartUpVec_.begin() + indicesToRemove.at(1), vPartUpVec_.begin() + kkUpm1);
                std::iter_swap(indexPartUpVec_.begin() + indicesToRemove.at(1), indexPartUpVec_.begin() + kkUpm1);

                std::iter_swap(vPartUpVec_.begin() + indicesToRemove.at(0), vPartUpVec_.begin() + kkUpm1 - 1);
                std::iter_swap(indexPartUpVec_.begin() + indicesToRemove.at(0), indexPartUpVec_.begin() + kkUpm1 - 1);
            }

            vPartUpVec_.pop_back();
            vPartUpVec_.pop_back();
            indexPartUpVec_.pop_back();
            indexPartUpVec_.pop_back();
        }
        else
        {
            assert(indexPartDownVec_.at(indicesToRemove.at(0)) + 1 == indexPartDownVec_.at(indicesToRemove.at(1)));

            if (indicesToRemove.at(0) == kkDownm1)
            {

                std::iter_swap(vPartDownVec_.begin() + indicesToRemove.at(1), vPartDownVec_.begin() + kkDownm1 - 1);
                std::iter_swap(indexPartDownVec_.begin() + indicesToRemove.at(1), indexPartDownVec_.begin() + kkDownm1 - 1);
            }
            else if (indicesToRemove.at(0) == kkDownm1 - 1)
            {

                std::iter_swap(vPartDownVec_.begin() + indicesToRemove.at(1), vPartDownVec_.begin() + kkDownm1);
                std::iter_swap(indexPartDownVec_.begin() + indicesToRemove.at(1), indexPartDownVec_.begin() + kkDownm1);
            }
            else if (indicesToRemove.at(1) == kkDownm1)
            {
                std::iter_swap(vPartDownVec_.begin() + indicesToRemove.at(0), vPartDownVec_.begin() + kkDownm1 - 1);
                std::iter_swap(indexPartDownVec_.begin() + indicesToRemove.at(0), indexPartDownVec_.begin() + kkDownm1 - 1);
            }
            else if (indicesToRemove.at(1) == kkDownm1 - 1)
            {
                std::iter_swap(vPartDownVec_.begin() + indicesToRemove.at(0), vPartDownVec_.begin() + kkDownm1);
                std::iter_swap(indexPartDownVec_.begin() + indicesToRemove.at(0), indexPartDownVec_.begin() + kkDownm1);
            }
            else
            {
                std::iter_swap(vPartDownVec_.begin() + indicesToRemove.at(1), vPartDownVec_.begin() + kkDownm1);
                std::iter_swap(indexPartDownVec_.begin() + indicesToRemove.at(1), indexPartDownVec_.begin() + kkDownm1);

                std::iter_swap(vPartDownVec_.begin() + indicesToRemove.at(0), vPartDownVec_.begin() + kkDownm1 - 1);
                std::iter_swap(indexPartDownVec_.begin() + indicesToRemove.at(0), indexPartDownVec_.begin() + kkDownm1 - 1);
            }

            vPartDownVec_.pop_back();
            vPartDownVec_.pop_back();
            indexPartDownVec_.pop_back();
            indexPartDownVec_.pop_back();
        }

        // std::cout << "After remove " << std::endl;

        // Print();
        // std::cout << "End removeTwoVertexParts " << std::endl;
    }

    /*
    Get the index for vPartUpVec and vPartDownVec corresponding to a given vertex Index (will be the same for only different spin interactions) 
    */
    std::vector<size_t>
    GetIndicesSpins(const size_t &pp, const FermionSpin_t &spin) const
    {
        const size_t vertexKey = verticesKeysVec_.at(pp);
        const VertexPart x = data_.at(pp).vStart();
        const VertexPart y = data_.at(pp).vEnd();
        //Print();
        std::vector<size_t> indices;
        if (spin == FermionSpin_t::Up)
        {
            //Find in order the vertexParts corresponding to the same vertex
            for (size_t ii = 0; ii < indexPartUpVec_.size(); ii++)
            {
                if (vertexKey == indexPartUpVec_.at(ii))
                {
                    indices.push_back(ii);
                }
            }

            for (size_t ii = 0; ii < indexPartUpVec_.size(); ii++)
            {
                if (vertexKey + 1 == indexPartUpVec_.at(ii))
                {
                    indices.push_back(ii);
                }
            }
        }
        else
        {
            for (size_t ii = 0; ii < indexPartDownVec_.size(); ii++)
            {
                if (vertexKey == indexPartDownVec_.at(ii))
                {
                    indices.push_back(ii);
                }
            }

            for (size_t ii = 0; ii < indexPartDownVec_.size(); ii++)
            {
                if (vertexKey + 1 == indexPartDownVec_.at(ii))
                {
                    indices.push_back(ii);
                }
            }
        }

        x.spin() != y.spin() ? assert(indices.size() == 1) : assert(indices.size() == 2);
        return indices;
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
            std::iter_swap(indexPartUpVec_.begin() + ii, indexPartUpVec_.begin() + jj);
        }
        else
        {
            std::iter_swap(vPartDownVec_.begin() + ii, vPartDownVec_.begin() + jj);
            std::iter_swap(indexPartDownVec_.begin() + ii, indexPartDownVec_.begin() + jj);
        }
    }

    // vois Save(const std::string& filename)
    // {
    //     for(size_t ii=0;ii<data_.size(); ii++)
    //     {
    //         std::cout << "o1, o2, spin1, spin2, tau1, tau2 = " <<
    //     }

    // }

    void SaveConfig(const std::string &fname)
    {
        std::ofstream fout(fname);
        for (size_t ii = 0; ii < size(); ii++)
        {
            const auto x = data_.at(ii).vStart();
            const auto y = data_.at(ii).vEnd();
            assert(x.tau() == y.tau());
            assert(x.aux() == y.aux());

            fout << static_cast<int>(x.aux()) << " " << x.site() << " " << x.tau() << " " << static_cast<int>(x.spin()) << " " << static_cast<int>(y.spin()) << " " << x.orbital() << " " << y.orbital() << " " << std::endl;
        }
    }

    //Getters
    size_t size() const { return data_.size(); };
    size_t NUp() const { return vPartUpVec_.size(); };
    size_t NDown() const { return vPartDownVec_.size(); };
    std::vector<size_t> verticesKeysVec() const { return verticesKeysVec_; }; // Each vertex has a unqique key identifying it

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
    std::vector<size_t> indexPartDownVec_; // Ex: The row and col 0 of Ndown_ is associtated to the vertexPart of the vertex given by the id of indexPartDownVec_.at(0)
    std::vector<size_t> verticesKeysVec_;  // Each vertex has a unqique key identifying it
    size_t key_;

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
                                                          NOrb_ * NOrb_ * 2 * 2 / 2 - NOrb_ // Pauli principale and dont double count pairs
                                                          ),
                                                      isOrbitalDiagonal_(jj["IsOrbitalDiagonal"].get<bool>())

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

        if ((o1 == o2) && (spin1 != spin2))
        {
            vertextype = VertexType::HubbardIntra;
            const VertexPart vStart(tau, site, FermionSpin_t::Up, o1, aux);
            const VertexPart vEnd(tau, site, FermionSpin_t::Down, o2, aux);
            // std::cout << "GetKxio1o2(HubbardIntra) = " << GetKxio1o2(vertextype) << std::endl;

            return Vertex(vertextype, vStart, vEnd, aux, GetKxio1o2(vertextype));
        }
        else if ((o1 != o2) && (spin1 != spin2))
        {
            // std::cout << "Here !" << std::endl;
            vertextype = VertexType::HubbardInter;
            const VertexPart vStart(tau, site, FermionSpin_t::Up, o1, aux);
            const VertexPart vEnd(tau, site, FermionSpin_t::Down, o2, aux);
            // std::cout << "GetKxio1o2(HubbardInter) = " << GetKxio1o2(vertextype) << std::endl;
            return Vertex(vertextype, vStart, vEnd, aux, GetKxio1o2(vertextype));
        }
        else if ((o1 != o2) && (spin1 == spin2))
        {
            vertextype = VertexType::HubbardInterSpin;
            const VertexPart vStart(tau, site, spin1, o1, aux);
            const VertexPart vEnd(tau, site, spin2, o2, aux);
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
            U_xio1o2 = isOrbitalDiagonal_ ? 0.0 : Utensor.UPrime();
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
                U_xio1o2 = (Utensor.UPrime() - Utensor.JH());
            }
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
    const bool isOrbitalDiagonal_;
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
