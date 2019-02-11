#pragma once

// External Libraries
#include <armadillo>
#include "../../deps/nlohmann_json/json.hpp"

using Json = nlohmann::json;

// Inspired by Patrick Sémon

using cd_t = std::complex<double>;
using Sign_t = int;
using Site_t = size_t;
using Tau_t = double;
using Orbital_t = size_t;
using UInt64_t = size_t; // unsigned long long int;

using SuperSite_t = std::pair<size_t, size_t>; // site, then orbital number

enum class AuxSpin_t
{
    Up,
    Down,
    Zero
};

enum class FermionSpin_t
{
    Up,
    Down
};

using SiteVectorCD_t = arma::cx_vec;
using SiteRowCD_t = arma::cx_rowvec;
using ClusterSitesCD_t = std::vector<arma::cx_vec>;
using ClusterMatrixCD_t = arma::cx_mat;
using ClusterCubeCD_t = arma::cx_cube;

using SiteVector_t = arma::vec;
using SiteRow_t = arma::rowvec;
using ClusterSites_t = std::vector<arma::vec>;
using ClusterMatrix_t = arma::mat;
using ClusterCube_t = arma::cube;

using UpdStats_t = std::map<std::string, std::valarray<size_t>>;

////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////

namespace Utilities
{


std::string GetSpinName(const FermionSpin_t &spin) { return (spin == FermionSpin_t::Up ? "Up" : "Down"); }

size_t GetIndepOrbitalIndex(const size_t &o1, const size_t &o2, const size_t &NOrb)
{

    assert(o1 < NOrb);
    assert(o2 < NOrb);

    size_t indepOrbitalIndex = 0;
    const std::pair<size_t, size_t> pairTarget = o1 < o2 ? std::make_pair(o1, o2) : std::make_pair(o2, o1);

    for (Orbital_t nu1 = 0; nu1 < NOrb; ++nu1)
    {
        for (Orbital_t nu2 = nu1; nu2 < NOrb; ++nu2)
        {

            if (pairTarget == std::make_pair(nu1, nu2))
            {
                return indepOrbitalIndex;
            }
            indepOrbitalIndex++;
        }
    }

    assert(o1 == 0);
    assert(o2 == 0);
    return 0;
}

} // namespace Utilities
