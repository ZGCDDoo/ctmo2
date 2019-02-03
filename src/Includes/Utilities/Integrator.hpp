
#pragma once

#include "../../deps/cubature/cubature.h"
#include "Utilities.hpp"

namespace Integrator
{
const size_t MAXEVALS = 2000000;
const size_t KXPTS = 128;

template <typename TFct> ClusterMatrixCD_t GridKTilde(TFct fct, size_t kxpts = KXPTS)
{
    const double fact = 1.0 / (static_cast<double>(kxpts * kxpts * kxpts));
    ClusterMatrixCD_t integral(fct.n_rows(), fct.n_cols());
    integral.zeros();

    for (size_t ii = 0; ii < kxpts; ii++)
    {
        double kx = static_cast<double>(ii) / static_cast<double>(kxpts) * 2.0 * M_PI / static_cast<double>(fct.Nx);
        for (size_t jj = 0; jj < kxpts; jj++)
        {
            const double ky = static_cast<double>(jj) / static_cast<double>(kxpts) * 2.0 * M_PI / static_cast<double>(fct.Ny);
            for (size_t kk = 0; kk < kxpts; kk++)
            {
                const double kz = static_cast<double>(kk) / static_cast<double>(kxpts) * 2.0 * M_PI / static_cast<double>(fct.Nz);

                integral += fct(kx, ky, kz);
            }
        }
    }

    return (fact * integral);
}

//===========================================================================
// Cubature, general integration for full matrix functions in three dimensions
//===========================================================================

// By default, the integration is done in three dimensions.
template <typename TFct> int IntegrandCubature(unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval)
{
    // TFct est un struct ou une classe qui implémente les éléments suivants:
    // TFct fct: fct(x, y), appelle a la fct, fct.n_cols, fct.n_rows, la dimensionalité est de 3
    size_t tmp_to_shutup_warning = fdim + ndim; // to shutup unused variable.
    tmp_to_shutup_warning++;

    TFct fct = *((TFct *)fdata);

    ClusterMatrixCD_t tmpMat = fct(x[0], x[1], x[2]);

    for (size_t i = 0; i < fct.n_rows(); i++)
    {
        for (size_t j = 0; j < fct.n_cols(); j++)
        {
            fval[i + fct.n_rows() * j] = tmpMat(i, j).real();
            fval[i + fct.n_rows() * (fct.n_cols() + j)] = tmpMat(i, j).imag();
        }
    }

    return 0; // success
}

template <typename TFct>
ClusterMatrixCD_t Cubature(TFct fct, double *xmin, double *xmax, size_t maxevals = MAXEVALS, double absError = 1.49e-6,
                           double relError = 1.49e-6)
{
    const unsigned nelem = fct.n_rows() * fct.n_cols();
    double *val = new double[2 * nelem]; // for complex values
    double *err = new double[2 * nelem];

    hcubature(2 * nelem, IntegrandCubature<TFct>, &fct, 3, xmin, xmax, maxevals, absError, relError, ERROR_INDIVIDUAL, val, err);

    ClusterMatrixCD_t result(fct.n_rows(), fct.n_cols());

    for (size_t i = 0; i < fct.n_rows(); i++)
    {
        for (size_t j = 0; j < fct.n_cols(); j++)
        {
            const double real = val[i + fct.n_rows() * j];
            const double imag = val[i + fct.n_rows() * (fct.n_cols() + j)];
            result(i, j) = cd_t(real, imag);
        }
    }

    return result;
}

template <typename TFct> ClusterMatrixCD_t CubatureKTilde(TFct fct, size_t maxevals = MAXEVALS)
{

    assert(fct.n_rows() == fct.n_cols());

    double xmin[3] = {-M_PI / (fct.Nx), -M_PI / (fct.Ny), -M_PI / (fct.Nz)};
    double xmax[3] = {M_PI / (fct.Nx), M_PI / (fct.Ny), M_PI / (fct.Nz)};

    const double fact = fct.n_rows() / (8.0 * M_PI * M_PI * M_PI);
    return (fact * Cubature(fct, xmin, xmax, maxevals));
}

// =======================================
//   For DCA, do diagonal integration.
//========================================

// By default, the integration is done in three dimensions.
template <typename TFct> int IntegrandCubatureDCA(unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval)
{
    // TFct est un struct ou une classe qui implémente les éléments suivants:
    // TFct fct: fct(x, y), appelle a la fct, fct.n_cols, fct.n_rows, la dimensionalité est de 3
    size_t tmp_to_shutup_warning = fdim + ndim; // to shutup unused variable.
    tmp_to_shutup_warning++;

    TFct fct = *((TFct *)fdata);

    ClusterMatrixCD_t tmpMat = fct(x[0], x[1], x[2]);
    for (size_t i = 0; i < fct.n_rows(); i++)
    {
        fval[2 * i] = tmpMat(i, i).real();
        fval[2 * i + 1] = tmpMat(i, i).imag();
    }

    return 0; // success
}

template <typename TFct>
ClusterMatrixCD_t CubatureDCA(TFct fct, double *xmin, double *xmax, size_t maxevals = MAXEVALS, double absError = 1.49e-8,
                              double relError = 1.49e-8)
{
    const unsigned nelem = fct.n_rows(); // We only evaluate the diagonal parts, because diagonal in big K
    double *val = new double[2 * nelem]; // for complex values
    double *err = new double[2 * nelem];

    hcubature(2 * nelem, IntegrandCubatureDCA<TFct>, &fct, 3, xmin, xmax, maxevals, absError, relError, ERROR_INDIVIDUAL, val, err);

    ClusterMatrixCD_t result(fct.n_rows(), fct.n_cols());
    result.zeros();

    for (size_t i = 0; i < fct.n_rows(); i++)
    {
        const double real = val[2 * i];
        const double imag = val[2 * i + 1];
        result(i, i) = cd_t(real, imag);
    }

    return result;
}

template <typename TFct> ClusterMatrixCD_t CubatureKTildeDCA(TFct fct, size_t maxevals = MAXEVALS)
{

    assert(fct.n_rows() == fct.n_cols());

    double xmin[3] = {-M_PI / (fct.Nx), -M_PI / (fct.Ny), -M_PI / (fct.Nz)};
    double xmax[3] = {M_PI / (fct.Nx), M_PI / (fct.Ny), M_PI / (fct.Nz)};

    const double fact = fct.n_rows() / (8.0 * M_PI * M_PI * M_PI);
    return (fact * CubatureDCA(fct, xmin, xmax, maxevals));
}

} // namespace Integrator