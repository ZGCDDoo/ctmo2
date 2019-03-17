#pragma once

#include "ctmo/Foundations/GreenMat.hpp"

namespace Fourier
{

cd_t MatToTau(const SiteVectorCD_t &greenMat, const double &tau,
              const double &beta) // Only for a "scalar green function, not a cluster green"
{
    cd_t greenTau = 0.0;
    const cd_t im{0.0, 1.0};

    for (size_t n = 0; n < greenMat.size(); n++)
    {
        const double wn = (2.0 * n + 1.0) * M_PI / beta;
        greenTau += std::exp(-im * wn * tau) * greenMat(n);
        greenTau += std::exp(im * wn * tau) * std::conj(greenMat(n));
    }

    return (greenTau / beta);
}

cd_t MatToTauAnalytic(SiteVectorCD_t greenMat, const double &tau, const double &beta, const cd_t &fm, const cd_t &sm,
                      const cd_t &tm)
{

    cd_t result{0.0, 0.0};

    // result+= les moments en tau calculés analytiquement
    result += -0.5 * fm;
    result += (tau / 2.0 - beta / 4.0) * sm;
    result += -1.0 / 4.0 * (tau * (tau - beta)) * tm;

    // On transforme la greenMat moins ses moments
    for (size_t n = 0; n < greenMat.n_elem; n++)
    {
        const cd_t iwn(0.0, (2.0 * n + 1.0) * M_PI / beta);
        greenMat(n) -= fm / iwn + sm / (iwn * iwn) + tm / (iwn * iwn * iwn);
    }

    result += MatToTau(greenMat, tau, beta);

    return result;
}

ClusterMatrixCD_t MatToTauCluster(const GreenMat::GreenCluster0Mat &greenCluster0Mat, const double &tau)
{

    ClusterCubeCD_t dataMat = greenCluster0Mat.data();
    const double beta = greenCluster0Mat.beta();
    ClusterMatrixCD_t result(dataMat.n_rows, dataMat.n_cols);
    result.zeros();

    // result+= les moments en tau calculés analytiquement
    result += -0.5 * greenCluster0Mat.fm();
    result += (tau / 2.0 - beta / 4.0) * greenCluster0Mat.sm();
    result += -1.0 / 4.0 * (tau * (tau - beta)) * greenCluster0Mat.tm();

    // On transforme la greenMat moins ses moments
    for (size_t n = 0; n < dataMat.n_slices; n++)
    {
        cd_t iwn(0.0, (2.0 * n + 1.0) * M_PI / beta);
        dataMat.slice(n) -= greenCluster0Mat.fm() / (iwn) + greenCluster0Mat.sm() / (iwn * iwn) +
                            greenCluster0Mat.tm() / (iwn * iwn * iwn);
    }

    for (size_t i = 0; i < dataMat.n_rows; i++)
    {
        for (size_t j = 0; j < dataMat.n_cols; j++)
        {
            result(i, j) += MatToTau(dataMat.tube(i, j), tau, beta);
        }
    }

    return result;
}

} // namespace Fourier