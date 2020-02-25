#include <algorithm>
#include "combinedKernel.h"


float LaplaciannnerProduct6( const float (&p)[XDIM][YDIM][ZDIM],  float (&z)[XDIM][YDIM][ZDIM])
{
    double result = 0.;
// #pragma omp parallel for 
#pragma omp parallel for reduction(+:result)
    for (int i = 1; i < XDIM - 1; i++)
        for (int j = 1; j < YDIM - 1; j++)
            for (int k = 1; k < ZDIM - 1; k++)
            {
                z[i][j][k] =
                -6 * p[i][j][k]
                + p[i+1][j][k]
                + p[i-1][j][k]
                + p[i][j+1][k]
                + p[i][j-1][k]
                + p[i][j][k+1]
                + p[i][j][k-1];
                result += (double)p[i][j][k] * (double)z[i][j][k];
            }

    return (float)result;
}

float SaxpyNorm8(const float (&z)[XDIM][YDIM][ZDIM],
                 float (&r)[XDIM][YDIM][ZDIM],
                 const float negAlpha)
{
    float result = 0.;
#pragma omp parallel for reduction(max:result)
    for (int i = 1; i < XDIM - 1; i++)
        for (int j = 1; j < YDIM - 1; j++)
            for (int k = 1; k < ZDIM - 1; k++)
            {
                r[i][j][k] = z[i][j][k] * negAlpha + r[i][j][k];
                result = std::max(result, std::abs(r[i][j][k]));
            }
    return result;
}

float CopyInnerProduct13( float (&z)[XDIM][YDIM][ZDIM], const float (&r)[XDIM][YDIM][ZDIM])
{
    double result = 0.;

#pragma omp parallel for reduction(+:result)
    for (int i = 1; i < XDIM - 1; i++)
        for (int j = 1; j < YDIM - 1; j++)
            for (int k = 1; k < ZDIM - 1; k++)
            {
                z[i][j][k] = r[i][j][k];
                result += (double)z[i][j][k] * (double)r[i][j][k];
            }

    return (float)result;
}

void SaxpyMerged(float (&p)[XDIM][YDIM][ZDIM],
                 float (&x)[XDIM][YDIM][ZDIM],
                 const float (&z)[XDIM][YDIM][ZDIM],
                 const float alpha,
                 const float beta)
{
// Should we use OpenMP parallel for here?
#pragma omp parallel for
    for (int i = 1; i < XDIM - 1; i++)
        for (int j = 1; j < YDIM - 1; j++)
            for (int k = 1; k < ZDIM - 1; k++)
            {
                x[i][j][k] = x[i][j][k] + alpha * p[i][j][k];
                p[i][j][k] = z[i][j][k] + beta * p[i][j][k];
            }
}
