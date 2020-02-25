#include "PointwiseOps.h"
#include <algorithm>

void Copy(const float (&x)[XDIM][YDIM][ZDIM], float (&y)[XDIM][YDIM][ZDIM])
{
#pragma omp parallel for
    for (int i = 1; i < XDIM - 1; i++)
        for (int j = 1; j < YDIM - 1; j++)
            for (int k = 1; k < ZDIM - 1; k++)
                y[i][j][k] = x[i][j][k];
}

void Saxpy(const float (&x)[XDIM][YDIM][ZDIM],
           const float (&y)[XDIM][YDIM][ZDIM],
           float (&z)[XDIM][YDIM][ZDIM],
           const float scale)
{
    // Should we use OpenMP parallel for here?
    for (int i = 1; i < XDIM - 1; i++)
        for (int j = 1; j < YDIM - 1; j++)
            for (int k = 1; k < ZDIM - 1; k++)
                z[i][j][k] = x[i][j][k] * scale + y[i][j][k];
}

// void SaxpyMerged(float (&p)[XDIM][YDIM][ZDIM],
//                  float (&x)[XDIM][YDIM][ZDIM],
//                  const float (&z)[XDIM][YDIM][ZDIM],
//                  const float alpha,
//                  const float beta)
// {
// // Should we use OpenMP parallel for here?
// #pragma omp parallel for
//     for (int i = 1; i < XDIM - 1; i++)
//         for (int j = 1; j < YDIM - 1; j++)
//             for (int k = 1; k < ZDIM - 1; k++)
//             {
//                 x[i][j][k] = x[i][j][k] + alpha * p[i][j][k];
//                 p[i][j][k] = z[i][j][k] + beta * p[i][j][k];
//             }
// }

// float SaxpyNorm8(const float (&z)[XDIM][YDIM][ZDIM],
//                 float (&r)[XDIM][YDIM][ZDIM],
//                 const float negAlpha)
// {
//     float result = 0.;
// #pragma omp parallel for reduction(max:result)
//     for (int i = 1; i < XDIM - 1; i++)
//         for (int j = 1; j < YDIM - 1; j++)
//             for (int k = 1; k < ZDIM - 1; k++)
//             {
//                 r[i][j][k] = z[i][j][k] * negAlpha + r[i][j][k];
//                 result = std::max(result, std::abs(r[i][j][k]));
//             }
//     return result;
// }