#pragma once

#include "Parameters.h"

// Copy array x into y
void Copy(const float (&x)[XDIM][YDIM][ZDIM], float (&y)[XDIM][YDIM][ZDIM]);

// Scale array x by given number, add y, and write result into z
void Saxpy(const float (&x)[XDIM][YDIM][ZDIM], const float (&y)[XDIM][YDIM][ZDIM],
           float (&z)[XDIM][YDIM][ZDIM], const float scale);
// void SaxpyMerged(float (&p)[XDIM][YDIM][ZDIM],
//                  float (&x)[XDIM][YDIM][ZDIM],
//                  const float (&z)[XDIM][YDIM][ZDIM],
//                  const float alpha,
//                  const float beta);
// float SaxpyNorm8(const float (&z)[XDIM][YDIM][ZDIM],
//                  float (&r)[XDIM][YDIM][ZDIM],
//                  const float negAlpha);