#pragma once
#include "Parameters.h"
float LaplaciannnerProduct6(const float (&p)[XDIM][YDIM][ZDIM], float (&z)[XDIM][YDIM][ZDIM]);
void SaxpyMerged(float (&p)[XDIM][YDIM][ZDIM],
                 float (&x)[XDIM][YDIM][ZDIM],
                 const float (&z)[XDIM][YDIM][ZDIM],
                 const float alpha,
                 const float beta);
float SaxpyNorm8(const float (&z)[XDIM][YDIM][ZDIM],
                 float (&r)[XDIM][YDIM][ZDIM],
                 const float negAlpha);
float CopyInnerProduct13(float (&z)[XDIM][YDIM][ZDIM],
                         const float (&r)[XDIM][YDIM][ZDIM]);