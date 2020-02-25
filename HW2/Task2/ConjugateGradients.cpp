#include "Laplacian.h"
#include "Parameters.h"
#include "PointwiseOps.h"
#include "Reductions.h"
#include "Utilities.h"
#include "Timer.h"

#include "combinedKernel.h"

#include <iostream>

extern Timer timerLaplacian;

extern Timer timerInnerProduct6;
extern Timer timerInnerProduct13;
extern Timer timerInnerProductOther;

extern Timer timerSaxpy8;
extern Timer timerSaxpy16a;
extern Timer timerSaxpy16b;
extern Timer timerSaxpyOther;

extern Timer timerNorm2;
extern Timer timerNorm8;

extern Timer timerCopy4;
extern Timer timerCopy13;

void ConjugateGradients(
    float (&x)[XDIM][YDIM][ZDIM],
    const float (&f)[XDIM][YDIM][ZDIM],
    float (&p)[XDIM][YDIM][ZDIM],
    float (&r)[XDIM][YDIM][ZDIM],
    float (&z)[XDIM][YDIM][ZDIM],
    const bool writeIterations)
{

    // Algorithm : Line 2
    timerLaplacian.Restart();ComputeLaplacian(x, z);timerLaplacian.Pause();
    timerSaxpyOther.Restart();Saxpy(z, f, r, -1);timerSaxpyOther.Pause();
    timerNorm2.Restart(); float nu = Norm(r); timerNorm2.Pause();

    // Algorithm : Line 3
    if (nu < nuMax)
        return;

    // Algorithm : Line 4
    timerCopy4.Restart();Copy(r, p);timerCopy4.Pause();
    timerInnerProductOther.Restart();float rho = InnerProduct(p, r);timerInnerProductOther.Pause();

    // Beginning of loop from Line 5
    for (int k = 0;; k++)
    {
        // std::cout << "Residual norm (nu) after " << k << " iterations = " << nu << std::endl;

        // Algorithm : Line 6
        // timerLaplacian.Restart(); ComputeLaplacian(p, z); timerLaplacian.Pause();
        // timerInnerProduct6.Restart(); float sigma = InnerProduct(p, z); timerInnerProduct6.Pause();
        timerInnerProduct6.Restart(); float sigma = LaplaciannnerProduct6(p, z); timerInnerProduct6.Pause();

        // Algorithm : Line 7
        float alpha = rho / sigma;

        // Algorithm : Line 8
        // timerSaxpy8.Restart(); Saxpy(z, r, r, -alpha); timerSaxpy8.Pause();
        // timerNorm8.Restart(); nu = Norm(r); timerNorm8.Pause();
        timerNorm8.Restart(); nu = SaxpyNorm8(z,r,-alpha); timerNorm8.Pause();

        // Algorithm : Lines 9-12
        if (nu < nuMax || k == kMax)
        {   
            timerSaxpyOther.Restart(); Saxpy(p, x, x, alpha); timerSaxpyOther.Pause();
            std::cout << "Conjugate Gradients terminated after " << k << " iterations; residual norm (nu) = " << nu << std::endl;
            // if (writeIterations) WriteAsImage("x", x, k, 0, 127);
            return;
        }

        // Algorithm : Line 13
        // timerCopy13.Restart(); Copy(r, z);timerCopy13.Pause();
        // timerInnerProduct13.Restart(); float rho_new = InnerProduct(z, r); timerInnerProduct13.Pause();
        timerInnerProduct13.Restart(); float rho_new = CopyInnerProduct13(z, r); timerInnerProduct13.Pause();

        // Algorithm : Line 14
        float beta = rho_new / rho;

        // Algorithm : Line 15
        rho = rho_new;

        // Algorithm : Line 16
        // timerSaxpy16a.Restart(); Saxpy(p, x, x, alpha); timerSaxpy16a.Pause();
        // timerSaxpy16b.Restart(); Saxpy(p, r, p, beta); timerSaxpy16b.Pause();
        timerSaxpy16b.Restart(); SaxpyMerged(p, x, r,alpha, beta); timerSaxpy16b.Pause();


        // if (writeIterations) WriteAsImage("x", x, k, 0, 127);
    }
}
