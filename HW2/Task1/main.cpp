#include "ConjugateGradients.h"
#include "Timer.h"
#include "Utilities.h"

Timer timerLaplacian;

Timer timerInnerProduct6;
Timer timerInnerProduct13;
Timer timerInnerProductOther;

Timer timerNorm2;
Timer timerNorm8;

Timer timerSaxpy8;
Timer timerSaxpy16a;
Timer timerSaxpy16b;
Timer timerSaxpyOther;

Timer timerCopy4;
Timer timerCopy13;

int main(int argc, char *argv[])
{
    using array_t = float(&)[XDIM][YDIM][ZDIM];

    float *xRaw = new float[XDIM * YDIM * ZDIM];
    float *fRaw = new float[XDIM * YDIM * ZDIM];
    float *pRaw = new float[XDIM * YDIM * ZDIM];
    float *rRaw = new float[XDIM * YDIM * ZDIM];
    float *zRaw = new float[XDIM * YDIM * ZDIM];

    array_t x = reinterpret_cast<array_t>(*xRaw);
    array_t f = reinterpret_cast<array_t>(*fRaw);
    array_t p = reinterpret_cast<array_t>(*pRaw);
    array_t r = reinterpret_cast<array_t>(*rRaw);
    array_t z = reinterpret_cast<array_t>(*zRaw);

    // Initialization
    {
        Timer timer;
        timer.Start();
        InitializeProblem(x, f);
        timer.Stop("Initialization : ");
    }
    Timer totalConjugateGradients;
    // Call Conjugate Gradients algorithm
    timerLaplacian.Reset();
    timerInnerProduct6.Reset();
    timerInnerProduct13.Reset();
    timerInnerProductOther.Reset();
    timerSaxpy8.Reset();
    timerSaxpy16a.Reset();
    timerSaxpy16b.Reset();
    timerSaxpyOther.Reset();
    timerNorm2.Reset();
    timerNorm8.Reset();
    timerCopy4.Reset();
    timerCopy13.Reset();

    totalConjugateGradients.Start();
    ConjugateGradients(x, f, p, r, z);
    totalConjugateGradients.Stop("Total ConjugateGradients time is : ");

    timerLaplacian.Print("Total Laplacian Time : ");

    timerInnerProductOther.Print("InnerProduct() on line 4, : Time = ");
    timerInnerProduct6.Print("InnerProduct() on line 6 : Time = ");
    timerInnerProduct13.Print("InnerProduct() on line 13 : Time = ");

    timerSaxpyOther.Print("Saxpy() on line 2, 9-12 : Time = ");
    timerSaxpy8.Print("Saxpy() on line 8 : Time = ");
    timerSaxpy16a.Print("1st Saxpy() on line 16 : Time = ");
    timerSaxpy16b.Print("2nd Saxpy() on line 16 : Time = ");

    timerNorm2.Print("Norm() on line 2 : Time = ");
    timerNorm8.Print("Norm() on line 8 : Time = ");

    timerCopy4.Print("Copy() on line 4 : Time = ");
    timerCopy13.Print("Copy() on line 13 : Time = ");


    return 0;
}
