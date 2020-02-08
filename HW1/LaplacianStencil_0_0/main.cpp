#include "Timer.h"
#include "Laplacian.h"

#include <iomanip>
#include <omp.h>

int main(int argc, char *argv[])
{
    int num_workers = 1; // default value
    if (argc > 1)
    {
        num_workers = atoi(argv[1]);
        omp_set_num_threads(num_workers);
    }

    using array_t = float(&)[XDIM][YDIM];

    float *uRaw = new float[XDIM * YDIM];
    float *LuRaw = new float[XDIM * YDIM];
    array_t u = reinterpret_cast<array_t>(*uRaw);
    array_t Lu = reinterpret_cast<array_t>(*LuRaw);

    Timer timer;

    for (int test = 1; test <= 10; test++)
    {
        std::cout << "Running test iteration " << std::setw(2) << test << " ";
        timer.Start();
        ComputeLaplacian(u, Lu, num_workers);
        timer.Stop("Elapsed time : ");
    }

    return 0;
}
