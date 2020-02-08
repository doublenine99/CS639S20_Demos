#pragma once

#define XDIM 16384*2
#define YDIM 16384*2

void ComputeLaplacian(const float (&u)[XDIM][YDIM], float (&Lu)[XDIM][YDIM]);
