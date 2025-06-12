#include "math_operations.h"


int fatorial (int n)
{
    int res = 1;

    while(n >= 1)
    {
        res *= n;
        n -= 1;
    }

    return res;
}



int calculateCombinations (int D, int j)
{
    if(j == 0 || j == D)
    {
        return 1;
    }

    return fatorial(D) / (fatorial(j) * fatorial(D-j));
}



double calculateBJForM (int j, int D, int L)
{
    return (double) ((D * L) / calculateCombinations(D,j));
}