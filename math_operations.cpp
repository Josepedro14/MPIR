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
    return (double)(D * L) / (double)calculateCombinations(D,j);
}


// Função que resolve um sistema de equações lineares utilizando o método de eliminação de Gauss
void calculateMatrixGauss (Eigen::MatrixXi A, int rows, int cols, int &det, bool sol)
{   
    for(int i = 0; i < cols - 1; i++)
    {
        if(A(i,i) == 0)
        {            
            for(int j = i + 1; j < rows; j++)
            {
                if(A(j,i) != 0)
                {
                    A.row(i).swap(A.row(j));
                    break;
                }

            }
        }


        int pivo = A(i,i);

        for(int j = i+1; j < rows; j++)
        {
            if(A(j,i) != 0)
            {
                int inverse_pivoFq = inverso_modFq(pivo);
                int scalar = mult_modFq(A(j,i), inverse_pivoFq);

                for(int k = 0; k < cols; k++)
                {
                    A(j,k) = sub_modFq(A(j,k), mult_modFq(scalar, A(i,k)));
                }
            }
        }
    }


    if(cols == rows)
    {
        for(int i = 0; i < rows; i++)
        {
            for(int j = 0; j < cols; j++)
            {
                if(i == j)
                {
                    det = mult_modFq(A(i,j),det);
                }
            }
        }
    }


    if(sol)
    {
        std::vector<int> solutions(cols - 1, 0); 

        for (int i = cols - 2; i >= 0; i--) 
        {
            int soma = 0;

            for (int j = i + 1; j < cols - 1; j++)
            {
                soma = add_modFq(soma, mult_modFq(A(i, j), solutions[j]));
            }

            int rhs = sub_modFq(A(i, cols - 1), soma); 
            int pivot = A(i, i);

            if (pivot != 0)
            {
                int inv = inverso_modFq(pivot);
                solutions[i] = mult_modFq(rhs, inv);
            }
 
        }

        std::cout << "\nSoluções encontradas:\n";
        for (size_t i = 0; i < solutions.size(); ++i)
        {
            std::cout << "s" << i + 1 << " = " << solutions[i] << std::endl;
        }
    }
}

