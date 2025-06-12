#include "random_operations.h"


std::vector <int> chooseAleatoryPair (Eigen::MatrixXd &matrixProb, int K, int D, std::vector<int> &pairChosenValues)
{
    // Gerar um número aleatório entre 0 e 1 para usar nas probabilidades
    double random_num = (double) rand() / RAND_MAX;
    double sum = 0.0;

     std::cout << "\nProbabilidade do número gerado: " << random_num << '\n';

    // Percorremos a matriz de probabilidades 
    for(int i = 0; i < (K-D) + 1; i++)
    {
        for(int j = 0; j < D; j++)
        {
            // Temos o valor de soma que comparamos com o número aleatório
            sum += matrixProb(i,j);

            if(random_num <= sum)
            {
                // Atribuímos i,j consoante o valor de sum no momento
                pairChosenValues.push_back(i);
                pairChosenValues.push_back(j+1);
                return pairChosenValues;
            }
        }
    }

    return pairChosenValues;
}
