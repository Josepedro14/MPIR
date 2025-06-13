#include "random_operations.h"


std::vector <int> chooseAleatoryPair (Eigen::MatrixXd &matrixProb, int K, int D, std::vector<int> &pairChosenValues, ZZ maxVal)
{
    // Gerar um número aleatório entre 0 e maxVal [0,maxVal) para usar nas probabilidades
    ZZ random_num = RandomBnd(maxVal);
    // Converter esse número para double (NTL não tem float nem double)
    double random = conv<double>(random_num) / conv<double>(maxVal);
    double sum = 0.0;

     std::cout << "\nProbabilidade do número gerado: " << random << '\n';

    // Percorremos a matriz de probabilidades 
    for(int i = 0; i < (K-D) + 1; i++)
    {
        for(int j = 0; j < D; j++)
        {
            // Temos o valor de soma que comparamos com o número aleatório
            sum += matrixProb(i,j);

            if(random <= sum)
            {
                // Atribuímos i,j consoante o valor de sum no momento
                pairChosenValues.push_back(i);
                pairChosenValues.push_back(j+1);
                return pairChosenValues;
            }
        }
    }

    // Valores default para i,j caso não encontre nenhum no for
    pairChosenValues.push_back(0);
    pairChosenValues.push_back(1);
    return pairChosenValues;
}