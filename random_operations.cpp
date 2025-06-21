#include "random_operations.h"


void chooseAleatoryPair (Eigen::MatrixXd &matrixProb, int K, int D, std::vector<int> &pairValues, ZZ maxVal)
{
    // Gerar um número aleatório entre 0 e maxVal [0,maxVal)
    ZZ random_num = RandomBnd(maxVal);
    // Obtemos um valor aleatório entre [0,1) através da divisão do valor gerado com o passado em parâmetro
    double random = conv<double>(random_num) / conv<double>(maxVal);
    // variável de controlo para escolha do par (i,j) consoante o 'random' gerado
    double sum = 0.0;

     std::cout << "\nProbabilidade do número gerado: " << random << '\n';

    // Percorremos a matriz de probabilidades 
    for(int i = 0; i < (K-D) + 1; i++)
    {
        for(int j = 0; j < D; j++)
        {
            
            sum += matrixProb(i,j);

            // Comparamos o valor de sum com random se for maior ou igual significa que já passamos, então na posição do último elemento que adicionámos
            // vamos buscar o valor da linha e colocamos em i e vamos buscar o valor da coluna + 1 e colocamos em j. 
            if(sum >= random)
            {
                pairValues.push_back(i);
                pairValues.push_back(j+1);
            }
        }
    }

    // Valores default para par (i,j) caso nenhum seja encontrado no for.
    pairValues.push_back(0);
    pairValues.push_back(1);
}