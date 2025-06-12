#pragma once

#include <eigen3/Eigen/Dense>



// Função para calcular a multiplicar uma matriz de tamanho D*D N vezes
Eigen::MatrixXd multiplyMatrixNTimes (Eigen::MatrixXd &matrix_M, int num_of_times, int D);


/*
    Função auxiliar usada na função a seguir que recebe um vetor gD-1 da Matriz G  calcular o vetor gD a partir do anterior
        Exemplo:
        Imaginemos que temos g1: 0 1 0 2 então o que aqui estaremos a fazer é: fazer uma rotação de cada elemento uma posição para a direita ficando: 
        g2: 2 0 1 0 
*/
Eigen::VectorXd buildCircularVector (Eigen::VectorXd &gDVec, int D);


// Função que preenche recursivamente a matriz Gmatrix de tamanho (DxD) com vetores calculados a partir da função anterior
void fillGMatrixInRecursive (Eigen::MatrixXd &Gmatrix, Eigen::VectorXd &g1Vec, int D, int line);

