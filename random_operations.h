#pragma once

#include <iostream>
#include <vector>
#include <eigen3/Eigen/Dense>
#include <NTL/ZZ.h>

using namespace NTL;


// Função para gerar aleatoriamente o par (i,j) segundo uma matriz de probabilidades (matrixProb).
// Geramos um número aleatório (obtida através de maxVal), vemos onde está posicionado na matriz e devolvemos o par (i,j+1).
//
// Exemplo:
// Suponhamos os seguintes parâmetros K = 4, D = 2 e a seguinte probabilidade: 0.449206.
// Conseguimos construir a seguinte matriz de probabilidades (matrixProb) :
//  [ 0.133333  0.066667 ]  <-- i = 0
//  [ 0.266667  0.266667 ]  <-- i = 1
//  [ 0.266667      0    ]  <-- i = 2
//    ^             ^  
//    |             |
//  j = 0         j = 1
// Como 0.133333 + 0.066667 + 0.266667 + 0.266667 >= 0.449206, 
// O par correspondente é i = 1, j = 2 e será colocado no vetor pairChosenValues nesta mesma ordem.
std::vector <int> chooseAleatoryPair (Eigen::MatrixXd &matrixProb, int K, int D, std::vector<int> &pairChosenValues, ZZ maxVal);