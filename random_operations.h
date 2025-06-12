#pragma once

#include <iostream>
#include <vector>
#include <eigen3/Eigen/Dense>


// Função para gerar aleatoriamente o par i,j nesta posição no vetor, para isto usamos a matriz de probabilidades, geramos um número aleatório e vemos onde está posicionado
std::vector <int> chooseAleatoryPair (Eigen::MatrixXd &matrixProb, int K, int D, std::vector<int> &pairChosenValues);