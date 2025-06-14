#pragma once

#include <eigen3/Eigen/Dense>


// Função que calcula a multiplicação de uma matriz (D,D) por ela mesma num_of_times vezes.
//
// Suponhamos D = 2, num_of_times = 2 e a matrix_M seguinte: [ 2   1 ]
//                                                           [ 1   2 ]
// Exemplo:
//  matrix_M X matrix_M = [ 2   1 ]  X  [ 2   1 ]  =   [ 5   4 ]
//                        [ 1   2 ]     [ 1   2 ]      [ 4   5 ]
//                        
Eigen::MatrixXd multiplyMatrixNTimes (Eigen::MatrixXd &matrix_M, int num_of_times, int D);


// Função auxiliar usada na função seguinte. Recebe um vetor g{D-1} da Matriz G e calcula o vetor g{D} a partir do anterior.
//
//  Exemplo:
//  Imaginemos que temos g1: 0 1 0 2. 
//  Esta função realiza uma rotação circular dos elementos uma posição à direita, obtendo então: 
//      g2: 2 0 1 0 
Eigen::VectorXd buildCircularVector (Eigen::VectorXd &gDVec, int D);


// Função que preenche recursivamente a matriz Gmatrix (D,D).
// Esta função preenche a linha atual (line) e chama-se recursivamente para a próxima linha (line + 1), com o novo vetor gerado a partir de g1Vec usando a função anterior.
//
// Exemplo:
//  Para D = 3, line = 0, e um g1Vec: 0 1 2.
//
//   Devolve a seguinte matriz:
//      [ 0 1 2 ]
//      [ 2 0 1 ]
//      [ 1 2 0 ]
void fillGMatrixInRecursive (Eigen::MatrixXd &Gmatrix, Eigen::VectorXd &g1Vec, int D, int line);

