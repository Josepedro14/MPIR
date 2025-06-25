#pragma once

#include <eigen3/Eigen/Dense>
#include "finite_field_operations.h"
#include "print_operations.h"

// Funções Matemáticas Auxiliares


// Função para calcular o fatorial de um determinado número.
//
//  Exemplo:
//      Ao receber número 5: 
//      Calcula  5*4*3*2*1
//      E devolve o seu resultado neste caso 120
int fatorial (int n);



// Função para calcular as Combinações.
//
//  Exemplo:
//      D -> representa n e j -> k na simbologia que utilizamos em matemática C(n,k)
//      Se D = 5 e j = 2 então fazemos C(5,2) = 10 
int calculateCombinations (int D, int j);



// Função para calcular βj := D * L / C(D,j).
//
//  Exemplo:
//      Dado j = 1, D = 5 e L = 2, temos:
//      β1 = 5 * 2 / C(5,1) = 2 
// Este elemento vai ser posteriormente utilizado no cálculo da Matriz M presente na página 2 do documento elemento (1) ( mpir_functions.cpp -> função: build_MatrixM_fgAndChosePairij ).
double calculateBJForM (int j, int D, int L);



// Função recursiva que utiliza o Teorema de Laplace para calcular o determinante de uma matriz usando operações sobre um Finite Field de ordem q
//
//  Exemplo:
//     Dado um Finite Field de ordem 7, um D = 2 e uma matriz Gmatrix(D,D) = [ 1   2 ], segundo o teorema de laplace aplicado à primeira linha, obtemos:
//                                                                            [ 3   4 ]
//
//      1 x (-1)^(0+0) x det | 4 | +  2 x (-1)^(0+1) x det | 3 | 
//    = 1 x 1 x 4 + 2 x (-1) x 3
//      Utilizando as operações sobre um finite field (demonstradas no ficheiro finite_field_operations.h)
//    = (1 x 1 mod 7) x 4 + (2 x (-1) mod 7) x 3
//    = 1 x 4 mod 7 + 5 x 3 mod 7
//    = 4 + 1 mod 7
//    = 5
//
// Onde mod representa o resto da divisão inteira
int calculateDetrec (Eigen::MatrixXi Gmatrix, int D);



void calculateSubpacketsGauss (Eigen::MatrixXi A, int rows, int cols, int L, int D, int K, int symbols_subpacket);
