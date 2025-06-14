#pragma once


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
