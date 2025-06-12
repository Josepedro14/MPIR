#pragma once

// Funções Matemáticas Auxiliares


/*
    Função para calcular o fatorial de um determinado número.
        Ex:
        Ao receber número 5: 
        Calcula  5*4*3*2*1
        E devolve o seu resultado neste caso 120
*/
int fatorial (int n);


/*
    Função para calcular as Combinações.
        Ex:
        D -> representa n e j -> k na simbologia que utilizamos em matemática C(n,k)
        Se D = 5 e j = 2 então fazemos C(5,2) = 10 
*/
int calculateCombinations (int D, int j);


/*
    Função para calcular Bj.
    Este elemento vai ser posteriormente utilizado no cálculo da Matriz M presente na página 2 do documento elemento (1) (função: std::vector <int> build_MatrixM_fgAndChosePairij -> mpir_functions.cpp).
*/
double calculateBJForM (int j, int D, int L);
