#pragma once

#include <eigen3/Eigen/Dense>
#include <NTL/ZZ_p.h>
#include <NTL/ZZ.h>

using namespace NTL;


//    Operações Elementares em Finite Fields de ordem q (Fq)
//    a,b -> são elementos do finite field de ordem q e q é a ordem do mesmo (neste caso >= 3 e primo).



// Adição em aritmética modular para um dado Finite Field de ordem q
// 
//  Para q = 5 (um número primo), os elementos pertencentes a esse finite field são: [0,1,2,3,4].
//  Um exemplo da operação de adição neste finite field de ordem 5 para os elementos 3 e 4 seria:
//          3 + 4 mod 5 -> 7 mod 5 = 2 
//  Onde 'mod' representa o resto da divisão inteira
int add_modFq (int a, int b);


// Subtração em aritmética modular para um dado Finite Field de ordem q
//
//  Para q = 7 (um número primo), os elementos pertencentes a esse finite field são: [0,1,2,3,4,5,6].
//  Um exemplo da operação de subtração neste finite field de ordem 7 para os elementos 6 e 2 seria:
//          6 - 2 mod 7 -> 4 mod 7 = 4 
//  Onde 'mod' representa o resto da divisão inteira
int sub_modFq (int a, int b);


// Multiplicação em aritmética modular para um dado Finite Field de ordem q
//
//  Para q = 3 (um número primo), os elementos pertencentes a esse finite field são: [0,1,2].
//  Um exemplo da operação de multiplicação neste finite field de ordem 3 para os elementos 1 e 2 seria:
//          1 * 2 mod 3 -> 2 mod 3 = 2 
//  Onde 'mod' representa o resto da divisão inteira
int mult_modFq (int a, int b);


// Adição em Fq para vetores Eigen::VectorXi
//
// Supondo dois vetores vec1 e vec2 do mesmo tamanho e que contém elementos de Fq (finite field de ordem q = 5).
// Esta função aplica a operação de adição (add_modFq) e retorna um novo vetor do mesmo tamanho de vec1 e vec2 com o resultado.
// Exemplo:
//   vec1 = (1, 2, 3, 4)
//   vec2 = (4, 2, 3, 1)
//   Output = (0, 4, 1, 0)
Eigen::VectorXi addVectorsFq (const Eigen::VectorXi &vec1, const Eigen::VectorXi &vec2);


// Multiplicação em Fq para vetor Eigen::VectorXi e um escalar
//
// Dado um vetor vec1 e um elemento val, ambos pertencentes ao conjunto de elementos de Fq (finite field de ordem q = 5).
// Esta função aplica a operação de multiplicação (mult_modFq) a cada elemento de vec1 com val e retorna um novo vetor do mesmo tamanho (vec1) com o resultado.
// Exemplo:
//   vec1 = (1, 2, 3, 4)
//   val = 2
//   Output = (2, 4, 1, 3)
Eigen::VectorXi multVectorXVal (const Eigen::VectorXi &vec1, int val);


// Subtração em Fq para vetores Eigen::VectorXi
//
// Supondo dois vetores vec1 e vec2 do mesmo tamanho e que contém elementos de Fq (finite field de ordem q = 5).
// Esta função aplica a operação de subtração (sub_modFq) e retorna um novo vetor do mesmo tamanho de vec1 e vec2 com o resultado.
// Exemplo:
//   vec1 = (4, 3, 2, 4)
//   vec2 = (1, 1, 2, 1)
//   Output = (3, 2, 0, 3)
Eigen::VectorXi subVectorsFq(const Eigen::VectorXi &vec1, const Eigen::VectorXi &vec2);


