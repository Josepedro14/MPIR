#pragma once

#include <vector>
#include <random>
#include <eigen3/Eigen/Dense>
#include "structures.h"
#include "math_operations.h"
#include "matrix_operations.h"
#include "random_operations.h"
#include "print_operations.h"


// Parâmetros globais usados nas funções abaixo:
// K -> número total de mensagens em cada servidor.
// D -> número de mensagens de interesse.
// L -> grau de Subpacketização (número de subpacotes por mensagem).
// N ->  número de servidores.
// symbols_subpackets -> número de elementos do finite field (Fq) por subpacote.
//
// Estes parâmetros têm o mesmo significado em todas as funções seguintes.



// Função para preencher o vetor mensagens com elementos do finite field de ordem q (Fq) já inicializado no main.
// Nesta função vamos percorrer todas as 'K' mensagens e para cada mensagem vamos preencher todos os 'L' subpacotes com 'symbols_subpackets' elementos aleatórios.
// Já com todos os subpacotes de todas as mensagens devidamente preenchidos vamos agora baralhar a ordem dos subpacotes de cada mensagem utilizando para isso um gerador aleatório.
//
// Parâmetros:
// std::vector<Message> messages  ->   vetor de mensagens a ser preenchido 
// std::mt19937 shuffle_random  ->  gerador aleatório, utilizado para baralhar os subpacotes
void buildShuffle_Subpackets(std::vector<Message> &messages, int K, int L, int symbols_subpacket, std::mt19937 shuffle_random);



// Função para calcular a Matriz M da página 2 do documento elemento (1).
// É uma matriz quadrada (D,D), que irá ser usada para calcular as funções f e g que irão ser usadas no cálculo de probabilidades.
// Após calculadas as probabilidades esta função retorna um par (i,j) escolhido aleatoriamente.
//
// Parâmetros:
// Matriz M  ->  Eigen::MatrixXd matrix_M
// ZZ maxVal ->  usado para o cálculo de um número aleatório entre [0,1)
std::vector <int> build_MatrixM_fgAndChosePairij (Eigen::MatrixXd &matrix_M, int K, int D, int L, ZZ maxVal);



// Função que cria:
// - um sparse-vector 'h' de tamanho (K-D) que contém i elementos não nulos aleatórios pertencentes a (Fq). 
// - uma matriz quadrada invertível 'G' (D,D) que contém j elementos não nulos aleatórios pertencentes a (Fq).
//
// Parâmetros:
// std::vector <int> pairIJ  ->  o par (i,j) calculado na função anterior e usado na construção dos dois elementos referidos acima.
// std::vector<Message> &messages  ->  as mensangens já preenchidas pela 1º função definida neste ficheiro.
// std::mt19937 shuffle_random  ->  usado para distribuir os elementos pertencentes a (Fq) de forma aleatória por 'h' e 'G'.
void buildSparseVectorHAndG (std::vector <int> &pairIJ, int D, int K, int L, int N, int symbols_subpacket, std::mt19937 shuffle_random, std::vector<Message> &messages);



// Função onde separamos as mensagens de interesse e de interferência. 
// Com as mesmas separadas construímos 'N' vetores designados por {vn} e {Yn}:
// - Y{n} -> combinação linear dos subpacotes das mensagens de (interesse ou interferência).
// - v{n} -> vetor de coeficientes duma combinação linear (consoante o valor de n).
//
// Nota:
// - Para n = 1 : os vetores acima correspondem às mensagens de interferência.
// - Para n = (l-1) * D + m + 1 : os vetores acima correspondem às mensagens de interesse, onde 1 <= l <= 'L' e 1 <= m <= 'D'.
//
//  De um conjunto 'K' de mensagens as primeiras 'D' são as de interesse e as restantes (K-D) são as de interferência. 
//
// Parâmetros:
// std::mt19937 shuffle_random  ->  gerador aleatório usado para baralhar os vetores v{n} e as respetivas respostas Y{n} pelos 'N' servidores.
// std::vector<Message> &messages  ->  as mensagens já preenchidas pela 1º função definida neste ficheiro.
// Eigen::VectorXd h  ->  sparse-vector h calculado na função anterior, utilizado para o cálculo dos 'N' vetores das mensagens de interferência.
// Eigen::MatrixXd Gmatrix  ->  matriz G também calculada na função anteiror, utilizado para o cálculo dos 'N' vetores das mensagens de interesse. 
void constructNVectors (int D, int K, int i_index, int L, int N, int symbols_subpacket, std::mt19937 shuffle_random, std::vector<Message> &messages,Eigen::VectorXd h,Eigen::MatrixXd Gmatrix);
