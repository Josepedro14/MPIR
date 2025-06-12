#pragma once

#include <vector>
#include <random>
#include <eigen3/Eigen/Dense>
#include "structures.h"
#include "math_operations.h"
#include "matrix_operations.h"
#include "random_operations.h"
#include "print_operations.h"

/*
Função para preencher o atributo mensagens com elementos; 
    Mais concretamente criamos um subpacote auxiliar Subpacket new_subpacket e preenchemo-lo
    com valores aleatórios gerados a partir de um finite field de ordem q, com o subpacket 
    preenchido podemos adicioná-lo aos subpacotes de uma determinada mensagem, por fim basta-nos
    agora já com todos os subpacotes devidamente preenchidos e colocados na mensagem baralhá-los,
    para isso usamos o método shuffle em que para uma determinada mensagem passamos-lhe o ínico dos 
    subpackets o fim e o shuffle_random para ele os embaralhar.
*/
void buildShuffle_Subpackets(std::vector<Message> &messages, int K, int L, int symbols_subpacket, int q, std::mt19937 shuffle_random);



/*
Função para calcular a Matriz M da página 2 do documento elemento (1).
    É uma matriz de tamanho DxD, que irá ser usada para calcular as funções f e g úteis para a escolha das probabilidades
    Cada elemento da matriz necessita do uso da função acima definida
*/
std::vector <int> build_MatrixM_fgAndChosePairij (Eigen::MatrixXd &matrix_M, int K, int D, int L);



/*
Função para dado um determinado par, criar um i-sparse vector h de tamanho K-D com valores do finite field de ordem q, e que tem i nonzero entries em posições random. 
Depos calculamos uma matrix G j-regular invertível (tem de ser quadrada e o determinante diferente de 0) de tamanho DxD em que G = [g1T,...,gDT] sobre o finite field de ordem q
Regras da matriz G são as seguintes:
    i) O vetor g1 possui exatamente j entradas não nulas em posições random
    ii) Para cada 2 <= m <= D, as posições da entradas diferentes de 0 no vetor gm são uma rotação circular das posições das entradas não nulas do vetor gm-1 (o anterior)
*/
void buildSparseVectorGAndVn (std::vector <int> &pairIJ, int D, int K, int q, int L, int N, int symbols_subpacket, std::mt19937 shuffle_random,std::vector<Message> &messages);



/*
Função para construir N (número de servers) vetores de tamanho K * L com valores aleatórios do finite field de ordem q
    Estes vetores são gerados um algoritmo random, baseado no conjunto de mensagens de interesse (W), as regras são as seguintes:
    - Se i = 0 (i -> i_index gerado de forma aleatória anteriormente) então v1 é um vetor de zeros, caso contrário v1 é o vetor coeficiente 
    correspondente a combinação linear Y1 definida como Y1 = h * [Xu1,1,...,XuK-D,1] ^T onde h é o sparse vector calculado anteriormente e u1,...,uK-D 
    são os índices das K-D mensagens de interferência numa ordem crescente ou decrescente mas fixa.
    - Para cada 1 <= l <= L e 1 <= m <= D o vetor v(l-1)*D+m+1 é o coeficente correspondente à combinação linear Y(l-1)*D+m+1 definida como:
    Y(l-1)*D+m+1 := Y1 + gm * [Xw,l,...XwD,l] ^T onde gm é um vetor da Gmatrix calculada anteriormente e w1,...,wD são os índices das D mensagens de interesse, numa ordem crescente ou decrescente mas fixa.
    Imaginemos que temos 4 mensagens a,b,c,d então pelo exemplo ilustrativo do documento -> a,b são mensagens de interesse e c,d são mensagens de interferência
*/
void constructNVectors (int D, int K,int q, int i_index, int L, int N, int symbols_subpacket, std::mt19937 shuffle_random, std::vector<Message> &messages,Eigen::VectorXd h,Eigen::MatrixXd Gmatrix);
