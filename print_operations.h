#pragma once

#include <iostream>
#include <vector>
#include <eigen3/Eigen/Dense>
#include "structures.h"


// Função para mostrar mensagens e subpacotes pela respetiva hierarquia.
// Percorre todas as mensagens, para cada mensagem, percorre todos os subpacotes e mostra os seus elementos.
void show_MessagesSubpackets (std::vector <Message> &messages, int K, int L, int symbols_subpacket);

// Função para mostrar uma matriz com valores 'double'.
// Recebe como parâmetros a matriz (matrix) e as suas dimensões linhas (rows) e colunas (cols).
void show_Matrix (Eigen::MatrixXd &matrix, int rows, int cols);

// Função para mostrar um vetor com valores 'double' (Eigen::VectorXd).
// Recebe como parâmetros o vetor (vec) e o seu tamanho (size).
void show_Vector (Eigen::VectorXd &vec, int size);

// Função para mostrar um vetor com valores 'int' (Eigen::VectorXi).
// Recebe como parâmetros o vetor (vec) e o seu tamanho (size).
void show_Vectorxi (Eigen::VectorXi &vec, int size);

// Função para mostrar vector com valores 'int' (std::vector <int>).
// Recebe como parâmetros o vetor (vec) e o seu tamanho (size).
void show_VectorSTD (std::vector <int> &vec, int size);

// Função para mostrar os subpacotes das mensagens de interesse e de interferência.
// Recebe como parâmetros um vetor (std::vector) de vetores (Eigen::VectorXi) messages, o seu tamanho (size) e o número de elementos por subpacote (symbols_subpacket).
void show_SubpacketsVectorXi (std::vector<Eigen::VectorXi> messages, size_t size, int symbols_subpacket);

// Função para mostrar os vetores (n_vectors) enviados como query (Qn^[W]) para cada servidor e as respetivas respostas (Y_vectors) como answer (An^[W]).
// Recebe como parâmetros um vetor (server_indexs) com os índices baralhados de forma aleatória, os (n_vectors) e (Y_vectors) para os servidores e algumas variáveis auxiliares.
// N -> Número de Servidores
// L -> Grau de Subpacketização
// K -> Número de mensangens em cada servidor
// symbols_subpacket -> número de elementos por subpacote
void show_QuerysAndAnswersServer (std::vector <int> &server_indexs, std::vector <Eigen::VectorXi> &n_vectors, std::vector <Eigen::VectorXi> &Y_vectors, int N, int L, int K, int symbols_subpacket);