#pragma once

#include <iostream>
#include <vector>
#include <eigen3/Eigen/Dense>
#include "structures.h"


// Definição de uma função para mostrar as mensagens e subpacotes pela respetiva hierarquia.
// Percorremos as mensagens todas, em cada mensagem percorremos todos os subpacotes e em cada subpacote mostramos todos os seus elementos.
void show_MessagesSubpackets (std::vector <Message> &messages, int K, int L, int symbols_subpacket);

// Função para mostrar uma matriz (double)
void show_Matrix (Eigen::MatrixXd &matrix, int rows, int cols);

// Função para mostrar um vetor (double)
void show_Vector (Eigen::VectorXd &vec, int size);

// Função para mostrar um vetor (int)
void show_Vectorxi (Eigen::VectorXi &vec, int size);

// Função para mostrar vector std::vector <int>
void show_VectorSTD (std::vector <int> &vec, int size);

// Função para mostrar subpackets das demand messages e interference messages
void show_SubpacketsVectorXi (std::vector<Eigen::VectorXi> messages, size_t size, int symbols_subpacket);

// Função para após termos preenchido e baralhado o vetor dos índices dos servers agora com esse vetor enviamos vn para os servers para simular uma query (Qn^[W]) ao server n, de seguida
// obtemos a answer (An^[W]) por parte daquele server que seria a combinação linear Y para aquele dado vetor vn.
void show_QuerysAndAnswersServer (std::vector <int> &server_indexs, std::vector <Eigen::VectorXi> &n_vectors, std::vector <Eigen::VectorXi> &Y_vectors, int N, int L, int D, int K);