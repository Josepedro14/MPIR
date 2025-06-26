#include <iostream>
#include <vector>
#include <random>
#include <algorithm>
#include <NTL/ZZ.h>

using namespace NTL;

// Estrutura que representa um conjunto e contém:
//  Índices das mensagens da base de dados
//  Informação sobre a presença do índice da mensagem a recuperar
struct Set
{
    int values [2];
    bool hasindexW = false;                                      
};                                                            


// Função para aplicar o XOR nos conjuntos das partições, na resposta do servidor
int applyXOR (int message1, int message2)
{
    return message1 ^ message2;
}


// Função para dar print às partições criadas
void print_partition (Set particao_aleatoria[], int K)
{
    std::cout << '\n';
    for(int i = 0; i < K/2; i++)
    {
        std::cout << "Print partition " << i+1 << " Values: " << particao_aleatoria[i].values[0] << " " << particao_aleatoria[i].values[1] << '\n';
    }
}

// Função principal que constrói a partição para os K/2 conjuntos
void build_Partition(int K, int W, int S, Set particao_aleatoria[], std::mt19937 shuffle_random)
{
    // vetor que contém os índices das mensagens presentes na BD, exceto o da mensagem a recuperar e da informação lateral
    std::vector<int> msg_indexs;

    // Preenche-lo com índices exceto W e S
    for(int i = 0; i < K; i++)
    {
        if(i != W && i != S)
        {
            msg_indexs.push_back(i);
        }
    }

    // Baralhar as posições dos índices para criar aleatoridade na junção de elementos num conjunto
    std::shuffle(msg_indexs.begin(),msg_indexs.end(),shuffle_random);

    // Definir conjunto que vai conter o índice da mensagem desejada e da informação lateral
    Set pairWS;
    // campo do conjunto usado para controlar a presença do índice W
    pairWS.hasindexW = true;
    pairWS.values[0] = W;
    pairWS.values[1] = S;

    // Dispor de forma aleatória o conjunto criado na partição e preencher os restantes conjuntos com elementos do vetor de índices de forma aleatória
    int index = 0;
    int random_index = conv<int>(RandomBnd((ZZ) (K/2)));
    for(int j = 0; j < (K/2); j++)
    {
        if(j == random_index)
        {
            particao_aleatoria[j] = pairWS;
        }

        else 
        {
            particao_aleatoria[j].values[0] = msg_indexs[(2*index)];
            particao_aleatoria[j].values[1] = msg_indexs[(2*index) + 1];
            index++;
        }
    }

    // Chamar função para printar a partição
    print_partition(particao_aleatoria,K);
}


// Função que devolve as respostas do server ao utilizador
// Para cada conjunto/par calcula o XOR das duas mensagens binárias na base de dados segundo os respetivos índices presentes no conjunto
void serverAnswer(int K, Set particao_aleatoria[], std::vector<ZZ> &answers, std::vector<ZZ> database_messages)
{
    for(int i = 0; i < (K/2); i++)
    {
        int val1 = conv<int>(database_messages[particao_aleatoria[i].values[0]]);
        int val2 = conv<int>(database_messages[particao_aleatoria[i].values[1]]);
       answers[i] = conv<ZZ>(applyXOR(val1,val2));
    }
}


// Função que recupera a mensagem pela resposta do servidor
// Verificamos com o campo hasindexW da estrutura Set se o conjunto atual contém o índice da mensagem desejada
// Se sim aplicamos o XOR à correspondente resposta do server e à informação lateral que possuímos
void recoverMsgAnswer (int K, int W, int S, Set particao_aleatoria [], std::vector<ZZ> answers, std::vector<ZZ> database_messages)
{

    for(int i = 0; i < (K/2); i++)
    {
        if(particao_aleatoria[i].hasindexW)
        {
           int answer = conv<int>(answers[i]);
           int XS = conv<int>(database_messages[S]); 
           int XW = applyXOR(answer,XS);

           std::cout << "\nValue from database " << conv<ZZ>(XW) << '\n';
           break;
        }
    }
}


int main ()
{
    /*
        Um servidor
        K --> Representa o número de mensagens na base de dados (> 0 e par)
        W --> Representa o índice da mensagem desejada pelo utilizador
        S --> Representa o índice da mensagem (informação lateral) que o utilizador já possui (W != S)
    */
    int K = 8, W = 2, S = 4;
    Set particao_aleatoria [(K/2)];
    int size_bits = 4;

    std::vector<ZZ> database_messages (K);
    std::vector<ZZ> answers (K/2);

    // Preencher as K mensagens da base de dados com valores em binário aleatórios
    for(int i = 0; i < K; i++)
    {
        database_messages[i] = RandomBits_ZZ(size_bits);
        std::cout << "\nMessage " << i+1 << ": " << database_messages[i] << '\n';
    }

    std::random_device rd;
    std::mt19937 shuffle_random(rd());

    build_Partition(K,W,S,particao_aleatoria,shuffle_random);
    serverAnswer(K,particao_aleatoria,answers,database_messages);

    // Dar print às respostas do servidor 
    std::cout << "\nResponse from server: " << '\n';
    for(size_t i = 0; i < answers.size(); i++)
    {
        std::cout << answers[i] << '\n';
    }

    recoverMsgAnswer(K,W,S,particao_aleatoria,answers,database_messages);

    return 0;
}

