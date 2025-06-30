#include "mpir_functions.h"


// Definição da função principal
//    K -> Número de mensagens em cada Server
//    N -> Número de Servers
//    D -> Número de mensagens pretendidas
//    q -> finite field de ordem q (q >= 3 e primo)
//    m -> tamanho de cada mensagem ( m >= 2 e par)
//    L -> Grau de Subpacketização (número de subpacotes)
// 155 520
int main ()
{
    int K = 4, N = 5, D = 2, m = 4;
    int L = (N-1) / D, symbols_subpacket = m/L;
    // Não pode ser mais que 31 bits (valor máximo para o int), pois estou a converter de ZZ_p para int nos vetores Eigen::VectorXi então não posso gerar um número ZZ_p maior que o máximo de int
    int max_bits = 3;

    // Valor auxiliar para gerar número aleatório no intervalo [0,1) na função chooseAleatoryPair (random_operations.cpp)
    ZZ maxVal = conv<ZZ>(1000000);

    // Ordem q do finite field
    ZZ q;   

    // Gerar um número primo aleatório q enquanto for < 3
    do{
        RandomPrime(q,max_bits,10);
    }while (conv<int>(q) < 3);
    
    // Inicializar o finite_field com a ordem q
    ZZ_p::init(q);

    std::cout << "\nPrint ordem do finite field: " << conv<int>(q) << '\n';

    // Declarar e inicializar gerador aleatório
    std::random_device rd;
    std::mt19937 shuffle_random(rd());

    // Criar vetor de mensagens, matriz M e vetor para armazenar par (i,j)
    std::vector <Message> messages(K);
    Eigen::MatrixXd matrix_M(D,D);
    std::vector <int> pairIJ;

    buildShuffle_Subpackets(messages, K, L, symbols_subpacket, shuffle_random);

    pairIJ = build_MatrixM_fgAndChosePairij(matrix_M, K, D, L, maxVal);
    //std::cout << "\nPrint index i value: " << pairIJ[0] << '\n';
    //std::cout << "\nPrint index j value: " << pairIJ[1] << '\n';

    buildSparseVectorHAndG(pairIJ, D, K, L, N, symbols_subpacket, shuffle_random, messages);    

    return 0;    
}