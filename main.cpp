#include "mpir_functions.h"

/*
Definição da função principal
    k -> Número de mensagens em cada Server
    N -> Número de Servers
    D -> Número de mensagens pretendidas
    q -> finite field de ordem q (q >= 3 && tem de ser primo)
    m -> tamanho de cada mensagem (acho (verificar), tem de ser par)
    L -> Grau de Subpacketização
*/
int main ()
{
    int K = 4, N = 5, D = 2, q = 5, m = 4;
    int L = (N-1) / D;
    int symbols_subpacket = m/L;

    srand(time(NULL));
    std::random_device rd;
    std::mt19937 shuffle_random(rd());

    std::vector <Message> messages(K);
    Eigen::MatrixXd matrix_M(D,D);
    std::vector <int> pairIJ;

    buildShuffle_Subpackets(messages,K,L,symbols_subpacket,q,shuffle_random);
    pairIJ = build_MatrixM_fgAndChosePairij(matrix_M,K,D,L);
    std::cout << "\nPrint index i value: " << pairIJ[0] << '\n';
    std::cout << "\nPrint index j value: " << pairIJ[1] << '\n';

    buildSparseVectorGAndVn(pairIJ, D, K, q, L, N, symbols_subpacket, shuffle_random,messages);    

    return 0;    
}