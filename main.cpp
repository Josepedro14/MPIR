#include "mpir_functions.h"


/*
Definição da função principal
    k -> Número de mensagens em cada Server
    N -> Número de Servers
    D -> Número de mensagens pretendidas
    q -> finite field de ordem q (q >= 3 && tem de ser primo)
    m -> tamanho de cada mensagem (tem de ser par)
    L -> Grau de Subpacketização
*/
int main ()
{
    int K = 4, N = 5, D = 2, m = 4;
    int L = (N-1) / D, symbols_subpacket = m/L;
    // Não pode ser mais que 31 bits (valor para o máximo para o int), pois estou a converter de ZZ_p para int nos vetores Eigen então não posso gerar um número ZZ_p maior que o máximo de int
    int max_bits = 10;

    ZZ maxVal = conv<ZZ>("1000000");

    ZZ q;
    RandomPrime(q,max_bits,10);

    ZZ_p::init(q);

    std::random_device rd;
    std::mt19937 shuffle_random(rd());

    std::vector <Message> messages(K);
    Eigen::MatrixXd matrix_M(D,D);
    std::vector <int> pairIJ;

    buildShuffle_Subpackets(messages, K, L, symbols_subpacket, shuffle_random);

    pairIJ = build_MatrixM_fgAndChosePairij(matrix_M, K, D, L, maxVal);
    std::cout << "\nPrint index i value: " << pairIJ[0] << '\n';
    std::cout << "\nPrint index j value: " << pairIJ[1] << '\n';

    buildSparseVectorGAndVn(pairIJ, D, K, L, N, symbols_subpacket, shuffle_random, messages);    

    return 0;    
}