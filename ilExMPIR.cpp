#include <iostream>
#include <vector>
#include <algorithm>
#include <random>
#include <Eigen/Dense>

// Definimos uma estrutura para um subpacket que é apenas um vetor de simbolos pertencentes ao Finite field de ordem q
struct Subpacket
{
    std::vector<int> numsR_finite_field;
};

// Definimos uma estrutura para as mensagens em que uma mensagem é um vetor de Subpackets
struct Message
{
    std::vector<Subpacket> subpackets;
};


// Função para mostrar as mensagens e os respetivos subpacotes
void show_MessagesSubpackets(std::vector<Message> messages, int K, int L, int symbols_subpacket)
{
    for(int i = 0; i < K; i++)
    {
        std::cout << "Message " << i+1 << ": " << '\n';

        for(int j = 0; j < L; j++)
        {
            std::cout << "Subpacket " << j+1 << ": " << '\n';

            for(int k = 0; k < symbols_subpacket; k++)
            {
                std::cout << "Value: " << messages[i].subpackets[j].numsR_finite_field[k] << '\n';
            }
        }

        std::cout << '\n';
    }
}

// Função para mostrar os valores presentes na matriz construída
void show_MatrixM (Eigen::MatrixXd matrix_Prob, int rows, int cols)
{
    std::cout << "Matriz: " << '\n';

    for(int line = 0; line < rows; line++)
    {
        for(int col = 0; col < cols; col++)
        {
            std::cout << matrix_Prob(line,col) << " ";
        }

        std::cout << '\n';
    }
}

void showVector (Eigen::VectorXd vec, int D)
{
    std::cout << "Print Vector:" << '\n';

    for(int i = 0; i < D; i++)
    {
        std::cout << vec(i) << " " << '\n';
    }
}


void buildShuffle_Subpackets(std::vector<Message> &messages, int K, int L, int symbols_subpacket, int q, std::mt19937 shuffle_random)
{
    // Preencher as mensagens (subpackets) com valores do finite field de ordem q
    for(int i = 0; i < K; i++)
    {
        for(int j = 0; j < L; j++)
        {
            Subpacket new_subpacket;

                for(int k = 0; k < symbols_subpacket; k++)
                {
                    int finite_field_num = rand() % q;
                    new_subpacket.numsR_finite_field.push_back(finite_field_num);
                }

                messages[i].subpackets.push_back(new_subpacket);
        }

        std::shuffle(messages[i].subpackets.begin(),messages[i].subpackets.end(),shuffle_random);
    }

    //show_MessagesSubpackets(messages,K,L,symbols_subpacket);
}


// Função para calcular o fatorial de um número especificado
int fatorial (int num)
{
    int res = 1;

    while(num >= 1)
    {
        res *= num;
        num--;
    }

    return res;
}


// Função para calcular as combinações
int calculateCombinacoes (int D, int j)
{
    if(j == 0 || j == D)
    {
        return 1;
    }

    return fatorial(D) / (fatorial(j) * fatorial(D-j));
}


// Função para calcular Bj para usar na matriz M
double calculateBjForM (int j, int D, int L)
{
    return (double) ((D * L) / calculateCombinacoes(D,j));
}


// Função para multiplicar matriz N vezes
Eigen::MatrixXd multiplyNMatrix (Eigen::MatrixXd &matrix_M, int num_of_times, int D)
{
    Eigen::MatrixXd matrixRes = Eigen::MatrixXd::Identity(D,D);

    for(int i = 0; i < num_of_times; i++)
    {
        matrixRes *= matrix_M;
    }

    return matrixRes;
}

std::vector<int> chooseAleatoryPair (Eigen::MatrixXd matrix_Prob, int K, int D, std::vector<int> &pairChosenValues)
{
    double random_num = (double) rand() / RAND_MAX;
    double sum=0.0;

    //std::cout << "Probabilidade: " << random_num << '\n';

    for(int i = 0; i < (K-D) + 1; i++)
    {   
        for(int j = 0; j < D; j++)
        {
            sum += matrix_Prob(i,j);

            if(random_num <= sum)
            {
                pairChosenValues.push_back(i);
                pairChosenValues.push_back(j+1);
                return pairChosenValues;
            }
        }
    }

    return pairChosenValues;
}


// Função para calcular a matriz M que irá ser usada para calcular as funções f e g utéis para as probabilidades
std::vector<int> build_MatrixM_PlusfgAndProb (Eigen::MatrixXd &matrix_M, int K, int D, int L)
{
    // Definir a matriz identidade DxD e o vetor coluna de 1's de tamanho D
    Eigen::MatrixXd identityMatriz = Eigen::MatrixXd::Identity(D,D);
    Eigen::VectorXd col_vectorof1s(D);

    for(int line = 0; line < D; line++)
    {
        for(int col = 0; col < D; col++)
        {
            if(line == 0)
            {
                matrix_M(line,col) = (double)  (1 / calculateBjForM(1,D,L));
            }

            else if(line == (col + 1))
            {
                matrix_M(line,col) = (double) (calculateBjForM(line,D,L) / calculateBjForM(line+1,D,L));
            }

            else 
            {
                matrix_M(line,col) = 0;
            }
        }

        col_vectorof1s(line) = 1;
    }

    //show_MatrixM(matrix_M,D,D);
    Eigen::MatrixXd matrixKMinusD = multiplyNMatrix(matrix_M, (K-D), D);
    Eigen::VectorXd fT = col_vectorof1s.transpose() * matrixKMinusD;
    Eigen::MatrixXd matrixMPlusIdM = matrix_M + identityMatriz;
    Eigen::MatrixXd matrixKMinusDWthIdM = multiplyNMatrix(matrixMPlusIdM, (K-D), D);
    Eigen::VectorXd gT = col_vectorof1s.transpose() * matrixKMinusDWthIdM;

    int indexJ = -1;
    double max_valuefDivg = -1.0;

    for(int k = 0; k < D; k++)
    {
        double value = (double) (fT(k) / gT(k));

       if(value  > max_valuefDivg)
       {
            max_valuefDivg = value;
            indexJ = k; 
       }
    }

    // Para as probabilidades e possíveis pares (i,j) temos -> 0 <= i <= K - D, 1 <= j <= D
    // Neste exemplo os conjuntos possíveis seriam: (0,1);(0,2);(1,1);(1,2);(2,1);(2,2)
    Eigen::MatrixXd matrix_Prob((K-D) + 1,D);
    Eigen::VectorXd probabilitiesPjD (D);
    std::vector<int> pairChosenValues;
    
    for(int i = 0; i < D; i++)
    {
        if(i == indexJ)
        {
            probabilitiesPjD(i) = (double) (1 / gT(indexJ));
        }

        else 
        {
            probabilitiesPjD(i) = 0;
        }

        matrix_Prob((K-D),i) = probabilitiesPjD(i);
    }

    for(int i = 0; i <= (K-D)-1; i++)
    {   
        Eigen::MatrixXd matrixAux = multiplyNMatrix(matrix_M,(K-D) - i,D);
        Eigen::VectorXd probabilitiesPiD = calculateCombinacoes((K-D),i) * matrixAux * probabilitiesPjD;

        for(int j = 0; j < D; j++)
        {
            matrix_Prob(i,j) = probabilitiesPiD(j);
        }
       
    }
    
    //show_MatrixM(matrix_Prob,(K-D) + 1, D);
    //std::cout << "\nSoma das probabilidades de todos os pares: " << matrix_Prob.sum() << '\n';

    // O par está na ordem i,j
    pairChosenValues = chooseAleatoryPair(matrix_Prob,K,D,pairChosenValues);

    return pairChosenValues;
}


void buildSparseVectorGAndVn(std::vector<int> pairIJ, int D, int K, int q, std::mt19937 shuffle_random)
{
    int i = pairIJ[0]; 
    std::vector<int> h(K - D, 0); 
    std::vector<int> h_index(K - D); 
    int nums_ToFill = i;
    
    for (int i = 0; i < (K - D); i++)
    {
        h_index[i] = i;
    }

    std::shuffle(h_index.begin(), h_index.end(), shuffle_random);
    
    for (int j = 0; j < nums_ToFill; j++)
    {
        int num_finite_field = rand() % q; 
        h[h_index[j]] = num_finite_field; 
    }

    for(int i = 0; i < (K-D); i++)
    {
        std::cout << h[i] << ' ';
    }

}



int main ()
{
    /*
     K -> Número de mensagens em cada Server
     N -> Número de Servers
     D -> Número de mensagens pertendidas
     q -> finite field de ordem q (q >= 3, e tem de ser primo)
     m -> tamanho de cada mensagem (talvez, tem de ser par)
     L -> grau de subpacketização
    */

    int K = 4, N = 5, D = 2, q = 5, m = 4;
    int L = (N-1) / D;
    int symbols_subpacket = m/L;

    srand(time(NULL));
    std::random_device rd;  
    std::mt19937 shuffle_random(rd());

    std::vector<Message> messages(K);
    Eigen::MatrixXd matrix_M(D,D);
    std::vector<int> pairIJ;

    buildShuffle_Subpackets(messages,K,L,symbols_subpacket,q,shuffle_random);
    pairIJ = build_MatrixM_PlusfgAndProb(matrix_M,K,D,L);

    buildSparseVectorGAndVn(pairIJ,D,K,q,shuffle_random);

    return 0; 
}