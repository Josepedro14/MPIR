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


std::vector<int> buildCircularVector (std::vector<int> &gDVec, int D, int lastVecValueInit, int &new_lastVecValue)
{
    std::vector<int> gDVecAux(D, 0);

    for(int i = 0; i < D; i++)
    {
        if(gDVec[i] > 0)
        {
            gDVecAux[i] = lastVecValueInit;
            lastVecValueInit = gDVec[i];
        }

        if(gDVecAux[i] > 0)
        {
            new_lastVecValue = gDVecAux[i];
        }
    }


    return gDVecAux;
}


void fillGMatrixInRecursive (Eigen::MatrixXd &Gmatrix, std::vector<int> &gDvec, int D, int lastVecValue, int line)
{

    if(line >= D)
    {
        return;
    }

    if(line == 0)
    {
        for(int col = 0; col < D; col++)
        {
            Gmatrix(line,col) = gDvec[col];
        }
        fillGMatrixInRecursive(Gmatrix,gDvec,D,lastVecValue,line + 1);
    }

    else 
    {
        int new_lastVecValue;
        std::vector<int> gDvecAux = buildCircularVector(gDvec,D,lastVecValue,new_lastVecValue);
        
        for(int col = 0; col < D; col++)
        {
            Gmatrix(line,col) = gDvecAux[col];
        }

        fillGMatrixInRecursive(Gmatrix,gDvecAux,D,new_lastVecValue,line + 1);
    }

}


void constructNVectors (int D, int K, int q, int i_index, std::mt19937 shuffle_random, std::vector<Message> &messages,std::vector<int> h,Eigen::MatrixXd Gmatrix, int L, int N, int symbols_subpacket)
{
    // Se i = 0 então v1 é preenchido apenas com zeros, caso contrário é dado por: Y1 = h * [Xu1,1,...,XuK-D,1]T onde u1,...,uK-D são os indices das K-D mensagens de interferência numa ordem crescente ou descendente mas fixa

    std::vector<std::vector<int>> n_vectors(N, std::vector<int>(K * L,0));
    std::vector<Eigen::VectorXi> demand_messages;
    std::vector<Eigen::VectorXi> interference_messages;
    Eigen::VectorXi Y1 = Eigen::VectorXi::Zero(symbols_subpacket);

    if(i_index == 0)
    {
        for(int i = 0; i < K*L; i++)
        {
            n_vectors[0][i] = 0;
        }
    }

    else 
    {

        //show_MessagesSubpackets(messages,K,L,symbols_subpacket);
        
        for(int i = 0; i < K; i++)
        {
            for(int j = 0; j < L; j++)
            {
                Eigen::VectorXi vec_aux(symbols_subpacket); 

                for(int k = 0; k < symbols_subpacket; k++)
                {
                    vec_aux(k) = messages[i].subpackets[j].numsR_finite_field[k];
                }

                if (i < (K - D)) {

                    demand_messages.push_back(vec_aux);
                } 

                else if (j == 0 && i >= (K - D)) {

                    interference_messages.push_back(vec_aux);
                }

            }
        }

        for(int i = 0; i < (K-D); i++)
        {
            Y1 += h[i] * interference_messages[i];
        }

        int index = 0;
        for(int k = (K-D); k < K; k++)
        {
            n_vectors[0][k * L] = Y1(index);
            index++;
        }

        /*
        std::cout << "\nPrint h: " << '\n';
        for(int j = 0; j < (K-D); j++)
        {
            std::cout << h[j] << ' ';
        }

        std::cout << "\nPrint Y1: " << '\n';
        for(int j = 0; j < (K-D); j++)
        {
            std::cout << Y1(j) << ' ';
        }

        std::cout << '\n';
        */
    }

    // Para 1 <= l <= L e 1 <= m <= D o vetor v(l-1)*D+m+1 := Y1 + gm * [Xw,l,...,XwD,l]T, onde w1,...,wD são os índices das mensagens requisitadas pelo utilizador
    
    std::vector <Eigen::VectorXi> yVectors;

    for(int l = 1; l < L; l++)
    {
        for(int m = 1; m < D; m++)
        {
            Eigen::VectorXi yVecX(D);
            Eigen::VectorXi yVecAux = Eigen::VectorXi::Zero(D);
            Eigen::VectorXd gmAux(D);
            std::vector <Eigen::VectorXi> demandMessagesAux;

            for(int col = 0; col < D; col ++)
            {
                gmAux(col) = Gmatrix(m-1,col);
            }

            for(size_t i = 0; i < demand_messages.size(); i++)
            {
                size_t index = i * L + (l-1);

                if(index < demand_messages.size())
                {
                    demandMessagesAux.push_back(demand_messages[index]);
                }
            }

            for(size_t j = 0; j < demandMessagesAux.size(); j++)
            {
                yVecAux += gmAux * demandMessagesAux[j]; 
            }

            yVecX = Y1 + yVecAux;
            yVectors.push_back(yVecX);
        }
    }
}


void buildSparseVectorGAndVn(std::vector<int> pairIJ, int D, int K, int q, std::mt19937 shuffle_random, std::vector<Message> &messages, int L, int N, int symbols_subpacket)
{
    int i_index = pairIJ[0], j_index = pairIJ[1]; 
    std::vector<int> h(K - D, 0); 
    std::vector<int> h_index(K - D); 
    int nums_ToFilli = i_index, nums_ToFillj = j_index, lastVecValue;

    //std::cout << "\nMostrar valor indice i: " <<  i_index <<  '\n';
    //std::cout << "\nMostrar valor indice j: " <<  j_index <<  '\n';
    
    for (int i = 0; i < (K - D); i++)
    {
        h_index[i] = i;
    }

    std::shuffle(h_index.begin(), h_index.end(), shuffle_random);
    
    for (int j = 0; j < nums_ToFilli; j++)
    {
        int num_finite_field = rand() % q;
        int num = num_finite_field != 0 ? num_finite_field : num_finite_field + 1; 
        h[h_index[j]] = num;
    }

   /*
    std::cout << "\nMostrar h sparse vector " << '\n';

    for(int i = 0; i < (K-D); i++)
    {
        std::cout << h[i] << ' ';
    }

    std::cout << '\n';
    
    */

    // Construir uma matriz j-regular (DxD) sobre Fq com duas condições:
    //      1. O vetor g1 possui exatamente j entradas não nulas em posições aleatórias
    //      2. Para cada 2 <= m <= D as posições das entradas não nulas no vetor gm são uma rotação circular das posições das entradas não nulas do vetor gm-1
    // Numa matriz regular ela tem de ser quadrada e o determinante diferente de 0

    Eigen::MatrixXd Gmatrix(D,D);
    std::vector<int> gD_index(D);
    std::vector<int> gDvec (D,0);
    int lastIndex = 0;

    for(int i = 0; i < D; i++)
    {
        gD_index[i] = i;
    }

    std::shuffle(gD_index.begin(),gD_index.end(), shuffle_random);

    for(int k = 0; k < nums_ToFillj; k++)
    {
        int num_finite_field = rand() % q;
        int num = num_finite_field != 0 ? num_finite_field : num_finite_field + 1;
        gDvec[gD_index[k]] = num;
        if(gD_index[k] >= lastIndex)
        {
            lastVecValue = num;
        }
        lastIndex = gD_index[k];
    }
    
    /*
    std::cout << "\nMostrar vector g1 " << '\n';

    for(int i = 0; i < D; i++)
    {
        std::cout << gDvec[i] << ' ';
    }

    std::cout << '\n';

    //Exemplo maior Gmatrix:

    std::vector<int> vecteste = {1,0,4,5,0};
    Eigen::MatrixXd Gmatrixteste(5,5);

    fillGMatrixInRecursive(Gmatrixteste,vecteste,5,lastVecValue,0);

    show_MatrixM(Gmatrixteste,5,5);
    */

    fillGMatrixInRecursive(Gmatrix,gDvec,D,lastVecValue,0);

    //show_MatrixM(Gmatrix,D,D);
    //std::cout << "O determinante da Gmatrix é: " << Gmatrix.determinant() << '\n'; 

    constructNVectors(D,K,q,i_index,shuffle_random, messages,h,Gmatrix,L,N,symbols_subpacket);

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

    buildSparseVectorGAndVn(pairIJ,D,K,q,shuffle_random,messages,L,N,symbols_subpacket);

    return 0; 
}