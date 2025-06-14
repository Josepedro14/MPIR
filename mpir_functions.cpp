#include "mpir_functions.h"
#include "finite_field_operations.h"

// Função que preenche os subpackets das mensagens com (symbols_subpacket) elementos.
void buildShuffle_Subpackets(std::vector<Message> &messages, int K, int L, int symbols_subpacket, std::mt19937 shuffle_random)
{
    // Iterar sobre as mensagens
    for(int i = 0; i < K; i++)
    {
        // Iterar sobre os subpackets da mensagem
        for(int j = 0; j < L; j++)
        {
            // Criar um subpacket para auxiliar o processo
            Subpacket new_subpacket;

            // Gerar números aleatórios do finite field de ordem q e adicioná-los ao subpacket
            for(int k = 0; k < symbols_subpacket; k++)
            {
                ZZ_p finite_field_num = random_ZZ_p();
                new_subpacket.numsR_finite_field.push_back(finite_field_num);
            }

            // Adicionar o subpacket à mensagem
            messages[i].subpackets.push_back(new_subpacket);
        }

        // Dar shuffle a um conjunto de subpackets de uma determinada mensagem
        std::shuffle(messages[i].subpackets.begin(),messages[i].subpackets.end(),shuffle_random);
    }

    std::cout << '\n';

    // Chamar função de print para as mensagens e respetivos subpacotes
    show_MessagesSubpackets(messages,K,L,symbols_subpacket);
}


// Função que calcula a matrix M, as funções f e g e uma matriz de probabilidades.
// Retorna um par(i,j) obtido aleatoriamente através das probabilidades.
std::vector <int> build_MatrixM_fgAndChosePairij (Eigen::MatrixXd &matrix_M, int K, int D, int L, ZZ maxVal)
{
    // Define uma matriz identidade (D,D) para auxiliar no cálculo das funções f e g
    Eigen::MatrixXd identityMatrix = Eigen::MatrixXd::Identity(D,D);
    // Definir um vetor de tamanho D quer irá ser preenchido por 1's e também será usado no cálculo das funções f e g
    Eigen::VectorXd colVec1s(D);

    // Iterar sobre a matriz para a preencher
    for(int line = 0; line < D; line++)
    {
        for(int col = 0; col < D; col++)
        {
            // Se estivermos na primeira linha da matriz então colocamos 1/β1
            if(line == 0)
            {
                matrix_M(line,col) = (double) (1 / calculateBJForM(1,D,L));
            }
            
            // Caso contrário estamos na segunda linha ou mais então iremos começar a preencher no formato da matriz identidade
            // Se o valor da linha atual for igual ao valor da coluna + 1, preenchemos com βline/βline+1
            else if (line == (col + 1))
            {
                matrix_M(line,col) = (double) (calculateBJForM(line,D,L) / calculateBJForM(line+1,D,L));
            }

            // Quando não pertence a nenhum dos dois acima preenchemos a zeros
            else 
            {
                matrix_M(line,col) = 0.0;
            }
        }

        // Preencher o vetor coluna com 1's
        colVec1s(line) = 1.0;
    }

    std::cout << "\nPrint Matriz M: " << '\n';
    show_Matrix(matrix_M,D,D);
    std::cout << '\n';
    show_Vector(colVec1s,D);

    // Agora iremos calcular as funções fT e gT definidas no documento na página 2, elementos (2) e (3) respetivamente
    // M ^ (K-D)
    Eigen::MatrixXd matrixKminusD = multiplyMatrixNTimes(matrix_M ,(K-D), D);
    // fT = 1^T * M ^(K-D)
    Eigen::VectorXd fT = colVec1s.transpose() * matrixKminusD;
    // M + I   
    Eigen::MatrixXd matrixMPlusI = matrix_M + identityMatrix;
    // (M+I) ^ (K-D)
    Eigen::MatrixXd matrixKminusDWthI = multiplyMatrixNTimes(matrixMPlusI, (K-D), D);
    // gT = 1^T * (M+I) ^(K-D)
    Eigen::VectorXd gT = colVec1s.transpose() * matrixKminusDWthI;

    std::cout << "\nPrint fT vector: " << '\n';
    show_Vector(fT,D);
    std::cout << "\nPrint gT vector: " << '\n';
    show_Vector(gT,D);

    // Aqui vamos tentar descobrir qual o valor do índice j* este será o índice j pertencente a [D] que maximiza fj/gj calculados anteriormente.
    int indexJ = -1;
    double maxValfDIVg = -1.0;

    // Percorremos os vetores fT e gT
    for(int j = 0; j < D; j++)
    {   
        double val = (double) (fT(j) / gT(j));

        // Se a sua divisão for maior então atualizamos o valor e o novo índice j
        if(val > maxValfDIVg)
        {
            maxValfDIVg = val;
            indexJ = j;
        }
    }

    std::cout << "\nPrint index j*: " << indexJ << " e o maxValfDIVg para esse j* é: " << maxValfDIVg << '\n';

    // Para as probabilidades e possíveis pares (i,j) temos -> (0,1);(0,2);(1,1);(1,2);(2,1);(2,2)
    Eigen::VectorXd probabilitiesPjD (D);
    // Definir matriz onde junto as probabilidades calculadas para depois sortear
    Eigen::MatrixXd matrixProb((K-D) + 1, D);

    for(int i = 0; i < D; i++)
    {
        // Se i = j* então PK-D,i = 1/gj* 
        if(i == indexJ)
        {
            probabilitiesPjD(i) = 1.0 / gT(indexJ);
        }

        // caso contrário PK-D,i = 0.0
        else 
        {
            probabilitiesPjD(i) = 0.0;
        }

        // Adicionamos à nossa matriz de probabilidades os valores calculados para a linha (K-D)
        matrixProb((K-D),i) = probabilitiesPjD(i);
    }

    std::cout << "\nMostrar probabilidades PK-D,i: " << '\n';
    show_Vector(probabilitiesPjD,D);

    std::cout << "\nMostrar Matriz Probabilidades pós PK-D,i: " << '\n';
    show_Matrix(matrixProb,(K-D) + 1, D);


    // Probabilidades Pi,D onde 0 <= i <= K-D-1
    // Para cada i calculamos a probabilidade de Pi,n onde 1 <= n <= D
    //  Pi,n = Ck-n,i * M^k-n-i * Pk-D,n
    for(int i = 0; i <= (K-D) - 1; i++)
    {   
        // Calcular M^K-D-i
        Eigen::MatrixXd matrixAux = multiplyMatrixNTimes(matrix_M, (K-D) - i, D);
        // Calcular probabilidades Pi,D
        double combinations = (double)calculateCombinations((K-D),i);
        Eigen::VectorXd probabilitiesPiD = combinations * matrixAux * probabilitiesPjD;
        
        // Passsar os valores para a matriz das probabilidades na respetiva linha
        for(int j = 0; j < D; j++)
        {
            matrixProb(i,j) = probabilitiesPiD(j);
        }

            std::cout << "\nMostrar probabilidades Pi,D na iteração i = " << i << '\n';
            show_Vector(probabilitiesPiD,D);

            std::cout << "\nMostrar Matriz Probabilidades fim na iteração i = " << i << '\n';
            show_Matrix(matrixProb,(K-D) + 1, D);
    }

     std::cout << "\nSoma de todas as probabilidades da Matriz: " << matrixProb.sum() << '\n'; 

    // Define vetor que conterá o par obtido pelas probabilidades e em que os seus 2 elementos estão dispostos pela ordem i,j
    std::vector <int> pairValues;
    
    // Chamamos a função que vai calcular o par (i,j) 
    chooseAleatoryPair(matrixProb, K, D, pairValues, maxVal);

    return pairValues;  
}


// Função para dado um determinado par, criar um i-sparse vector h de tamanho K-D com valores do finite field de ordem q, e que tem i nonzero entries em posições random. 
// Depos calculamos uma matrix G j-regular invertível (tem de ser quadrada e o determinante diferente de 0) de tamanho DxD em que G = [g1T,...,gDT] sobre o finite field de ordem q
// Regras da matriz G são as seguintes:
//     i) O vetor g1 possui exatamente j entradas não nulas em posições random
//     ii) Para cada 2 <= m <= D, as posições das entradas diferentes de 0 no vetor gm são uma rotação circular das posições das entradas não nulas do vetor gm-1 (o anterior)
void buildSparseVectorHAndG (std::vector <int> &pairIJ, int D, int K, int L, int N, int symbols_subpacket, std::mt19937 shuffle_random,std::vector<Message> &messages)
{
    // Obter o par (i,j) calculados na função acima (build_MatrixM_fgAndChosePairij)
    int i_index = pairIJ[0], j_index = pairIJ[1];
    // Construir o i-sparse vector h inicialmente preenchido com 0's
    Eigen::VectorXd h = Eigen::VectorXd::Zero(K - D);
    // Array do mesmo tamanho que o acima que vai possuir os índices de h para os baralharmos
    Eigen::VectorXi h_index (K - D);
    // Número de campos com valores do finite field de ordem q a preencher no i-sparse vector h e no vetor g1 respetivamente
    int nums_ToFill_I = i_index, nums_ToFill_J = j_index;

    // Preencher com os índices
    for(int i = 0; i < (K-D); i++)
    {
        h_index(i) = i;
    }

    // Dar shuffle no array de índices para depois colocar os elementos de forma aleatória
    std::shuffle(h_index.begin(), h_index.end(), shuffle_random);

    std::cout << "\nMostrar vetor de índices baralhado h: " << '\n';
    show_Vectorxi(h_index, K-D);

    // Preencher o i-sparse vector h com i elementos do finite field de ordem q (nums_ToFill_I == i)
    for(int j = 0; j < nums_ToFill_I; j++)
    {
        ZZ_p num_finite_field = random_ZZ_p();

        int num = conv<int>(num_finite_field);
        if(num == 0)
        {
            num = conv<int>(num_finite_field + ZZ_p(1));
        }

        // Colocar os elementos de forma aleatória
        h(h_index(j)) = num;
    }

    std::cout << "\nMostrar sparse vector-i (h) com i entradas não nulas em posições aleatórias" << '\n';
    show_Vector(h,K-D);

    // Construir matriz G j-regular e invertível de tamanho (D,D)
    Eigen::MatrixXd Gmatrix(D,D);
    
    while(true)
    {
        // Temos de construir o g1 vector de tamanho D com j entradas não nulas em posições aleatórias cada uma delas com um valor do finite field de ordem q
        Eigen::VectorXd g1vec = Eigen::VectorXd::Zero(D);
        // Array do mesmo tamanho que o acima e que vai conter os índices baralhados
        Eigen::VectorXi g1_index (D);

       // Preencher com os índices
        for(int i = 0; i < D; i++)
        {
            g1_index(i) = i;
        }  

        // Dar shuffle no array de índices para inserirmos os números gerados por Fq em posições aleatórias
        std::shuffle(g1_index.begin(),g1_index.end(),shuffle_random);

        std::cout << "\nMostrar vetor de índices baralhado g1: " << '\n';
        show_Vectorxi(g1_index, D);

        // Preencher o vector g1 pertencente á matrix G com j elementos do finite field de ordem q (nums_ToFill_J == j)
        for(int k = 0; k < nums_ToFill_J; k++)
        {
            ZZ_p num_finitefield = random_ZZ_p();

            int val= conv<int>(num_finitefield);
            if(val == 0)
            {
                val = conv<int>(num_finitefield + ZZ_p(1));
            }

            // Colocar os elementos de forma aleatória
            g1vec(g1_index(k)) = val;
        }

        std::cout << "\nMostrar g1 vector: " << '\n';
        show_Vector(g1vec, D);

        // Chamada de função para preencher a matriz (Gmatrix) de forma recursiva
        fillGMatrixInRecursive(Gmatrix,g1vec,D,0);

        std::cout << "\nPrint Gmatrix: " << '\n';
        show_Matrix(Gmatrix,D,D);

        std::cout << "\nPrint determinante Gmatrix: " << Gmatrix.determinant() << '\n';

        // Se o determinante da matriz G for maior que 0 então a matriz é invertível e podemos avançar caso contrário repetimos o processo de gerar a matriz G (Gmatrix)
        if(Gmatrix.determinant() != 0)
        {
            break;
        }

    }
    
    // Chamada de função para construir os N vectors e simular a comunicação com os servers
    constructNVectors(D,K,i_index,L,N,symbols_subpacket,shuffle_random,messages,h,Gmatrix);
}


// Função para construir N vetores de tamanho K * L.
// Estes vetores são gerados um algoritmo random, baseado no conjunto de mensagens de interesse (W), as regras são as seguintes:
//      - Se i = 0 (i -> valor do parij calculado anteriormente) então v1 é um vetor de zeros, caso contrário v1 é o vetor coeficiente 
//        correspondente a combinação linear Y1 definida como Y1 = h * [Xu1,1,...,XuK-D,1] ^T onde h é o sparse vector calculado anteriormente e u1,...,uK-D 
//        são os índices das K-D mensagens de interferência numa ordem crescente ou decrescente mas fixa.
//      - Para cada 1 <= l <= L e 1 <= m <= D o vetor v(l-1)*D+m+1 é o coeficente correspondente à combinação linear Y(l-1)*D+m+1 definida como:
//        Y(l-1)*D+m+1 := Y1 + gm * [Xw,l,...XwD,l] ^T onde gm é um vetor da Gmatrix calculada anteriormente e w1,...,wD são os índices das D mensagens de interesse, numa ordem crescente ou decrescente mas fixa.
void constructNVectors (int D, int K, int i_index, int L, int N, int symbols_subpacket, std::mt19937 shuffle_random, std::vector<Message> &messages,Eigen::VectorXd h,Eigen::MatrixXd Gmatrix)
{
    // Definir vetor que vai conter todos os N vetores de tamanho K*L que correspondem aos coeficientes das combinações lineares calculadas a seguir
    std::vector<Eigen::VectorXi> n_vectors(N, Eigen::VectorXi::Zero(K * L));
    // Vetor que irá armazenar os subpacotes 1 das (K-D) mensagens de interferência
    std::vector<Eigen::VectorXi> interference_messages;
    // Vetor que irá armazenar todos os subpacotes das D mensagens de interesse
    std::vector<Eigen::VectorXi> demand_messages;
    // Combinação linear das mensagens de interferência
    Eigen::VectorXi Y1 = Eigen::VectorXi::Zero(symbols_subpacket);
    // Vetor que vai conter todos os vetores Y(l-1)*D+m+1 a ser calculados em que 1 <= l <= L e 1 <= m <= D (também vai ter o Y1)
    std::vector <Eigen::VectorXi> Y_vectors;
    // Vetor que vai conter todos os vetores Zn-1 := Yn - Y1, onde 2 <= n <= N
    std::vector <Eigen::VectorXi> Z_vectors;


    // Criar um vetor que vai conter o sparse-vector h mas agora com elementos int em vez de double
    Eigen::VectorXi h_int = h.cast<int> ();

        // Construir um vetor de vetores Eigen::VectorXi para as interference messages e para as demand messages cada um vai conter os subpacotes correspondetes ao tipo de mensagem
        for(int i = 0; i < K; i++)
        {
            for(int j = 0; j < L; j++)
            {
                Eigen::VectorXi vecAux(symbols_subpacket);

                for(int k = 0; k < symbols_subpacket; k++)
                {
                     vecAux(k) = conv<int>(messages[i].subpackets[j].numsR_finite_field[k]);
                }

                // Armazenar nas interference messages o 1º subpacote de cada mensagem de interferência 
                if(i >= D && j == 0)
                {
                    interference_messages.push_back(vecAux);
                }

                // Armazenar em demand messages os subpacotes de cada mensagem de interesse
                else if(i < D)
                {
                    demand_messages.push_back(vecAux);
                }
            }
        }
        
        std::cout << "\nPrint interference messages (subpacket 1) messages (K-D...K): " << '\n';
        show_SubpacketsVectorXi(interference_messages,interference_messages.size(),symbols_subpacket);
        std::cout << '\n';
        std::cout << "\nPrint demand messages (subpacket 0...D) messages (0...D): " << '\n';
        show_SubpacketsVectorXi(demand_messages,demand_messages.size(),symbols_subpacket);
        std::cout << '\n';

        // Calcular o vetor Y1 que vai servir para esconder as inteções do utilizador uma vez que é construído com os elementos das mensagens que não lhe interessam (mensagens de interferência)
        // Utilizamos as funções de adição de vetores e multiplicação de vetor por um escalar sobre Fq definidas no ficheiro (finite_field_operations.cpp)
        for (int k = 0; k < K-D; k++)
        {
            Y1 = addVectorsFq(Y1, multVectorXVal(interference_messages[k], h_int(k))); 
        }

        // Construir vetor de coeficientes de Y1 (vn1), no lugar dos subpacotes 1 das mensagens 3 e 4 irá colocar os respetivos coeficientes usados para calcular Y1, tendo como base esta distribuição de pacotes no vetor neste caso ( X1,1 , X1,2 , X2,1 , X2,2 , X3,1 , X3,2 , X4,1 , X4,2 ). 
        // Imaginemos vn1: 0 0 0 0 0 0 0 0 assim no princípio então passará a estar assim: 0 0 0 0 2 0 3 0 supondo que os valores dos coeficientes pertencentes ao finite field são 2 e 3
        int index = 0;
        for(int i = D; i < K; i++)
        {
            n_vectors[0](i*L) = h_int(index);
            index++;
        }

        // Guardar a combinação linear Y1 em Y_vectors
        Y_vectors.push_back(Y1);


        int num_nvec = 0;
        // Iteramos sobre os subpacotes
        for(int l = 0; l < L; l++)
        {
            // Iteramos sobre cada vetor da matriz Gmatrix
            for(int m = 0; m < D; m++)
            {
                Eigen::VectorXi gm(D);
                Eigen::VectorXi Yvec = Eigen::VectorXi::Zero(symbols_subpacket);
                num_nvec++;

                // Aqui estamos a ir buscar o vetor gm à Gmatrix para depois usarmos na multiplicação
                for(int j = 0; j < D; j++)
                {
                    gm(j) = (int) Gmatrix(m, j);
                }

                // Percorremos as mensagens de interesse para cada uma vamos buscar o subpacote correspondente (o atual l), multiplicamos pelo elemento do vetor gm e adicionamos ao Y atual que é dado pela fórmula (l-1)*D+m+1
                // Para além disso construímos os n_vectors que faltam para os restantes Y's 
                for(int i = 0; i < D; i++)
                {
                    Eigen::VectorXi subpacketAux = demand_messages[i*L + l];
                    n_vectors[num_nvec][i * L + l] = gm(i);
   
                    Yvec = addVectorsFq(Yvec,multVectorXVal(subpacketAux,gm(i)));                    
                }

                // Adicionamos-lhe o Y1 que corresponde às mensagens de interferência
                Yvec = addVectorsFq(Yvec,Y1);
                // Guardamos o vetor Y(l-1)*D+m+1 atual em Y_vectors
                Y_vectors.push_back(Yvec);
            }
        }

        // Print Y's
        for(int i = 0; i < N; i++)
        {
            std::cout << "\nPrint Y" << i+1 << ": " << '\n';
            show_Vectorxi(Y_vectors[i],Y_vectors[i].size());
        }

        // Print Vn's
        for(int j = 0; j < N; j++)
        {
            std::cout << "\nPrint vector vn" << j+1 << ": " << '\n';
            show_Vectorxi(n_vectors[j], K*L);
        }

        // Criar uma permutação aleatória π: [N] -> [N]
        // Definir vetor de índices para baralhar de forma a manter a privacidade e esconder dos servidores as mensagens de interesse
        std::vector <int> server_indexs (N);
        for(int k = 0; k < N; k++)
        {
            server_indexs[k] = k;
        } 

        // Baralhar o vetor de índices para os servers
        std::shuffle(server_indexs.begin(),server_indexs.end(),shuffle_random);

        std::cout << "\n-------------------------------------------------COMUNICATIONS WITH SERVERS-------------------------------------------------" << '\n';

        std::cout << "\nPrint vetor server_indexs baralhado: ";
        show_VectorSTD(server_indexs,N);
        
        // Chamar função para simular conversas com os servers
        show_QuerysAndAnswersServer(server_indexs, n_vectors, Y_vectors, N, L, K, symbols_subpacket);

}