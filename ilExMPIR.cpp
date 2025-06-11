#include <iostream>
#include <vector>
#include <algorithm>
#include <random>
#include <eigen3/Eigen/Dense>

/* 
Definir uma estrutura para subpacotes;
   Um subpacote é um elemento que contém um conjunto de números pertencentes ao finite field de ordem q
   Então ele pode ser representado como um vetor de inteiros
*/
struct Subpacket
{
    std::vector <int> numsR_finite_field;
};


/*
Definir uma estrutura de uma mensagem;
    Uma mensagem é um conjunto de subpacotes em que cada subpacote como descrito acima é um conjunto de m/L símbolos
    - L -> representa o grau de supacketização (número de pacotes em que é subdividido cada mensagem)
    - m -> tamanho de cada mensagem (nº de elementos (acho eu confirmar))    
*/
struct Message
{
    std::vector <Subpacket> subpackets;
};

/*
Definição de uma função para mostrar as mensagens e subpacotes pela respetiva hierarquia
    Basicamente percorremos as mensagens todas em cada mensagem percorremos todos os subpacotes e em cada subpacote mostramos todos os seus elementos
*/
void show_MessagesSubpackets (std::vector <Message> &messages, int K, int L, int symbols_subpacket)
{
    // Iterar sobre as mensagens 
    for(int i = 0; i < K; i++)
    {
        std::cout << "Message " << i+1 << ": " << '\n';
        //  Iterar sobre os subpacotes
        for(int j = 0; j < L; j++)
        {
            std::cout << "Subpacket " << j+1 << ": " << '\n';
            // Iterar sobre os elementos dos subpacotes
            for(int k = 0; k < symbols_subpacket; k++)
            {
                std::cout << "Value: " << messages[i].subpackets[j].numsR_finite_field[k] << '\n';
            }
        }

        std::cout << '\n';
    }
}



// Função para mostrar uma matriz
void show_Matrix (Eigen::MatrixXd &matrix, int rows, int cols)
{
    std::cout << "Matriz:" << '\n';

    for(int line = 0; line < rows; line++)
    {
        for(int col = 0; col < cols; col++)
        {
            std::cout << matrix(line,col) << " ";
        }

        std::cout << '\n';
    }
}


// Função para mostrar um vetor
void show_Vector (Eigen::VectorXd &vec, int size)
{
    std::cout << "Vector: " << '\n';

    for(int i = 0; i < size; i++)
    {
        std::cout << vec(i) << " " << '\n';
    }
}


// Função para mostrar um vetor
void show_Vectorxi (Eigen::VectorXi &vec, int size)
{
    std::cout << "Vector: " << '\n';

    for(int i = 0; i < size; i++)
    {
        std::cout << vec(i) << " ";
    }

    std::cout << '\n';
}


// Função para mostrar vector std::vector <int>
void show_VectorSTD (std::vector <int> &vec, int size)
{
    std::cout << "\nVector STD: " << '\n';

    for(int i = 0; i < size; i++)
    {
        std::cout << vec[i] << " ";
    }

    std::cout << '\n';
}


// Função para mostrar subpackets das demand messages e interference messages
void show_SubpacketsVectorXi (std::vector<Eigen::VectorXi> messages, size_t size, int symbols_subpacket)
{
     for(size_t i = 0; i < messages.size(); i++)
        {
            std::cout << "\nPrint subpacotes das mensagens: " << '\n';
            for(int j = 0; j < symbols_subpacket; j++)
            {
                std::cout << messages[i](j) << " ";
            }
        }
}


/*
    Função para após termos preenchido e baralhado o vetor dos índices dos servers agora com esse vetor enviamos vn para os servers para simular uma query (Qn^[W]) ao server n, de seguida
    obtemos a answer (An^[W]) por parte daquele server que seria a combinação linear Y para aquele dado vetor vn.
*/
void show_QuerysAndAnswersServer (std::vector <int> &server_indexs, std::vector <Eigen::VectorXi> &n_vectors, std::vector <Eigen::VectorXi> &Y_vectors, int N, int L, int D, int K)
{
    for(int i = 0; i < N; i++)
    {
        std::cout << "\n--------------------------------------------------------- SERVER " << i+1 << " ---------------------------------------------------------\n";
        // Aqui estamos a enviar para o server indicado o vetor vn como o Query (Qn^[W]) ao server (n -> server_indexs[i])
        std::cout << "\nQuerry: " << '\n';
        show_Vectorxi(n_vectors[server_indexs[i]], K*L);
        // Aqui o server (n -> server_indexs[i]) vai devolver/calcular uma Answer (An^[W]) que é a combinação linear Yn para o vetor vn recebido
        std::cout << "\nAnswer: " << '\n';
        show_Vectorxi(Y_vectors[server_indexs[i]], D);
        std::cout << "\n----------------------------------------------------------------------------------------------------------------------------\n";
    }
}


/*
Função para calcular o fatorial de um determinado número
    Ex:
    Ao receber número 5: 
    Calcula  5*4*3*2*1
    E devolve o seu resultado
*/
int fatorial (int n)
{
    int res = 1;

    while(n >= 1)
    {
        res *= n;
        n -= 1;
    }

    return res;
}


// Função para calcular as Combinações em Matemática
int calculateCombinations (int D, int j)
{
    if(j == 0 || j == D)
    {
        return 1;
    }

    return fatorial(D) / (fatorial(j) * fatorial(D-j));
}


/*
Função para calcular Bj;
    Este elemento vai ser posteriormente utilizado no cálculo da Matriz M presente na página 2 do documento elemento (1)
*/
double calculateBJForM (int j, int D, int L)
{
    return (double) ((D * L) / calculateCombinations(D,j));
}


/*
Função para calcular a multiplicar uma matriz de tamanho D*D N vezes
*/
Eigen::MatrixXd multiplyMatrixNTimes (Eigen::MatrixXd &matrix_M, int num_of_times, int D)
{
    Eigen::MatrixXd matrixRes = Eigen::MatrixXd::Identity(D,D);

    for(int i = 0; i < num_of_times; i++)
    {
        matrixRes *= matrix_M;
    }

    return matrixRes;
}


/*
Função para gerar aleatoriamente o par i,j nesta posição no vetor, para isto usamos
a matriz de probabilidades, geramos um número aleatório e vemos onde está posicionado
*/
std::vector <int> chooseAleatoryPair (Eigen::MatrixXd &matrixProb, int K, int D, std::vector<int> &pairChosenValues)
{
    // Gerar um número aleatório entre 0 e 1 para usar nas probabilidades
    double random_num = (double) rand() / RAND_MAX;
    double sum = 0.0;

     std::cout << "\nProbabilidade do número gerado: " << random_num << '\n';

    // Percorremos a matriz de probabilidades 
    for(int i = 0; i < (K-D) + 1; i++)
    {
        for(int j = 0; j < D; j++)
        {
            // Temos o valor de soma que comparamos com o número aleatório
            sum += matrixProb(i,j);

            if(random_num <= sum)
            {
                // Atribuímos i,j consoante o valor de sum no momento
                pairChosenValues.push_back(i);
                pairChosenValues.push_back(j+1);
                return pairChosenValues;
            }
        }
    }

    return pairChosenValues;
}



/*
Função que recebe um vetor gD-1 da Matriz G  calcular o vetor gD a partir do anterior
*/
Eigen::VectorXd buildCircularVector (Eigen::VectorXd &gDVec, int D)
{
    // Vetor auxiliar que vai conter a rotação das posições
    Eigen::VectorXd gDVecAux = Eigen::VectorXd::Zero(D);

   /*
   Exemplo:
   Imaginemos que temos g1: 0 1 0 2 então o que aqui estaremos a fazer é: fazer uma rotação de cada elemento uma posição para a direita ficando: 
   g2: 2 0 1 0 
   */

    for(int i = 0; i < D; i++)
    {
        int position = i + 1;

        if(position > D-1)
        {
            position = 0;
        }

        gDVecAux(position) = gDVec(i);
    }

    return gDVecAux;
}



void fillGMatrixInRecursive (Eigen::MatrixXd &Gmatrix, Eigen::VectorXd &g1Vec, int D, int line)
{
    // Se a linha atual for maior que o tamanho da Matriz então retorna a Matriz obtida
    if(line >= D)
    {
        return;
    }

    // Se a linha for a primeira vamos usar o vetor g1 com j entradas não nulas que geramos anteriormente para preencher na primeira linha da Matriz
    if(line == 0)
    {
        for(int col = 0; col < D; col++)
        {
            Gmatrix(line,col) = g1Vec(col);
        }

        // Preencher o resto da Matriz recursivamente com os vetores circulares (gD) gerados a partir do anterior (gD-1)
        fillGMatrixInRecursive(Gmatrix, g1Vec, D, line + 1);
    }

    // Se não for nenhum dos casos acima então vamos preencher o resto da matriz com o auxílio da função buildCircularVector descrita acima
    else 
    {   
        // Vetor circular gD obtido através do anterior gD-1
        Eigen::VectorXd gDvecAux = buildCircularVector(g1Vec, D);

        // Preencher a próxima linha da matriz (line) com os valores do vetor gerado 
        for(int col = 0; col < D; col++)
        {
            Gmatrix(line,col) = gDvecAux(col);
        }

        // Preencher o resto da Matriz recursivamente ou devolvê-la preenchida
        fillGMatrixInRecursive(Gmatrix, gDvecAux, D, line + 1);
    }
}


/*
Função para preencher o atributo mensagens com elementos; 
    Mais concretamente criamos um subpacote auxiliar Subpacket new_subpacket e preenchemo-lo
    com valores aleatórios gerados a partir de um finite field de ordem q, com o subpacket 
    preenchido podemos adicioná-lo aos subpacotes de uma determinada mensagem, por fim basta-nos
    agora já com todos os subpacotes devidamente preenchidos e colocados na mensagem baralhá-los,
    para isso usamos o método shuffle em que para uma determinada mensagem passamos-lhe o ínico dos 
    subpackets o fim e o shuffle_random para ele os embaralhar.
*/
void buildShuffle_Subpackets(std::vector<Message> &messages, int K, int L, int symbols_subpacket, int q, std::mt19937 shuffle_random)
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
                int finite_field_num = rand() % q;
                new_subpacket.numsR_finite_field.push_back(finite_field_num);
            }

            // Adicionar o subpacket ao conjunto de subpackets de uma mensagem
            messages[i].subpackets.push_back(new_subpacket);
        }

        // Dar shuffle a um conjunto de subpackets de uma determinada mensagem
        std::shuffle(messages[i].subpackets.begin(),messages[i].subpackets.end(),shuffle_random);
    }

    // Chamar função de print para as mensagens e respetivos subpacotes
    show_MessagesSubpackets(messages,K,L,symbols_subpacket);
}



/*
Função para calcular a Matriz M da página 2 do documento elemento (1)
    É uma matriz de tamanho DxD, que irá ser usada para calcular as funções f e g úteis para a escolha das probabilidades
    Cada elemento da matriz necessita do uso da função acima definida
*/
std::vector <int> build_MatrixM_fgAndChosePairij (Eigen::MatrixXd &matrix_M, int K, int D, int L)
{
     // Define uma matriz identidade DxD para auxilidar no cálculo das funções f e g
    Eigen::MatrixXd identityMatrix = Eigen::MatrixXd::Identity(D,D);
    // Definir um vetor de tamanho D quer irá ser preenchido por 1's e também será usado no cálculo das funções f e g
    Eigen::VectorXd colVec1s(D);

    // Iterar sobre a matriz para a preencher
    for(int line = 0; line < D; line++)
    {
        for(int col = 0; col < D; col++)
        {
            // Se estivermos na primeira linha da matriz então colocamos 1/B1
            if(line == 0)
            {
                matrix_M(line,col) = (double) (1 / calculateBJForM(1,D,L));
            }
            
            // Caso contrário estamos na segunda linha ou mais então iremos começar a preencher no formato da matriz identidade
            // Ou seja se o valor da linha atual for igual ao valor da coluna + 1, preenchemos com Bline/Bline+1
            else if (line == (col + 1))
            {
                matrix_M(line,col) = (double) (calculateBJForM(line,D,L) / calculateBJForM(line+1,D,L));
            }

            // Quando não pertence a nenhum dos dois acima preenche a zeros
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

    for(int i = 0; i < D; i++)
    {
        double val = (double) (fT(i) / gT(i));

        if(val > maxValfDIVg)
        {
            maxValfDIVg = val;
            indexJ = i;
        }
    }

    std::cout << "\nPrint index j*: " << indexJ << " e o maxValfDIVg para esse j* é: " << maxValfDIVg << '\n';

    // Para as probabilidades e possíveis pares (i,j) temos -> 
    // Neste exemplo os conjuntos possíveis seriam: (0,1);(0,2);(1,1);(1,2);(2,1);(2,2)
    Eigen::VectorXd probabilitiesPjD (D);
    /// Definir matriz onde junto as probabilidades calculadas para depois sortear
    Eigen::MatrixXd matrixProb((K-D) + 1, D);

    for(int i = 0; i < D; i++)
    {
        // Se i = j* então Pk-d,i = 1/gj* 
        if(i == indexJ)
        {
            probabilitiesPjD(i) = (double) (1 / gT(indexJ));
        }

        // caso contrário Pk-d,i = 0.0
        else 
        {
            probabilitiesPjD(i) = 0.0;
        }

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
        // Calcular M^k-d-i
        Eigen::MatrixXd matrixAux = multiplyMatrixNTimes(matrix_M, (K-D) - i, D);
        //
        Eigen::VectorXd probabilitiesPiD = calculateCombinations((K-D),i) * matrixAux * probabilitiesPjD;
        
        // Passsar os valores para a matriz das probabilidades
        for(int j = 0; j < D; j++)
        {
            matrixProb(i,j) = probabilitiesPiD(j);
        }

            std::cout << "\nMostrar probabilidades Pi,D na iteração i = " << i << '\n';
            show_Vector(probabilitiesPiD,D);

            std::cout << "\nMostrar Matriz Probabilidades fim na iteração i = " << i << '\n';
            show_Matrix(matrixProb,(K-D) + 1, D);
    }

     std::cout << "\nSoma das probabilidades de todos os pares: " << matrixProb.sum() << '\n'; 

    // Define vetor que conterá o par obtido pelas probabilidades e está disposto pela ordem i,j
    std::vector <int> pairChosenValues;
    
    pairChosenValues = chooseAleatoryPair(matrixProb, K, D, pairChosenValues);

    return pairChosenValues;  
}



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
void constructNVectors (int D, int K,int q, int i_index, int L, int N, int symbols_subpacket, std::mt19937 shuffle_random, std::vector<Message> &messages,Eigen::VectorXd h,Eigen::MatrixXd Gmatrix)
{

    std::vector<Eigen::VectorXi> n_vectors(N, Eigen::VectorXi::Zero(K * L));
    std::vector<Eigen::VectorXi> interference_messages;
    std::vector<Eigen::VectorXi> demand_messages;
    Eigen::VectorXi Y1 = Eigen::VectorXi::Zero(symbols_subpacket);
    // Array que vai conter todos os vetores Y(l-1)*D+m+1 a ser calculados em que 1 <= l <= L e 1 <= m <= D 
    std::vector <Eigen::VectorXi> Y_vectors;
    // Array que vai conter todos os Zn vetores em que 2 <= n <= N onde Z1,...ZN-1 são combinações lineares independentes dos D("2")*L("2") subpacotes das mensagens de interesse
    std::vector <Eigen::VectorXi> Z_vectors;

    Eigen::VectorXi h_int = h.cast<int> ();

        // Construir um vetor de Eigen::VectorXi para as interference messages e para as demand messages cada um vai  conter os subpacotes correspondetes ao tipo de mensagem
        for(int i = 0; i < K; i++)
        {
            for(int j = 0; j < L; j++)
            {
                Eigen::VectorXi vecAux(symbols_subpacket);

                for(int k = 0; k < symbols_subpacket; k++)
                {
                    vecAux(k) = messages[i].subpackets[j].numsR_finite_field[k];
                }

                // Armazenar nas interference messages o 1º subpacote de cada mensagem de interferência 
                if(i >= K-D && j == 0)
                {
                    interference_messages.push_back(vecAux);
                }

                // Armazenar em demand messages os subpacotes de cada mensagem requerida
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

        // Calcular o vetor Y1 que vai servir para esconder as inteções do utilizador uma vez que é construído com os elementos 
        // das mensagens que não lhe interessam (mensagens de interferência)
        for (int k = 0; k < K-D; k++)
        {
            Y1 += h_int(k) * interference_messages[k]; 
        }

        // Construir vetor de coeficientes de Y1 (vn1), no lugar dos subpacotes 1 das mensagens 3 e 4 irá colocar os respetivos coeficientes usados para calcular Y1, tendo como base esta distribuição de pacotes no vetor neste caso ( X1,1 , X1,2 , X2,1 , X2,2 , X3,1 , X3,2 , X4,1 , X4,2 ). 
        // Imaginemos vn1: 0 0 0 0 0 0 0 0 assim no princípio então passará a estar assim: 0 0 0 0 2 0 3 0 supondo que os valores dos coeficientes pertencentes ao finite field são 2 e 3
        int index = 0;
        for(int i = (K-D); i < K; i++)
        {
            n_vectors[0](i*L) = h_int(index);
            index++;
        }

        // Guardar a combinação linear Y1 no array de Y_vectors
        Y_vectors.push_back(Y1);


        int num_nvec = 0;
        // Iteramos sobre os subpacotes
        for(int l = 0; l < L; l++)
        {
            // Iteramos sobre cada vetor da matriz Gmatrix
            for(int m = 0; m < D; m++)
            {
                Eigen::VectorXi gm(D);
                Eigen::VectorXi Yvec = Eigen::VectorXi::Zero(D);
                num_nvec++;

                // Aqui estamos a ir buscar o vetor gm à matrix Gmatrix para depois usarmos na multiplicação
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
   
                    Yvec += gm(i) * subpacketAux;                    
                }

                // Adicionamos-lhe o Y1 que corresponde às mensagens de interferência
                Yvec += Y1;
                // Guardamos o vetor Y atual no array de vetores Y
                Y_vectors.push_back(Yvec);
            }
        }

        // Print Y's
        for(int i = 0; i < N; i++)
        {
            std::cout << "\nPrint Y" << i+1 << ": " << '\n';
            show_Vectorxi(Y_vectors[i],symbols_subpacket);
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

        for(int i = 1; i < N; i++)
        {
            Eigen::VectorXi Zn;

            Zn = Y_vectors[i] - Y_vectors[0];
            Z_vectors.push_back(Zn);
        }

        show_QuerysAndAnswersServer(server_indexs, n_vectors, Y_vectors, N, L, D, K);

        std::cout << '\n';
        for(int i = 0; i < N-1; i++)
        {
            std::cout << "\nPrint Z" << i+2 << ": " << '\n';
            show_Vectorxi(Z_vectors[i],Z_vectors[i].size()); 
        }


}




/*
Função para dado um determinado par, criar um i-sparse vector h de tamanho K-D com valores do finite field de ordem q,
e que tem i nonzero entries em posições random. 
Depos calculamos uma matrix G j-regular invertível (tem de ser quadrada e o determinante diferente de 0) de tamanho DxD em que 
G = [g1T,...,gDT] sobre o finite field de ordem q
Regras da matriz G são as seguintes:
    i) O vetor g1 possui exatamente j entradas não nulas em posições random
    ii) Para cada 2 <= m <= D, as posições da entradas diferentes de 0 no vetor gm são uma rotação circular das posições das entradas não nulas do vetor gm-1 (o anterior)
*/
void buildSparseVectorGAndVn (std::vector <int> &pairIJ, int D, int K, int q, int L, int N, int symbols_subpacket, std::mt19937 shuffle_random,std::vector<Message> &messages)
{
    // Obter o par i,j calculados na função anterior (build_MatrixM_fgAndChosePairij)
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

    // Dar shuffle no array de índices para depois usar de modo a colocar os elementos de forma aleatória
    std::shuffle(h_index.begin(), h_index.end(), shuffle_random);

    std::cout << "\nMostrar vetor de índices baralhado h: " << '\n';
    show_Vectorxi(h_index, K-D);

    // Preencher o i-sparse vector h com elementos do finite field de ordem q (mais especificamente nums_ToFill_I elementos == i_index)
    for(int j = 0; j < nums_ToFill_I; j++)
    {
        int num_finite_field = rand() % q;
        int num = num_finite_field != 0 ? num_finite_field : num_finite_field + 1;
        // Colocar os elementos de forma aleatória
        h(h_index(j)) = num;
    }

    std::cout << "\nMostrar sparse vector-i (h) com i entradas não nulas em posições aleatórias" << '\n';
    show_Vector(h,K-D);

    // Construir matriz G j-regular e invertível de tamanho DxD
    Eigen::MatrixXd Gmatrix(D,D);
    bool invertible = false;

    while(!invertible)
    {
        // Temos de construir o g1 vector de tamanho D com j entradas não nulas em posições aleatórias cada uma delas com valor do finite field de ordem q
        Eigen::VectorXd g1vec = Eigen::VectorXd::Zero(D);
        // Array com o mesmo tamanho com o acima que vai conter os índices e que vai servir para inserir os números gerados em posições aleatórias
        Eigen::VectorXi g1_index (D);
        // Valores auxiliares para guardarmo último valor e índice diferente de no vetor g1
        /*
        Vão funcionar da seguinte forma: 
        Primeiro g_index leva shuffle para as posições aleatórias;
        Com a execução do loop imaginemos que nums_ToFill_J = 2, 
        
        */

       // Preencher com os índices
        for(int i = 0; i < D; i++)
        {
            g1_index(i) = i;
        }  

        // Dar shuffle no array de índices para depois usar de modo a colocar os elementos de forma aleatória
        std::shuffle(g1_index.begin(),g1_index.end(),shuffle_random);

        std::cout << "\nMostrar vetor de índices baralhado g1: " << '\n';
        show_Vectorxi(g1_index, D);

        // Preencher o vector g1 pertencente á matrix G com elementos do finite field de ordem q (mais especificamente nums_ToFill_J elementos == j_index)
        for(int k = 0; k < nums_ToFill_J; k++)
        {
            int num_finitefield = rand() % q;
            int val = num_finitefield != 0 ? num_finitefield : num_finitefield + 1;
            // Colocar os elementos de forma aleatória
            g1vec(g1_index(k)) = val;
        }

        std::cout << "\nMostrar g1 vector: " << '\n';
        show_Vector(g1vec, D);

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
    
    // Chamada de função para construir os N vectors
    constructNVectors(D,K,q,i_index,L,N,symbols_subpacket,shuffle_random,messages,h,Gmatrix);
}




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
