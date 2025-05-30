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


// Função para mostrar vector std::vector <int>
void show_VectorSTD (std::vector <int> &vec, int size)
{
    std::cout << "\nVector STD: " << '\n';

    for(int i = 0; i < size; i++)
    {
        std::cout << vec[i] << " " << '\n';
    }

    std::cout << '\n';
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
    //show_MessagesSubpackets(messages,K,L,symbols_subpacket);
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

    // std::cout << "Probabilidade: " << random_num << '\n';

    // Percorremos a matriz de probabilidades 
    for(int i = 0; i < (K-D); i++)
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
            }
        }
    }

    return pairChosenValues;
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

    //show_Matrix(matrix_M,D,D);
    //show_Vector(colVec1s,D);


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

    // Probabilidades Pi,D onde 0 <= i <= K-D-1
    // Para cada i calculamos a probabilidade de Pi,n onde 1 <= n <=D
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
    }

    // show_Matrix(matrixProb,(K-D) + 1, D);
    // std::cout << "\nSoma das probabilidades de todos os pares: " << matrixProb.sum() << '\n'; 

    // Define vetor que conterá o par obtido pelas probabilidades e está disposto pela ordem i,j
    std::vector <int> pairChosenValues;
    
    pairChosenValues = chooseAleatoryPair(matrixProb, K, D, pairChosenValues);

    return pairChosenValues;    
}

/*
Função que recebe um vetor gD-1 da Matriz G  calcular o vetor gD e recebe duas variáveis que representam respetivamente o último valor maior que 0 do vetor g e o que vai armazenar o último novo valor
*/
std::vector <int> buildCircularVector (std::vector <int> &gDVec, int D, int lastVecValueInit, int &new_lastVecValue)
{
    // Vetor auxiliar que vai conter a rotação das posições
    std::vector <int> gDVecAux(D,0);

    // Percorremos o vetor fornecido e nas posições em que o valor atual for maior que 0 o que fazemos é no novo vetor na mesma posição 
    // que o anterior colocamos o último valor guardado que está em lastVecValueInit e atualizamos o último valor e assim sucessivamente  
    /*
        Imaginando que temos o seguinte vetor g1: 0 2 0 0 3 0 7 0
        O lastVecValueInit tem o valor 7, o que está a acontecer é: 
        Na 1º posição o valor é 0 não fazemos nada
        Na 2º posição o valor é > 0 (2) então o que fazemos é no novo vetor na mesma posição colocar o lastVecValueInit (7), e guardar no lastVecValueInit guardar o valor atual de gDVec(2),
        fazemos isto até ao fim do vetor e é suposto obtermos o seguinte vetor circular g2 através do anterior.
        g2: 0 7 0 0 2 0 3 0
        Para além disso guardamos o último valor do novo vetor (g2) na variável new_lastVecValue neste caso -> 3, para utlizarmos na próxima iteração recursiva da função fillGMatrixInRecursive no lugar de lastVecValue
    */
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


void fillGMatrixInRecursive (Eigen::MatrixXd &Gmatrix, std::vector <int> &g1Vec, int D, int lastVecValue, int line)
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
            Gmatrix(line,col) = g1Vec[col];
        }

        // Preencher o resto da Matriz recursivamente com os vetores circulares (gD) gerados a partir do anterior (gD-1)
        fillGMatrixInRecursive(Gmatrix, g1Vec, D, lastVecValue, line + 1);
    }

    // Se não for nenhum dos casos acima então vamos preencher o resto da matriz com o auxílio da função buildCircularVector descrita acima
    else 
    {   
        // Definimos uma variável que vai possuir o novo último valor (> 0) do vetor
        int new_lastVecValue;
        // Vetor circular gD obtido através do anterior gD-1
        std::vector <int> gDvecAux = buildCircularVector(g1Vec, D, lastVecValue, new_lastVecValue);

        // Preencher a próxima linha da matriz (line) com os valores do vetor gerado 
        for(int col = 0; col < D; col++)
        {
            Gmatrix(line,col) = gDvecAux[col];
        }

        // Preencher o resto da Matriz recursivamente ou devolvê-la preenchida
        fillGMatrixInRecursive(Gmatrix, gDvecAux, D, new_lastVecValue, line + 1);
    }
}


/*
Função para construir N (número de servers) vetores de tamanho K * L com valores aleatórios do finite field de ordem q
    Estes vetores são gerados um algoritmo random, baseado no conjunto de mensagens requeridas (W), as regras são as seguintes:
    - Se i = 0 (i -> i_index gerado de forma aleatória anteriormente) então v1 é um vetor de zeros, caso contrário v1 é o vetor coeficiente 
    correspondente a combinação linear Y1 definida como Y1 = h * [Xu1,1,...,XuK-D,1] ^T onde h é o sparse vector calculado anteriormente e u1,...,uK-D 
    são os índices das K-D mensagens de interferência numa ordem crescente ou decrescente mas fixa.
    - Para cada 1 <= l <= L e 1 <= m <= D o vetor v(l-1)*D+m+1 é o coeficente correspondente à combinação linear Y(l-1)*D+m+1 definida como:
    Y(l-1)*D+m+1 := Y1 + gm * [Xw,l,...XwD,l] ^T onde gm é um vetor da Gmatrix calculada anteriormente e w1,...,wD são os índices das D mensagens requeridas, numa ordem crescente ou decrescente mas fixa.
*/
void constructNVectors (int D, int K,int q, int i_index, int L, int N, int symbols_subpacket, std::mt19937 shuffle_random, std::vector<Message> &messages,std::vector<int> h,Eigen::MatrixXd Gmatrix)
{
    
    std::vector <Eigen::VectorXi> n_vectors (N, Eigen::VectorXi (K * L));
    std::vector<Eigen::VectorXi> interference_messages;
    std::vector<Eigen::VectorXi> demand_messages;
    Eigen::VectorXi Y1 = Eigen::VectorXi::Zero(symbols_subpacket);

    // Variável auxiliar com o valor de h sparse vector utilizada para facilitar nas operações
    Eigen::VectorXi hvec(K-D);
    // Copiar o valor de h sparse vector para a nova variável
    for(int i = 0; i < K-D; i++)
    {
        hvec(i) = h[i];
    }

   // Se i = 0 então v1 é preenchido apenas com zeros, caso contrário é dado por: Y1 = h * [Xu1,1,...,XuK-D,1] ^T
    if(i_index == 0)
    {
        for(int i = 0; i < K*L; i++)
        {
            n_vectors[0](i) = 0;
        }
    }

    else 
    {   
        // Construir um vetor de Eigen::VectorXi para as interference messages e para as demand messages cada uma vai 
        // conter os subpacotes correspondetes ao tipo de mensagem
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
                if(i >= D && j == 0)
                {
                    interference_messages.push_back(vecAux);
                }

                // Armazenar em demand messages os subpacotes de cada mensagem requerida
                else 
                {
                    demand_messages.push_back(vecAux);
                }
            }
        }

        // Calcular o vetor Y1 que vai servir para esconder as inteções do utilizador uma vez que é construído com os elementos 
        // das mensagens que não lhe interessam (mensagens de interferência)
        for (int k = 0; k < K-D; k++)
        {
            Y1 += hvec(k) * interference_messages[k]; 
        }
        
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
    // Obter o par i,j calculados na função anterior
    int i_index = pairIJ[0], j_index = pairIJ[1];
    // Construir o i-sparse vector h inicialmente preenchido com 0's
    std::vector <int> h(K-D, 0);
    // Array do mesmo tamanho que o acima que vai possuir os índices de h para os baralharmos
    std::vector <int> h_index(K-D);
    // Número de campos com valores do finite field de ordem q a preencher no i-sparse vector h e no vetor g1 respetivamente
    int nums_ToFill_I = i_index, nums_ToFill_J = j_index;

    // Preencher com os índices
    for(int i = 0; i < (K-D); i++)
    {
        h_index[i] = i;
    }

    // Dar shuffle no array de índices para depois usar de modo a colocar os elementos de forma aleatória
    std::shuffle(h_index.begin(), h_index.end(), shuffle_random);

    // Preencher o i-sparse vector h com elementos do finite field de ordem q (mais especificamente nums_ToFill_I elementos == i_index)
    for(int j = 0; j < nums_ToFill_I; j++)
    {
        int num_finite_field = rand() % q;
        int num = num_finite_field != 0 ? num_finite_field : num_finite_field + 1;
        // Colocar os elementos de forma aleatória
        h[h_index[j]] = num;
    }

    // Mostar sparse vector
    //show_VectorSTD(h,(K-D));

    // Comstruir matriz G j-regular e invertível de tamanho DxD
    Eigen::MatrixXd Gmatrix(D,D);
    bool invertible = false;

    while(!invertible) {
    // Temos de construir o g1 vector de tamanho D com j entradas não nulas em posições aleatórias cada uma delas com valor do finite field de ordem q
    std::vector <int> g1vec (D,0);
    // Array com o mesmo tamanho com o acima que vai conter os índices e que vai servir para inserir os números gerados em posições aleatórias
    std::vector <int> g1_index(D);
    // Valores auxiliares para guardarmo último valor e índice diferente de no vetor g1
    /*
    Vão funcionar da seguinte forma: 
    Primeiro g_index leva shuffle para as posições aleatórias;
    Com a execução do loop imaginemos que nums_ToFill_I = 3, então irá pegar nos primeiros três elementos de g_index, imaginemos agora que estes são 2,3,1
    nesta ordem então vamos ver se o primeiro elemento de g_index neste caso 2 é maior que o lastIndex iniciado a 0, como é então guardamos o valor val na variável lastVecValue
    e atualizamos o lastIndex que passará a ser 2 neste caso, na segunda iteração do loop vemos se g1_index[1] = 3 é maior ou igual que o lastIndex = 2 como é então atualizamos o
    lastVecValue com o número que guardaríamos em g1vec[3], e atualizamos o lastIndex = 3, na terceira e última iteração do loop temos g1_index[2] = 1 e vemos se g1_index[2] > lastIndex que é 3
    como não é então não atualizamos e temos que o lastIndex é 3 e o valor correspondente que ele tem. 
    Nota que isto é útil pois obtemos o maior índice onde o valor é diferente de 0 em g1vec, tendo os índices em g1vec ex: 0,1,2,3,4,.. rapidamente percebemos então que temos a última posição do vetor
    com um valor diferente de 0, útil para utilizarmos na rotação dos seguintes vetores.
    */   
    int lastIndex = 0, lastVecValue;
    // Preencher com os índices
    for(int i = 0; i < D; i++)
    {
        g1_index[i] = i;
    }

    // Dar shuffle no array de índices para depois usar de modo a colocar os elementos de forma aleatória
    std::shuffle(g1_index.begin(),g1_index.end(),shuffle_random);

    // Preencher o vector g1 pertencente á matrix G com elementos do finite field de ordem q (mais especificamente nums_ToFill_J elementos == j_index)
    for(int k = 0; k < nums_ToFill_J; k++)
    {
        int num_finitefield = rand() % q;
        int val = num_finitefield != 0 ? num_finitefield : num_finitefield + 1;
        // Colocar os elementos de forma aleatória
        g1vec[g1_index[k]] = val;

        // Se o índice atual for maior que o último
        // Então guardo o valor  associado a esse índice
        if(g1_index[k] >= lastIndex)
        {   
            lastVecValue = val;
        }

        // Atualizo o último indíce
        lastIndex = g1_index[k];    
    }

    // std::cout << "\nMostrar g1 vector " << '\n';
    // show_VectorSTD(g1vec,D);
    // std::cout << '\n';

    // Chamada de função para construir a matriz G através da rotação e de forma recursiva
    fillGMatrixInRecursive(Gmatrix, g1vec, D, lastVecValue,0);

    // Se o determinante da matriz G for maior que 0 então a matriz é invertível e podemos avançar caso contrário repetimos o processo de gerar a matriz G
    if(Gmatrix.determinant() != 0)
    {
        invertible = true;
    }
}

    // Chamda de função para construir os N vectors
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
    // std::cout << "Print i value: " << pairIJ[0] << '\n';
    // std::cout << "Print j value: " << pairIJ[1] << '\n';

    buildSparseVectorGAndVn(pairIJ, D, K, q, L, N, symbols_subpacket, shuffle_random,messages);    

    return 0;
    
}