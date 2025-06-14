#include "matrix_operations.h"

// Função que multiplica uma matrix (matrix_M) de tamanho (D,D) por ela mesma (num_of_times) vezes
// Retorna a matriz: (matrix_M) ^ (num_of_times)
Eigen::MatrixXd multiplyMatrixNTimes (Eigen::MatrixXd &matrix_M, int num_of_times, int D)
{
    // Inicializamos a matriz (matrixRes) como matriz identidade (D,D)
    // Multiplicar uma matriz pela identidade resulta na própria matriz
    Eigen::MatrixXd matrixRes = Eigen::MatrixXd::Identity(D,D);

    for(int i = 0; i < num_of_times; i++)
    {
        matrixRes *= matrix_M;
    }

    return matrixRes;
}


// Função que realiza uma rotação circular à direita nos elementos do vetor de tamanho D fornecido como parâmetro.
// Exemplo: Imaginemos que temos um valor 2 na última posição do vetor então no novo vetor este elemento 2 estará na primeira posição.
Eigen::VectorXd buildCircularVector (Eigen::VectorXd &gDVec, int D)
{
    // Vetor auxiliar inicializado com zeros e que vai conter a rotação à direita das posições do vetor anterior
    Eigen::VectorXd gDVecAux = Eigen::VectorXd::Zero(D);

    for(int i = 0; i < D; i++)
    {   
        // Calculamos a nova posição para o elemento gDvec(i) no novo vetor gDVecAux
        int position = i + 1;

        // Se a nova posição calculada for maior que o tamanho do vetor então a nova posição é 0 (a primeira)
        if(position > D-1)
        {
            position = 0;
        }

        // Se o valor for diferente de zero aplicamos a rotação
        if(gDVec(i) > 0)
        {
            gDVecAux(position) = gDVec(i);
        }
        
    }

    // Retorna o vetor com a rotação
    return gDVecAux;
}


// Função que preenche uma matriz quadrada Gmatrix (D,D) de forma recursiva com vetores gerados por rotações circulares à direita dos anteriores.
// A primeira linha é preenchida com o vetor (g1Vec) passado como parâmetro.
void fillGMatrixInRecursive (Eigen::MatrixXd &Gmatrix, Eigen::VectorXd &g1Vec, int D, int line)
{   
    // Se a linha atual for maior que o tamanho da Matriz então retorna a Matriz obtida
    if(line >= D)
    {
        return;
    }

    // Se a linha for a primeira vamos usar o vetor (g1Vec) com j entradas não nulas para preencher a primeira linha da Matriz
    if(line == 0)
    {
        for(int col = 0; col < D; col++)
        {
            Gmatrix(line,col) = g1Vec(col);
        }

        // Continuar a preencher a Matriz recursivamente aplicando rotações circulares à direita sobre o vetor anterior.
        fillGMatrixInRecursive(Gmatrix, g1Vec, D, line + 1);
    }

    // Se não for nenhum dos casos acima então vamos preencher o resto da matriz com o auxílio da função acima (buildCircularVector) 
    else 
    {   
        // Vetor circular gDvecAux obtido através do anterior g1Vec
        Eigen::VectorXd gDvecAux = buildCircularVector(g1Vec, D);

        // Preencher a próxima linha da matriz (line) com os valores do vetor gerado 
        for(int col = 0; col < D; col++)
        {
            Gmatrix(line,col) = gDvecAux(col);
        }

        // Preencher o resto da Matriz recursivamente ou devolvê-la preenchida
        // No caso de existir mais uma linha a preencher na Matriz o vetor a ser usado agora (gDvecAux) é o que sofreu rotação a partir do original (g1Vec)
        fillGMatrixInRecursive(Gmatrix, gDvecAux, D, line + 1);
    }
}