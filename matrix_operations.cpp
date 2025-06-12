#include "matrix_operations.h"


Eigen::MatrixXd multiplyMatrixNTimes (Eigen::MatrixXd &matrix_M, int num_of_times, int D)
{
    Eigen::MatrixXd matrixRes = Eigen::MatrixXd::Identity(D,D);

    for(int i = 0; i < num_of_times; i++)
    {
        matrixRes *= matrix_M;
    }

    return matrixRes;
}



Eigen::VectorXd buildCircularVector (Eigen::VectorXd &gDVec, int D)
{
    // Vetor auxiliar que vai conter a rotação das posições
    Eigen::VectorXd gDVecAux = Eigen::VectorXd::Zero(D);

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