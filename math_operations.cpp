#include "math_operations.h"


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



int calculateCombinations (int D, int j)
{
    if(j == 0 || j == D)
    {
        return 1;
    }

    return fatorial(D) / (fatorial(j) * fatorial(D-j));
}



double calculateBJForM (int j, int D, int L)
{
    return (double)(D * L) / (double)calculateCombinations(D,j);
}


// Função que calcula o cofator para a coluna atual (line é fixo pois o teorema de laplace está a ser aplicado sempre à linha 0)
int calculateCofator (int line, int col)
{
    return ((line + col) % 2 == 0) ? 1 : -1;

}


// Função recursiva que calcula o determinante de uma matriz (D,D)
// Utiliza o Teorema de Laplace e é sempre aplicado à primeira linha da matriz/submatriz
int calculateDetrec (Eigen::MatrixXi Gmatrix, int D)
{   
    // Se o tamanho da matriz/submatriz for 1 então o determinante é o próprio valor
    if(D == 1)
    {
        return Gmatrix(0,0);
    }

    // Variável que vai conter o valor do determinante 
    int determinante = 0;
    // Submatriz de tamanho (D-1,D-1) e que contém todos os elementos da matriz inicial, exceto os da coluna atual (col) e linha 0 (primeira).
    Eigen::MatrixXi subMatrix(D-1,D-1);

    // Iteramos pelos elementos da primeira linha da matriz
    for(int col = 0; col < D; col++)
    {
        //  Para o elemento atual se for diferente de 0 então calculamos o determinante, caso contrário não vale a pena pois segundo o teorema de Laplace
        //  o determinante correspondente a esse valor seria também 0
        if(Gmatrix(0,col) != 0)
        {
            // variável que contém a posição da linha atual da submatriz
            int poslinesub = 0;
            // Vamos percorrer todos os elementos da matriz Gmatrix
            for(int linesub = 0; linesub < D; linesub++)
            {   
                if(linesub != 0)
                {   
                    // variável que contém a posição da coluna atual da submatriz
                    int poscolsub = 0;
                    for(int colsub = 0; colsub < D; colsub++)
                    {
                        // Se a linha atual não for a primeira que é onde está a ser aplicado o teorema de laplace 
                        // Se a coluna atual também não coincidir com a coluna do elemento da primeira linha escolhido (col)
                        if(colsub != col)
                        {   
                            // Colocar o elemento na submatriz
                            subMatrix(poslinesub,poscolsub) = Gmatrix(linesub,colsub);
                            poscolsub++;
                        }
                    }
                    poslinesub++;
                }
            }

        // Calcular o cofator atual
        int cofator = calculateCofator(0,col);
        // Chamar recursivamente para calcular o determinante da submatriz
        int determinateSubMatrix = calculateDetrec(subMatrix, D-1);
        // Calcular o valor do determinante segundo as operações do finite field de ordem q
        determinante = add_modFq(determinante,mult_modFq(cofator,mult_modFq(Gmatrix(0,col),determinateSubMatrix)));

        }
    }

    return determinante;
}

// Função que resolve um sistema de equações lineares utilizando o método de eliminação de Gauss
void calculateSubpacketsGauss (Eigen::MatrixXi A, int rows, int cols, int L, int D, int K, int symbols_subpacket)
{   
    std::cout << "\nPrint matriz inicio: " << '\n';
    show_Matrix_xi(A,rows,cols);

    for(int i = 0; i < cols - 1; i++)
    {
        if(A(i,i) == 0)
        {            
            for(int j = i + 1; j < rows; j++)
            {
                if(A(j,i) != 0)
                {
                    A.row(i).swap(A.row(j));
                    std::cout << "\nPrint A pos Swap: " << '\n';
                    show_Matrix_xi(A,rows,cols);
                    break;
                }

            }
        }


        int pivo = A(i,i);

        for(int j = i+1; j < rows; j++)
        {
            if(A(j,i) != 0)
            {
                int inverse_pivoFq = inverso_modFq(pivo);
                std::cout << "\nPrint inverso: " << inverse_pivoFq << '\n';
                int scalar = mult_modFq(A(j,i), inverse_pivoFq);
                std::cout << "\nPrint scalar: " << scalar << '\n';

                for(int k = 0; k < cols; k++)
                {
                    A(j,k) = sub_modFq(A(j,k), mult_modFq(scalar, A(i,k)));
                }

                std::cout << "\nPrint matrix pos mult escalar: " << '\n';
                show_Matrix_xi(A,rows,cols);
            }
        }
    }
 
    std::cout << "\nPrint matrix fim: " << '\n';
    show_Matrix_xi(A,rows,cols);


    std::vector<int> solution(cols - 1, 0); 

    for (int i = cols - 2; i >= 0; i--) 
    {
        int soma = 0;

        for (int j = i + 1; j < cols - 1; j++)
        {
            soma = add_modFq(soma, mult_modFq(A(i, j), solution[j]));
        }

        int rhs = sub_modFq(A(i, cols - 1), soma); 
        int pivot = A(i, i);

        if (pivot != 0)
        {
            int inv = inverso_modFq(pivot);
            solution[i] = mult_modFq(rhs, inv);
        }
 
    }

    std::cout << "\nSoluções encontradas:\n";
    for (size_t i = 0; i < solution.size(); ++i)
    {
        std::cout << "s" << i + 1 << " = " << solution[i] << std::endl;
    }

    
}