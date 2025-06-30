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



// Função que resolve um sistema de equações lineares utilizando o método de eliminação de Gauss, transformando a matriz recebida como parâmetro na sua forma escalonada 
void calculateMatrixGauss (Eigen::MatrixXi A, int rows, int cols, int &det, bool sol)
{   
    // Iterar pelas linhas da matriz
    for(int i = 0; i < cols - 1; i++)
    {
        // Se o pivo da linha atual for 0
        if(A(i,i) == 0)
        {           
            // Percorremos a coluna para ver se existe algum elemento diferente de 0 
            for(int j = i + 1; j < rows; j++)
            {
                if(A(j,i) != 0)
                {
                    // Se existir trocamos a linha do pivo com a linha do elemento encontrado
                    A.row(i).swap(A.row(j));
                    det = mult_modFq(det,-1);
                    break;
                }

            }
        }

        // Obter o valor do pivo atual
        int pivo = A(i,i);

        // Se encontrarmos algum elemento na coluna do pivo diferente de 0, vamos eliminá-lo
        for(int j = i+1; j < rows; j++)
        {
            if(A(j,i) != 0)
            {
                // Calcular o valor do inverso para o pivo atual
                int inverse_pivoFq = inverso_modFq(pivo);
                // O escalar representa o valor da divisão entre o elemento encontrado e o pivo, como a divisão não está definida multiplica-se pelo inverso do pivo
                int scalar = mult_modFq(A(j,i), inverse_pivoFq);

                // Eliminamos esse elemento através da subtração da sua linha com a linha do pivo multiplicada pelo valor escalar
                for(int k = 0; k < cols; k++)
                {
                    A(j,k) = sub_modFq(A(j,k), mult_modFq(scalar, A(i,k)));
                }
            }
        }
    }


    // Se a matriz for quadrada calcular o determinante
    if(cols == rows)
    {
        // Como a matriz já está na forma escalonada para obter o determinante basta multiplicar os elementos da diagonal principal
        for(int i = 0; i < rows; i++)
        {
            for(int j = 0; j < cols; j++)
            {
                if(i == j)
                {
                    det = mult_modFq(A(i,j),det);
                }
            }
        }

        std::cout << "\nGmatrix Escalonada " << '\n';
        std::cout << A << '\n';
    }


    if(sol)
    {
        // Vetor que vai conter as soluções do sistema (valores dos subpacotes)
        std::vector<int> solutions(cols - 1, 0); 

        // Percorremos todas as linhas da matriz a partir da última
        for (int line = cols - 2; line >= 0; line--) 
        {
            int soma = 0;

            // Para cada coluna a seguir ao pivo (da linha atual), multiplicar o elemento da coluna pela respetiva solução da coluna e adicionar à variável soma
            for (int col = line + 1; col < cols - 1; col++)
            {
                soma = add_modFq(soma, mult_modFq(A(line, col), solutions[col]));
            }

            // Obter o valor do pivo da linha atual
            int pivo = A(line, line);
            // Subtrair ao valor da última coluna (valor do vetor Z_n) a soma obtida na operação anterior
            int valZn = sub_modFq(A(line, cols - 1), soma); 

            if (pivo != 0)
            {
                // Calcular o inverso do pivo, para fazer a multiplicação (a divisão não é aplicada então fazemos a multiplicação pelo inverso)
                int inv = inverso_modFq(pivo);
                //  multiplicar o inverso pelo valor de Z_n na última coluna e guardar no vetor solutions para a coluna respetiva
                solutions[line] = mult_modFq(valZn, inv);
            }
 
        }

        // Dar print às soluções dos elementos dos subpacotes encontrados
        std::cout << "\nSoluções encontradas:\n";
        for (size_t i = 0; i < solutions.size(); ++i)
        {
            std::cout << "s" << i + 1 << " = " << solutions[i] << std::endl;
        }
    }
}

