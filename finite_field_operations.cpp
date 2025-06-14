#include "finite_field_operations.h"

// Funções elementares (aritmética modular) em Finite Fields (Fq), com q já inicializado no main


// Função que realiza a adição de dois inteiros (a,b) em Fq
// Retorna o resultado da adição em Fq como inteiro 
int add_modFq (int a, int b)
{   
    // Converter do tipo int para elementos de Fq (tipo ZZ_p) e realizar a adição 
    ZZ_p res = conv<ZZ_p>(a) + conv<ZZ_p>(b);

    // Converter o resultado para int
    return conv<int>(res);
}


// Função que realiza a subtração de dois inteiros (a,b) em Fq
// Retorna o resultado da subtração em Fq como inteiro 
int sub_modFq (int a, int b)
{   
    // Converter do tipo int para elementos de Fq (tipo ZZ_p) e realizar a subtração
    ZZ_p res = conv<ZZ_p>(a) - conv<ZZ_p>(b);

    // Converter o resultado para int
    return conv<int>(res); 
}


// Função que realiza a multiplicação de dois inteiros (a,b) em Fq
// Retorna o resultado da multiplicação em Fq como inteiro 
int mult_modFq (int a, int b)
{
    // Converter do tipo int para elementos de Fq (tipo ZZ_p) e realizar a multiplicação
    ZZ_p res = conv<ZZ_p>(a) * conv<ZZ_p>(b);

    // Converter o resultado para int
    return conv<int>(res);
}


// Função que faz a adição de dois vetores (Eigen::VectorXi) sobre Fq
// Cada elemento do novo vetor (vec_res) é obtido através da adição sobre Fq e de dois elementos um de vec1 e um de vec2
// Esta função é aplicada em vetores com o mesmo tamanho (vec1.size() == vec2.size())
Eigen::VectorXi addVectorsFq (const Eigen::VectorXi &vec1, const Eigen::VectorXi &vec2)
{
    Eigen::VectorXi vec_res (vec2.size());

    // Aplicamos a função definida acima (add_modFq)
    for(int i = 0; i < vec1.size(); i++)
    {
        vec_res(i) = add_modFq(vec1(i),vec2(i));
    }

    return vec_res;
}


// Função que faz a multiplicação de um escalar por um vetor (Eigen::VectorXi) sobre Fq (o escalar pertence a Fq)
// Cada elemento do novo vetor (vec_res) é obtido através da multiplicação sobre Fq entre o escalar e o elemento vec1(i) 
Eigen::VectorXi multVectorXVal (const Eigen::VectorXi &vec1, int val)
{
    Eigen::VectorXi vec_res(vec1.size());

    // Aplicamos a função definida acima (mult_modFq) 
    for(int i = 0; i < vec1.size(); i++)
    {
        vec_res(i) = mult_modFq(vec1(i), val);
    }

    return vec_res;
}
