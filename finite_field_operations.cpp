#include "finite_field_operations.h"

// Funções elementares (aritmética modular) em Finite Fields (Fq)

// Adição
int add_modFq (int a, int b)
{
    ZZ_p res = conv<ZZ_p>(a) + conv<ZZ_p>(b);

    return conv<int>(res);
}

// Subtração
int sub_modFq (int a, int b)
{
    ZZ_p res = conv<ZZ_p>(a) - conv<ZZ_p>(b);

    return conv<int>(res); 
}

// Multiplicação
int mult_modFq (int a, int b)
{
    ZZ_p res = conv<ZZ_p>(a) * conv<ZZ_p>(b);

    return conv<int>(res);
}

// Adição em Fq sobre vetores Eigen usando como auxiliar as funções acima 
Eigen::VectorXi addVectorsFq (const Eigen::VectorXi &vec1, const Eigen::VectorXi &vec2)
{
    Eigen::VectorXi vec_res (vec2.size());

    for(int i = 0; i < vec1.size(); i++)
    {
        vec_res(i) = add_modFq(vec1(i),vec2(i));
    }

    return vec_res;
}

// Multiplicação em Fq sobre vetores Eigen usando como auxiliar as funções acima 
Eigen::VectorXi multVectorXVal (const Eigen::VectorXi &vec1, int val)
{
    Eigen::VectorXi vec_res(vec1.size());

    for(int i = 0; i < vec1.size(); i++)
    {
        vec_res(i) = mult_modFq(vec1(i), val);
    }

    return vec_res;
}
