#include "finite_field_operations.h"

// Funções elementares (aritmética modular) em Finite Fields (Fq)

// Adição
int add_modFq (int a, int b, int q)
{
    return (a + b) % q;
}

// Subtração
int sub_modFq (int a, int b, int q)
{
    return (a - b + q) % q; 
}

// Multiplicação
int mult_modFq (int a, int b, int q)
{

    return (a * b) % q;
}