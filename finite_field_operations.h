#pragma once

#include <eigen3/Eigen/Dense>
#include <NTL/ZZ_p.h>
#include <NTL/ZZ.h>

using namespace NTL;

/*
    Operações Elementares em Finite Fields de ordem q (Fq)
    a,b -> são elementos do finite field de ordem q e q é a ordem do mesmo (neste caso >= 3 e tem de ser primo)
*/

// Adição em aritmética modular para um dado Finite Field de ordem q
int add_modFq (int a, int b);

// Subtração em aritmética modular para um dado Finite Field de ordem q
int sub_modFq (int a, int b);

// Multiplicação em aritmética modular para um dado Finite Field de ordem q
int mult_modFq (int a, int b);

// Adição em Fq para cada elemento pertencente aos vetores Eigen apresentados 
Eigen::VectorXi addVectorsFq (const Eigen::VectorXi &vec1, const Eigen::VectorXi &vec2);

// Multiplicação em Fq por um valor para cada elemento pertencente aos vetores Eigen apresentados 
Eigen::VectorXi multVectorXVal (const Eigen::VectorXi &vec1, int val);

