#pragma once

/*
    Operações Elementares em Finite Fields de ordem q (Fq)
    a,b -> são elementos do finite field de ordem q e q é a ordem do mesmo (neste caso >= 3 e tem de ser primo)
*/

// Adição em aritmética modular para um dado Finite Field de ordem q
int add_modFq (int a, int b, int q);

// Subtração em aritmética modular para um dado Finite Field de ordem q
int sub_modFq (int a, int b, int q);

// Multiplicação em aritmética modular para um dado Finite Field de ordem q
int mult_modFq (int a, int b, int q);