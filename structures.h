#pragma once

#include <vector>
#include "finite_field_operations.h"

/* 
Definir uma estrutura para subpacotes;
   Um subpacote é um elemento que contém um conjunto de números pertencentes ao finite field de ordem q
   Então ele pode ser representado como um vetor de inteiros
*/
struct Subpacket
{
    std::vector <ZZ_p> numsR_finite_field;
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
