#pragma once

#include <vector>
#include "finite_field_operations.h"

 
// Estrutura que representa um subpacote.
//   Um subpacote contém um conjunto de números pertencentes ao finite field de ordem q.
//   Então ele pode ser representado como um vetor de elementos ZZ_p.
struct Subpacket
{
    std::vector <ZZ_p> numsR_finite_field;
};



// Estrutura que representa uma mensagem.
//    Uma mensagem é composta por um conjunto de subpacotes, em que cada subpacote tem a estrutura descrita acima (Subpacket).
//    - L -> representa o grau de supacketização (número de subpacotes em que é dividida cada mensagem)
//    - m -> tamanho de cada mensagem
//
// Cada subpacote contém m/L símbolos = symbols_subpacket
struct Message
{
    std::vector <Subpacket> subpackets;
};
