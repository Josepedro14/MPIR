
/*
- Um servidor remoto 
- Uma base de dados (composta por um número par de mensagens binárias de igual comprimento representadas por: X1,...,Xk)
- O utilizador deseja descarregar uma mensagem do servidor (Objetivo: não revelar a mensagem desejada ao servidor)
- O utilizador possui uma mensagem Xs (1 <= s <= K) como informação lateral escolhida aleatoriamente entre as outras K mensagens e desconhecida pelo servidor
- Seja Xw a mensagem que o utilizador deseja e Xs a mensagem que possui com W,S ∈ {1,...,K} e W ≠ S. 

- Toy Example --> [8] Private Information Retrieval With Side Information The Single Server Case 
( Exemplo 1: Esquema PIR de Partição e Codificação )
*/

#include <iostream>
#include <vector>
#include <bitset>

// Conjuntos de tamanho 2 da partição aleatória
struct Set
{
    int values [2];                                      
};                                                            

// Estrutura das mensagens presentes na base de dados
struct Message
{
    int value;
    bool posChanged;
};

std::vector<Message> messages= {{0b00000001 , false}, {0b00000010, false}, {0b00000011 , false}, {0b00000100 , false}, {0b00000101, false}, {0b00000110, false}, {0b00000111, false}, {0b00001000, false}};

// Função para aplicar o XOR nos elementos dos conjuntos na resposta do servidor
int applyXOR (int message1, int message2)
{
    return message1 ^ message2;
}

void build_Particao (int K,int W,int S, Set particao_aleatoria[])
{
    // Retirar a mensagem desejada e a informação lateral do conjunto de mensagens para formar um par, e eliminá-las do conjunto de mensagens
    int desired_message = messages[W - 1].value;
    int known_message = messages[S - 1].value;

    if(W > S)
    {
        messages.erase(messages.begin() + (W - 1));
        messages.erase(messages.begin() + (S - 1));
    }

    else 
    {
        messages.erase(messages.begin() + (S - 1));
        messages.erase(messages.begin() + (W - 1));
    }
    

    int n_size = K - 2;
    // Inicializar um novo array sem as duas mensagens anteriormente retiradas
    std::vector<int> aleatory_messages(n_size);

    // Inicializar o gerador de números aleatórios sendo a semente o tempo atual em segundos
    srand(time(NULL));

    int num = 0;

    // Dispor aleatoriamente os elementos das mensagens no novo array criado para depois formarmos os conjuntos
    while(n_size > 0)
    {
        int random_index = rand() % (K - 2);

        //std::cout << "Random index " << random_index << '\n';

        if(!messages[random_index].posChanged)
        {
            aleatory_messages[random_index] = messages[num].value;

            //std::cout << messages[random_index].value << " to pos " << random_index << '\n';
            
            messages[random_index].posChanged = true;
            num += 1;
            n_size -= 1;
            
        }
    }


    // Formar os conjuntos
    int index = 0;
    int random_set_index = rand() % (K/2);

    for(int j = 0; j < K/2; j++)
    {
        if(j == random_set_index)
        {
            particao_aleatoria[j].values[0] = desired_message;
            particao_aleatoria[j].values[1] = known_message;
        }

        else 
        {
            particao_aleatoria[j].values[0] = aleatory_messages[(2 * index)];
            particao_aleatoria[j].values[1] = aleatory_messages[(2 * index) + 1];
            index += 1;
        }
    }

    
    // Ver conjuntos Criados
    
    /*

    for (int k = 0; k < K/2; k++)
    {
        std::cout  << k + 1 << "  -->  " << particao_aleatoria[k].values[0] << "," << particao_aleatoria[k].values[1] << "\n";
    }
        
    */
    


}

void serverResponse (int K, Set particao_aleatoria[])
{
    for(int i = 0; i < K/2; i++)
    {
        int xOR_set_res = applyXOR(particao_aleatoria[i].values[0], particao_aleatoria[i].values[1]);
        std::cout << "Server response " << i + 1 << ": " << std::bitset<8>(xOR_set_res) << '\n';
    }
}

int main ()
{
    /*
        K --> Representa o número de mensagens na base de dados
        W --> Representa a mensagem desejada pelo utilizador
        S --> Representa a mensagem (informação lateral) que o utilizador já possui e o servidor desconhece
    */
    int K = 8, W = 2, S = 4;
    Set particao_aleatoria [(K/2)];

    build_Particao(K,W,S,particao_aleatoria);
    
    serverResponse(K,particao_aleatoria);

    // Para ver o conteúdo da mensagem desejada basta fazermos XOR das respostas com a mensagem que já possuíamos e ver qual resultado é a mensagem desejada.


    return 0;
}

