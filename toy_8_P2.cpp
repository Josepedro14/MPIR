    #include <iostream>
    #include <vector>
    #include <math.h>
    #include <random>
    #include <NTL/ZZ.h>
    #include <algorithm>
    #include <NTL/ZZ_p.h>

    using namespace NTL;

    // Variáveis globais 
    int idx, remainingSetSizes;
    int indexC;

    // Estrutura que representa um conjunto e contém:
    //  Índices das mensagens da base de dados 
    //  Valores de probabilidades do conjunto
    //  Informação sobre a presença do índice da mensagem a recuperar
    struct Set 
    {
        std::vector<int> values;
        double begin_Prob;
        double end_Prob;
        bool hasindexW = false;
    };


    // Função que dá print das probabilidades de cada conjunto e o índice do conjunto escolhido aleatoriamente
    void mostrarProbabilidadesEEscolhido(std::vector<Set> &aleatory_Partition, int g, int chosen_index) {
        
        for (int j = 0; j < g; j++) {
            std::cout << "Conjunto " << j + 1 << " -> Probabilidade: [" 
                    << aleatory_Partition[j].begin_Prob << ", " 
                    << aleatory_Partition[j].end_Prob << "]\n";
        }
        
        std::cout << "\nÍndice do conjunto escolhido: " << chosen_index << "\n";
    }


    //  Função que dá print dos elementos de um conjunto
    void print_Sets(std::vector<Set> &aleatory_Partition, int g) {

        for (int i = 0; i < g; i++) {
            
            std::cout << "Partition " << i + 1 << ":\n";
            std::cout << "  Values: ";
            
            for (size_t j = 0; j < aleatory_Partition[i].values.size(); j++) 
            {
                std::cout << aleatory_Partition[i].values[j] << " ";
            }
            std::cout << "\n----------------------\n";
        }
    }


    // Função utilizada para mostrar as mensagens da base de dados
    void print_Messages (std::vector<int> message)
    {
        for(size_t i = 0; i < message.size(); i++)
        {
            std::cout << message[i] << " ";
        }

        std::cout << '\n';
    }


    // Função para aleatorizar a posição dos índices das mensagens da BD num vetor
    void shuffleMessages(std::vector<int> &messages_index, std::mt19937 shuffle_random)
    {

        std::shuffle(messages_index.begin(),messages_index.end(), shuffle_random);
    }

    
    // Função para retornar o índice de um conjunto que vai ser escolhido aleatoriamente 
    int chooseRandomNumPlusIndex (std::vector<Set> aleatory_Partition, int g, ZZ maxVal)
    {   
        // Gerar um número aleatório entre 0 e maxVal [0,maxVal)
        ZZ random = RandomBnd(maxVal);
        // Obtemos um valor aleatório entre [0,1) através da divisão do valor gerado com o passado em parâmetro
        double random_num = conv<double>(random) / conv<double>(maxVal);

        //std::cout << "Probabilidade: " << random_num << '\n';

        for(int i = 0; i <= g-1; i++)
        {
            if(random_num >= aleatory_Partition[i].begin_Prob && random_num <= aleatory_Partition[i].end_Prob)
            {
                return i;
            }
        }

        return 0;
    }


    // Função que calcula as respostas do servidor para os g conjuntos
    void serverAnswers(std::vector<Set> aleatory_Partition, int g, std::vector<int> &answers, std::vector<int> database_messages)
    {
        // Cada resposta corresponde a um produto interno da posição dos índices do conjunto g pelas respetivas mensagens da base de dados
        for(int i = 0; i < g; i++)
        {
            for(size_t j = 0; j < aleatory_Partition[i].values.size(); j++)
            {
                answers[i] += database_messages[aleatory_Partition[i].values[j]]; 
            }
        }

    }


    // Função que constrói os conjuntos/partições segundo um conjunto de regras
    void build_Partition (int K, int M, int W, int g, bool isSpecialCase, std::vector<int> S, std::vector<Set> &aleatory_Partition, std::mt19937 shuffle_random, ZZ maxVal)
    {
        // vetor que contém os índices das mensagens presentes na BD
        std::vector<int> messages_index;

        for(int i = 0; i < K; i++)
        {
            messages_index.push_back(i);
        }

        // Caso especial então o conjunto P1 = W U S, e o resto dos g-1 conjuntos são preenchidos de forma aleatória
        if(isSpecialCase)
        {
            // campo do conjunto usado para controlar a presença do índice W
            aleatory_Partition[0].hasindexW = true;

            // Adicionamos o índice W e removemo-lo do conjunto de índices
            aleatory_Partition[0].values.push_back(W);
            messages_index.erase(std::remove(messages_index.begin(),messages_index.end(),W), messages_index.end());

            // Preenchemos o resto do conjunto P1 com os índices de side information S
            for(size_t i = 0; i < S.size(); i++)
            {
                aleatory_Partition[0].values.push_back(S[i]);
                messages_index.erase(std::remove(messages_index.begin(),messages_index.end(),S[i]), messages_index.end());
            }

            // baralhar o resto dos índices nas mensagens
            shuffleMessages(messages_index,shuffle_random);

            // variável de controlo do número de mensagens restantes
            remainingSetSizes = K - (1 + S.size());
            idx = 0;


            for(int i = 1; i <= g-1; i++)
            {
                // Tamanho do conjunto i atual
                int size_Set = std::min(remainingSetSizes,(M+1));
                remainingSetSizes -= size_Set;

                // Preencher o conjunto com índices aleatórios das mensagens
                for(int j = 0; j < size_Set; j++)
                {
                    aleatory_Partition[i].values.push_back(messages_index[idx]);
                    idx++;
                }

            }
        }

        else 
        {
            // Caso geral em que os conjuntos P1,...,Pg têm tamanho (M+1) e Pg tem tamanho K-(g-1)(M+1)
            // Probabilidade: P1,...,Pg-1 --> (M+1)/ K
            // Probabilidade: Pg --> (K - (g-1)(M+1)) / K

            // Definir as probabilidades para o primeiro conjunto
            aleatory_Partition[0].begin_Prob = 0.0;
            aleatory_Partition[0].end_Prob = (double) (M+1) / K;

            // Obter as probabilidades dos conjuntos usando a probabilidade já calculada no primeiro
            for(int j = 1; j <= g-1; j++)
            {   
                if(j != (g-1))
                {
                    aleatory_Partition[j].begin_Prob = aleatory_Partition[j-1].end_Prob;
                    aleatory_Partition[j].end_Prob = aleatory_Partition[j].begin_Prob +  (double) (M+1) / K;
                }

                // Definir a probabilidade do último conjunto
                else 
                {
                    aleatory_Partition[g-1].begin_Prob =  1.0 - (double) (K - (g-1) * (M+1)) / K;
                    aleatory_Partition[g-1].end_Prob = 1.0;
                }
            }

            // O utilizador escolhe aleatoriamente um desses conjuntos com base nas probabilidades de um número aleatório gerado entre [0,1) (através da função chooseRandomlyPlusIndex)
            int chosen_index = chooseRandomNumPlusIndex(aleatory_Partition, g, maxVal);
            // O indice W esta no conjunto chosen_index
            aleatory_Partition[chosen_index].hasindexW = true;
            indexC = chosen_index + 1;

            // Se o conjunto escolhido for Pg então o utilizador coloca o índice W dentro desse conjunto e preenche os restantes K - (g-1)(M+1) - 1 espaços com os índices da side information
            // O resto dos conjuntos são preenchidos com os restantes índices de forma aleatória
            if(chosen_index == g-1)
            {
                // Adicionamos o índice W e removemo-lo do conjunto de índices
                aleatory_Partition[g-1].values.push_back(W);
                messages_index.erase(std::remove(messages_index.begin(),messages_index.end(),W), messages_index.end());

                // Preenchemos as restantes posições desse conjunto com índices de side information escolhidos de forma aleatória
                // Após serem adicionados são removidos do conjunto de side information S, e do conjunto de índices
                while(aleatory_Partition[g-1].values.size() != (size_t) (K - (g-1) * (M+1)))
                {
                    int random_indexS = conv<int>(RandomBnd((ZZ) S.size()));
                    aleatory_Partition[g-1].values.push_back(S[random_indexS]);
                    messages_index.erase(std::remove(messages_index.begin(), messages_index.end(), S[random_indexS]), messages_index.end());
                    S.erase(std::remove(S.begin(), S.end(), S[random_indexS]), S.end());
                }

                // baralhar o resto dos índices nas mensagens
                shuffleMessages(messages_index,shuffle_random);

                // variável de controlo do número de mensagens restantes
                remainingSetSizes = K - aleatory_Partition[g-1].values.size();
                idx = 0;


                for(int i = 0; i <= g-2; i++)
                {
                    // Tamanho do conjunto i atual
                    int size_Set = std::min(remainingSetSizes, (M+1));
                    remainingSetSizes -= size_Set;

                    // Preencher os restantes conjuntos com o conjunto de índices restantes
                    for(int j = 0; j < size_Set; j++)
                    {
                        aleatory_Partition[i].values.push_back(messages_index[idx]);
                        idx++;
                    }
                }
                
            }

            // Se o conjunto escolhido for um P ∈ P1,...,Pg-1, então o utilizador coloca o índice da mensagem desejada W e os índices de side information nesse conjunto.
            // O resto dos conjuntos são preenchidos com os restantes índices de forma aleatória
            else 
            {
                // Adicionamos o índice W e removemo-lo do conjunto de índices
                aleatory_Partition[chosen_index].values.push_back(W);
                messages_index.erase(std::remove(messages_index.begin(),messages_index.end(),W), messages_index.end());

                // Preenchemos o resto do conjunto P com os índices de side information S e removemo-los do vetor de índices
                for(size_t i = 0; i < S.size(); i++)
                {
                    aleatory_Partition[chosen_index].values.push_back(S[i]);
                    messages_index.erase(std::remove(messages_index.begin(),messages_index.end(),S[i]), messages_index.end());
                }

                // baralhar o resto dos índices nas mensagens
                shuffleMessages(messages_index,shuffle_random);

                // variável de controlo do número de mensagens restantes
                remainingSetSizes = (int) messages_index.size();
                idx = 0;

                // Preencher os restantes conjuntos que não sejam o escolhido com os índices restantes
                for(int i = 0; i <= g-1; i++)
                {
                    if(i != chosen_index)
                    {   
                        // Tamanho do conjunto i atual
                        int size_Set = std::min(remainingSetSizes, (M+1));
                        remainingSetSizes -= size_Set;

                        for(int j = 0; j < size_Set; j++)
                        {
                            aleatory_Partition[i].values.push_back(messages_index[idx]);
                            idx++;
                        }
                    }
                }

            }
        }

        //mostrarProbabilidadesEEscolhido(aleatory_Partition,g, indexC);

        // O utilizador dá shuffle á ordem dos conjuntos antes de enviá-los para o server, para
        // este não descobrir qual contém o índice da mensagem desejada
        std::shuffle(aleatory_Partition.begin(), aleatory_Partition.end(), shuffle_random);
        
        print_Sets(aleatory_Partition,g);
    }


    // Função que recupera a mensagem desejada pelo user do conjunto de mensagens da base de dados
    void recoverMsgAnsServer(std::vector<Set> aleatory_Partition,std::vector<int> answers, std::vector<int> S, std::vector <int> database_messages ,int g, int W)
    {
        int idxP,val;

        // Descobrir qual o conjunto onde está o índice da mensagem desejada
        for(int i = 0; i < g; i++)
        {
            if(aleatory_Partition[i].hasindexW)
            {
                idxP = i;
            }
        }

        val = answers[idxP];

        // Nesse conjunto comparar os seus índices com os índices da side information e se for igual remover o seu respetivo valor da correspondente resposta do server
        for(size_t i = 0; i < aleatory_Partition[idxP].values.size(); i++)
        {
            for(size_t j = 0; j < S.size(); j++)
            {
                if(aleatory_Partition[idxP].values[i] == S[j])
                {
                    val -= database_messages[S[j]];
                }
            }
        }

        std::cout << "\nValue from database with index " << W << " is: " << val << '\n';
    }



/*
    K -> Número de mensagens
    M -> Tamanho do conjunto de mensagens de side information
    W -> Índice da mensagem a recuperar
    S -> Conjunto de índices das mensagens de side information (não contém W)
    g -> Número de conjuntos 
*/

    int main ()
    {
        int K = 12,M = 2, W = 2, g;
        std::vector<int> S = {4,6};
        bool isSpecialCase = false;
        ZZ maxVal = conv<ZZ>(10000);

        // Calcular número de conjuntos
        g = (int) round((double) K / (M + 1));
        
        // Averiguar se estamos no caso especial ou base
        // K == M+1
        if(K % (M+1) == 0)
            isSpecialCase = true;

        // Preencher mensagens com valores aleatórios
        std::vector<int> database_messages (K);
        for(int i = 0; i < K; i++)
        {
            ZZ random = RandomBnd(maxVal);
            database_messages[i] = conv<int>(random);
        }

        std::cout << "\nPrint database messages: " << '\n';
        print_Messages(database_messages);

        std::random_device rd;  
        std::mt19937 shuffle_random(rd());
        std::vector<Set> aleatory_Partition(g);
        std::vector<int> answers (g);

        std::cout << '\n';
        // Função utilizada para construir os g conjuntos
        build_Partition(K,M,W,g,isSpecialCase,S,aleatory_Partition, shuffle_random,maxVal);
        
        // Função que calcula as respostas dos servers
        serverAnswers(aleatory_Partition,g,answers,database_messages);
        std::cout << "\nServer Responses: " << '\n';
        for(int i = 0; i < g; i++)
        {
            std::cout << "  " << answers[i] << '\n';
        }

        recoverMsgAnsServer(aleatory_Partition,answers,S,database_messages,g,W);

        return 0;
    }