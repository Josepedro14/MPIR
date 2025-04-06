#include <iostream>
#include <vector>
#include <math.h>
#include <random>
#include <algorithm>

int idx, remainingSetSizes;
int indexC;

// Formato dos conjuntos
struct Set 
{
    std::vector<int> values;
    double begin_Prob;
    double end_Prob;
};


// Função para retornar o índice de um conjunto que vai ser escolhido aleatoriamente 
int chooseRandomNumPlusIndex (std::vector<Set> &aleatory_Partition, int g)
{
    double random_num = (double) rand() / RAND_MAX;

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


void print_Partitions(std::vector<Set> &aleatory_Partition, int g) {
    
    for (int i = 0; i < g; i++) {
        
        std::cout << "Partition " << i + 1 << ":" << '\n';
        std::cout << " Values: ";
        
        for (size_t j = 0; j < aleatory_Partition[i].values.size(); j++) {
                
            std::cout << aleatory_Partition[i].values[j] << " ";
        }

        std::cout << '\n';
    }

    std::cout << '\n';
}


void print_Messages (std::vector<int> &message)
{
    for(size_t i = 0; i < message.size(); i++)
    {
        std::cout << '\n' << message[i] << '\n';
    }
}


void shuffleMessages(std::vector<int> &messages, std::mt19937 shuffle_random)
{

    std::shuffle(messages.begin(),messages.end(), shuffle_random);
}


void serverAnswers(std::vector<Set> &aleatory_Partition, int g, std::vector<int> &answers)
{

    for(int i = 0; i < g; i++)
    {
        Set set = aleatory_Partition[i];

        for(size_t j = 0; j < set.values.size(); j++)
        {
            answers[i] += set.values[j] * 1; 
        }
    }
}



void build_Partition (int K, int M, int W, int g, bool is1stCase, std::vector<int> &S, std::vector<Set> &aleatory_Partition, std::mt19937 shuffle_random)
{
    std::vector<int> messages;

    for(int i = 0; i < K; i++)
    {
        messages.push_back(i + 1);
    }

    srand(time(NULL));


    if(is1stCase)
    {
        // Caso especial então P1 = W U S, e o resto dos g-1 conjuntos preenchidos de forma aleatória

        aleatory_Partition[0].values.push_back(W);
        messages.erase(std::remove(messages.begin(),messages.end(),W), messages.end());

        for(size_t i = 0; i < S.size(); i++)
        {
            aleatory_Partition[0].values.push_back(S[i]);
            messages.erase(std::remove(messages.begin(),messages.end(),S[i]), messages.end());
        }

        shuffleMessages(messages,shuffle_random);

        remainingSetSizes = K - (1 + S.size());
        idx = 0;

        for(int i = 1; i <= g-1; i++)
        {
            int size_Sets = std::min(remainingSetSizes,(M+1));
            remainingSetSizes -= size_Sets;

            for(int j = 0; j < size_Sets; j++)
            {
                aleatory_Partition[i].values.push_back(messages[idx]);
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

        // O utilizador escolhe aleatoriamente um desses conjuntos com base nas probabilidades de um número aleatório gerado entre 0 e 1 (através da função chooseRandomlyPlusIndex)
        int chosen_index = chooseRandomNumPlusIndex(aleatory_Partition, g);
        indexC = chosen_index + 1;

        // Se o conjunto escolhido for Pg então o utilizador coloca W dentro desse conjunto e preenche os restantes K - (g-1)(M+1) - 1 espaços com o conjunto de informação lateral
        // O resto dos conjuntos tem elementos escolhidos de forma aleatória dentro dos restantes elementos

        if(chosen_index == g-1)
        {

            //std::cout << "\nÍndice == g-1 " << "\n\n";

            aleatory_Partition[g-1].values.push_back(W);
            messages.erase(std::remove(messages.begin(),messages.end(),W), messages.end());

            while(aleatory_Partition[g-1].values.size() != (size_t) (K - (g-1) * (M+1)))
            {
                int random_indexS = rand() % S.size();
                aleatory_Partition[g-1].values.push_back(S[random_indexS]);
                messages.erase(std::remove(messages.begin(), messages.end(), S[random_indexS]), messages.end());
                S.erase(std::remove(S.begin(), S.end(), S[random_indexS]), S.end());
            }

            shuffleMessages(messages,shuffle_random);
            remainingSetSizes = K - aleatory_Partition[g-1].values.size();
            idx = 0;

            for(int i = 0; i <= g-2; i++)
            {
                int size_Sets = std::min(remainingSetSizes, (M+1));
                remainingSetSizes -= size_Sets;

                for(int j = 0; j < size_Sets; j++)
                {
                    aleatory_Partition[i].values.push_back(messages[idx]);
                    idx++;
                }
            }
            
        }

        else 
        {
            // Se o conjunto escolhido for um P ∈ P1,...,Pg-1, então o utilizador coloca a mensagem desejada W e o conjunto de informação lateral nesse conjunto;
            // Depois preenche os outros conjuntos aleatoriamente até que todas as mensagens sejam usadas

            //std::cout << "\nÍndice /= g-1 " <<"\n";

            aleatory_Partition[chosen_index].values.push_back(W);
            messages.erase(std::remove(messages.begin(),messages.end(),W), messages.end());

            for(size_t i = 0; i < S.size(); i++)
            {
                aleatory_Partition[chosen_index].values.push_back(S[i]);
                messages.erase(std::remove(messages.begin(),messages.end(),S[i]), messages.end());
            }

            shuffleMessages(messages,shuffle_random);

            remainingSetSizes = (int) messages.size();
            idx = 0;

            for(int i = 0; i <= g-1; i++)
            {
                if(i != chosen_index)
                {
                    int size_Sets = std::min(remainingSetSizes, (M+1));
                    remainingSetSizes -= size_Sets;

                    for(int j = 0; j < size_Sets; j++)
                    {
                        aleatory_Partition[i].values.push_back(messages[idx]);
                        idx++;
                    }
                }
            }

        }
    }

    // O utilizador dá shuffle á ordem dos conjuntos antes de enviá-los para o server
    std::shuffle(aleatory_Partition.begin(), aleatory_Partition.end(), shuffle_random);

    print_Partitions(aleatory_Partition,g);
    
}


int main ()
{
    int K = 8,M = 2, W = 2, g;
    std::vector<int> S = {4,6};
    bool is1stCase = false;

    g = (int) round((double) K / (M + 1));
    if(K == M+1)
        is1stCase = true;

    std::random_device rd;  
    std::mt19937 shuffle_random(rd());
    std::vector<Set> aleatory_Partition(g);
    std::vector<int> answers (g);

    build_Partition(K,M,W,g,is1stCase,S,aleatory_Partition, shuffle_random);
    
    serverAnswers(aleatory_Partition,g,answers);

    for(int i = 0;i < g; i++)
    {
        std::cout << "Server Response " << i+1 << " : " << answers[i] << '\n';
    }

    return 0;
}