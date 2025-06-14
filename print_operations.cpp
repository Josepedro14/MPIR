#include "print_operations.h"


void show_MessagesSubpackets (std::vector <Message> &messages, int K, int L, int symbols_subpacket)
{
    // Percorre as mensagens 
    for(int i = 0; i < K; i++)
    {
        std::cout << "Message " << i+1 << ": " << '\n';
        //  Itera sobre os subpacotes da mensagem
        for(int j = 0; j < L; j++)
        {
            std::cout << "Subpacket " << j+1 << ": " << '\n';
            // Printa os elementos dos subpacotes da mensagem atual
            for(int k = 0; k < symbols_subpacket; k++)
            {
                std::cout << "Value: " << messages[i].subpackets[j].numsR_finite_field[k] << '\n';
            }
        }

        std::cout << '\n';
    }
}



void show_Matrix (Eigen::MatrixXd &matrix, int rows, int cols)
{
    std::cout << "Matriz :" << '\n';

    for(int line = 0; line < rows; line++)
    {
        for(int col = 0; col < cols; col++)
        {
            std::cout << matrix(line,col) << " ";
        }

        std::cout << '\n';
    }
}



void show_Vector (Eigen::VectorXd &vec, int size)
{
    std::cout << "Vector: " << '\n';

    for(int i = 0; i < size; i++)
    {
        std::cout << vec(i) << " " << '\n';
    }
}



void show_Vectorxi (Eigen::VectorXi &vec, int size)
{
    std::cout << "Vector: " << '\n';

    for(int i = 0; i < size; i++)
    {
        std::cout << vec(i) << " ";
    }

    std::cout << '\n';
}



void show_VectorSTD (std::vector <int> &vec, int size)
{
    std::cout << "\nVector std: " << '\n';

    for(int i = 0; i < size; i++)
    {
        std::cout << vec[i] << " ";
    }

    std::cout << '\n';
}



void show_SubpacketsVectorXi (std::vector <Eigen::VectorXi> messages, size_t size, int symbols_subpacket)
{   
    // Iteramos sobre os subpacotes das mensagens de interesse ou interferência e damos print aos valores que os mesmos contêm
     for(size_t i = 0; i < messages.size(); i++)
        {
            std::cout << "\nPrint subpacotes das mensagens: " << '\n';
            for(int j = 0; j < symbols_subpacket; j++)
            {
                std::cout << messages[i](j) << " ";
            }
        }
}



void show_QuerysAndAnswersServer (std::vector <int> &server_indexs, std::vector <Eigen::VectorXi> &n_vectors, std::vector <Eigen::VectorXi> &Y_vectors, int N, int L, int K, int symbols_subpacket)
{
    // Print comunicações com os servers
    for(int i = 0; i < N; i++)
    {
        std::cout << "\n--------------------------------------------------------- SERVER " << i+1 << " ---------------------------------------------------------\n";
        // Aqui estamos a enviar para o server indicado o vetor vn como o Query (Qn^[W]) ao server (n -> server_indexs[i])
        std::cout << "\nQuerry: " << '\n';
        show_Vectorxi(n_vectors[server_indexs[i]], K*L);
        // Aqui o server (n -> server_indexs[i]) vai devolver/calcular uma Answer (An^[W]) que é a combinação linear Yn para o vetor vn recebido
        std::cout << "\nAnswer: " << '\n';
        show_Vectorxi(Y_vectors[server_indexs[i]], symbols_subpacket);
        std::cout << "\n----------------------------------------------------------------------------------------------------------------------------\n";
    }
}