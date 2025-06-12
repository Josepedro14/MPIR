# Compilador
CXX = g++

# Flags do Compilador
CXXFLAGS = -Wall -I /usr/include/eigen3

# Source Files
CPPSRCS = main.cpp math_operations.cpp finite_field_operations.cpp matrix_operations.cpp mpir_functions.cpp print_operations.cpp random_operations.cpp


# Flags de Linkagem
# NTL FLAGS
# LDFLAGS = -lntl -lgmp 

# Object Files
OBJS = $(CPPSRCS:.cpp=.o)
	
# Executável
TARGET = mpir

all: $(TARGET)
	./$(TARGET)


$(TARGET) : $(OBJS)
	$(CXX) $(CXXFLAGS) $(OBJS) -o $(TARGET) 

	
# Compilar source files em object files
%.o : %.cpp 
	$(CXX) $(CXXFLAGS) -c $< -o $@


# Eliminar Object files e executável
clean:
	rm -rf $(TARGET) $(OBJS)
