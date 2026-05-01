CXX = g++
CXXFLAGS = -O0 -pg  -mfpmath=sse -fstack-protector-all -g -W -Wcast-align -Werror -pedantic -pedantic-errors -Wfloat-equal -Wpointer-arith -Wformat-security -Wmissing-format-attribute -Wformat=1 -Wwrite-strings -Wcast-align -Wno-long-long -Woverloaded-virtual -Wnon-virtual-dtor -Wcast-qual -Wno-suggest-attribute=format
CXXFLAGS += -Wall -Wextra -Wunused-parameter -Wno-error=unused-parameter
LDFLAGS = -lm -pthread

TARGET = a.out
SOURCES =  main/main.cpp \
main/thread.cpp \
main/time.cpp \
main/reduce_sum.cpp \
fill_msr_matrix/fill_msr_matrix.cpp \
fill_msr_matrix/calculation_сourrant_basis.cpp \
fill_msr_matrix/fill_right_side.cpp \
residual/claculate_residual.cpp \
residual/functions.cpp \
residual/calculate_approximation_fun.cpp \
check_alg/check_symm_and_rov_sum.cpp \
solve_msr_system/minimal_errors_method.cpp \
solve_msr_system/operations_matrix_or_vec.cpp \
solve_msr_system/symm_seidel_preconditioner.cpp
		
OBJECTS = $(SOURCES:.cpp=.o)
HEADERS = header.h

all: $(TARGET)

$(TARGET): $(OBJECTS)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(OBJECTS) $(LDFLAGS)

%.o: %.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(OBJECTS) $(TARGET)

rebuild: clean all

debug: CXXFLAGS += -DDEBUG
debug: rebuild


.PHONY: all clean rebuild run debug
