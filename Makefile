MAINS1 = SPINSHALF_SSE
CXX = g++ 
CXX_FLAGS = -std=c++17 -O3 -Wall -Wextra -g 
CPP_FILES = $(wildcard src/*.cpp)
OBJ_FILES = $(addprefix obj/, $(notdir $(CPP_FILES:.cpp=.o)))
INCLUDE = -I./inc/

$(MAINS1): $(OBJ_FILES)
	$(CXX) $(CXX_FLAGS) $(INCLUDE) -o $@ $^

obj/%.o: src/%.cpp 
	$(CXX) $(CXX_FLAGS) $(INCLUDE) -c -o $@ $<

clean:	
	$(RM) $(OBJ_FILES) $(MAIN)


