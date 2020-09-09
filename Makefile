CYTNX_INC := $(shell python -c "exec(\"import cytnx\nprint(cytnx.__cpp_include__)\")")
CYTNX_LDFLAGS := $(shell python -c "exec(\"import cytnx\nprint(cytnx.__cpp_linkflags__)\")")
CYTNX_LIB := $(shell python -c "exec(\"import cytnx\nprint(cytnx.__cpp_lib__)\")")/libcytnx.a
CYTNX_CXXFLAGS := $(shell python -c "exec(\"import cytnx\nprint(cytnx.__cpp_flags__)\")")
CC := g++
CCFLAGS := -std=c++11 -I$(CYTNX_INC) ${CYTNX_CXXFLAGS}
all: dqmc

dqmc: dqmc.o Parser.o
	$(CC) $^ -o $@ $(CYTNX_LIB) $(CYTNX_LDFLAGS)


dqmc.o: dqmc.cpp dqmc.hpp
	$(CC) -c $(CCFLAGS) $< 

Parser.o: Parser.cpp Parser.hpp
	$(CC) -c $(CCFLAGS) $<

.phony: clean

clean:
	rm dqmc *.o

