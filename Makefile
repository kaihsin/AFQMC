CYTNX_PATH := /home/kaywu/Cytnx_test
CYTNX_LDFLAGS := -lmkl_intel_ilp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -ldl -lm 
CC := g++
CCFLAGS := -std=c++11 -I$(CYTNX_PATH)/include 
LDFLAGS := -L$(CYTNX_PATH)/lib $(CYTNX_LDFLAGS)
all: dqmc

dqmc: dqmc.o Parser.o
	$(CC) $^ $(CYTNX_PATH)/lib/libcytnx.a -o $@ $(LDFLAGS)


dqmc.o: dqmc.cpp dqmc.hpp
	$(CC) -c $(CCFLAGS) $< 

Parser.o: Parser.cpp Parser.hpp
	$(CC) -c $(CCFLAGS) $<

.phony: clean

clean:
	rm dqmc *.o

