CC = g++
CFLAGS = `root-config --cflags` -std=c++11
LDFLAGS = `root-config --libs`
objects := $(patsubst %.cpp,%.o,$(patsubst ./src/%,./bin/%,$(wildcard ./src/*.cpp)))

all: CFLAGS += -O3
all: rootAnalyzer

rootAnalyzer: rootAnalyzer.cpp $(objects)
	$(CC) $(CFLAGS) -o rootAnalyzer rootAnalyzer.cpp $(objects) $(LDFLAGS)

bin/%.o: src/%.cpp
	$(CC) $(CFLAGS) -c -o $(patsubst ./src/%,./bin/%,$@) $<

clean:
	rm rootAnalyzer bin/*.o

