CC=g++
LDFLAGS=-lpthread
CFLAGS=-O3 -D NDEBUG -Wall $(CONFIGO)

# Configure COIN-OR libs
CONFIGO=$(shell PKG_CONFIG_PATH=$$CBC_PATH/lib/pkgconfig:$$CBC_PATH/share/pkgconfig: /usr/bin/pkg-config --cflags cbc)

SRC= $(wildcard ../src/*.cpp)
OBJT= $(SRC:.cpp=.o)
OBJ= $(subst ../src/,,$(OBJT))

all: $(OBJ)

#debug: CFLAGS=-g -D DEBUG -lpthread
debug: CFLAGS=-g $(CONFIGO) -D DEBUG -D DEBUG_FUNC
debug: $(OBJ)

prof: CFLAGS=-g -pg $(CONFIGO)
prof: $(OBJ)

%.o: ../src/%.cpp
	$(CC) -c $< $(CFLAGS) -o $@

clean:
	rm -f $(OBJ) *.o
