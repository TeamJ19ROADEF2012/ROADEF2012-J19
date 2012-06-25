CC=g++
LDFLAGS=-lpthread $(CONFIGEX)
CFLAGS=-O3 -D NDEBUG -Wall $(CONFIGO)

# Configure COIN-OR libs
CONFIGO=$(shell PKG_CONFIG_PATH=$$CBC_PATH/lib/pkgconfig:$$CBC_PATH/share/pkgconfig: /usr/bin/pkg-config --cflags cbc)
CONFIGEX=$(shell PKG_CONFIG_PATH=$$CBC_PATH/lib/pkgconfig:$$CBC_PATH/share/pkgconfig: /usr/bin/pkg-config --libs cbc)

SRC= $(wildcard ../src/*.cpp)
OBJT= $(SRC:.cpp=.o)
OBJ= $(subst /src/,/obj/,$(OBJT))
OBJN= $(subst ../obj/machineReassignment.o,,$(OBJ))
OBJM= $(subst ../obj/test_neighborhoodReassignment.o,,$(OBJ))

all: $(OBJ) machineReassignment neighborhoodReassignment

debug: CFLAGS=-g -D DEBUG -D DEBUG_FUNC $(CONFIGO)
debug: $(OBJ) neighborhoodReassignment machineReassignment

prof: CFLAGS=-g -pg $(CONFIGO)
prof: $(OBJ) machineReassignment neighborhoodReassignment

%.o: ../src/%.cpp
	cd ../obj/; make $@

machineReassignment: $(OBJ)
	$(CC) $(CFLAGS) -o $@ $(OBJM) $(LDFLAGS)

neighborhoodReassignment: $(OBJ)
	$(CC) $(CFLAGS) -o $@ $(OBJN) $(LDFLAGS)

clean:
	rm -f neighborhoodReassignment machineReassignment
