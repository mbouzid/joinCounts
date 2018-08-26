##

CXX	    = g++
CFLAGS	    = -std=c++11 -Wall -Wno-unused-function -fopenmp
OFLAGS	    = 
INC	    = 
LIB	    = -pthread -lpthread -lz
#PROF	    = -pg

##

OBJ	    = obj/main.o obj/tool.o
TGT	    = bin/joinCounts

##

all:  $(TGT)
	@echo "building done"

##

$(TGT): $(OBJ)
	$(CXX) $(PROF) $(CFLAGS) $(OFLAGS) -o $@ $^ $(LIB)

##

obj/%.o:  src/%.cpp
	$(CXX) $(PROF) $(INC) -c $(CFLAGS) $(OFLAGS) -o $@ $< $(LIB)



##

clean:	$(TGT) $(OBJ)
	rm $^
