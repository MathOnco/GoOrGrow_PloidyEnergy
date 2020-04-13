# As usual, you can build this software by typing:
#
#     make


# ----- Make Macros -----

CXXFLAGS =	-O2
CXX	=	g++

TARGETS =	go_or_grow
OBJECTS =	go_or_grow.cpp
REMOVE =	go_or_grow.o
LINKS =	-std=c++11 -larmadillo
#INCLUDE = -I/usr/include -I/usr/local/opt/lapack/include

# ----- Make Rules -----

all:	$(TARGETS)

go_or_grow:	$(OBJECTS)
	$(CXX) $(CXXFLAGS) $(INCLUDE) -o $(TARGETS) $(OBJECTS) $(LINKS)

clean:
	rm -f $(TARGETS)
