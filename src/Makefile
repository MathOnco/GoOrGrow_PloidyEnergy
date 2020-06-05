# As usual, you can build this software by typing:
#
#     make


# ----- Make Macros -----

CXXFLAGS =	-O2
CXX	=	g++

TARGETS =	go_or_grow
OBJECTS =	go_or_grow_V2.cpp
REMOVE =	go_or_grow.o
LINKS =	-std=c++11 -larmadillo

# ----- Make Rules -----

all:	$(TARGETS)

go_or_grow:	$(OBJECTS)
	$(CXX) $(CXXFLAGS) -o $(TARGETS) $(OBJECTS) $(LINKS)

clean:
	rm -f $(TARGETS)
