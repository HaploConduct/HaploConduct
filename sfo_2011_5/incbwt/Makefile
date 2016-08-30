CC = g++
CPPFLAGS = -Wall -O3
OBJS = rlcsa.o rlcsa_builder.o sasamples.o bits/bitvector.o bits/deltavector.o bits/rlevector.o bits/vectors.o misc/parameters.o misc/utils.o qsufsort/qsufsort.o


all: library

library: $(OBJS)
	ar rc rlcsa.a $(OBJS)

depend:
	g++ -MM *.cpp bits/*.cpp misc/*.cpp qsufsort/*.c > dependencies.mk

update:
	cp ../rlcsa.* ../rlcsa_builder.* ../sasamples.* .
	cp ../bits/* bits/
	cp ../misc/* misc/
	cp ../qsufsort/* qsufsort/

clean:
	rm -f rlcsa.a
	rm -f $(OBJS)


include dependencies.mk
