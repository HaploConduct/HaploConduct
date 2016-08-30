rlcsa.o: rlcsa.cpp rlcsa.h bits/vectors.h bits/deltavector.h \
 bits/bitvector.h bits/../misc/definitions.h bits/bitbuffer.h \
 bits/rlevector.h sasamples.h misc/definitions.h bits/bitbuffer.h \
 bits/deltavector.h misc/parameters.h misc/definitions.h misc/utils.h \
 qsufsort/qsufsort.h qsufsort/../misc/definitions.h
rlcsa_builder.o: rlcsa_builder.cpp rlcsa_builder.h rlcsa.h bits/vectors.h \
 bits/deltavector.h bits/bitvector.h bits/../misc/definitions.h \
 bits/bitbuffer.h bits/rlevector.h sasamples.h misc/definitions.h \
 bits/bitbuffer.h bits/deltavector.h misc/parameters.h misc/definitions.h
sasamples.o: sasamples.cpp sasamples.h misc/definitions.h \
 bits/bitbuffer.h bits/../misc/definitions.h bits/deltavector.h \
 bits/bitvector.h bits/bitbuffer.h misc/utils.h misc/definitions.h
bitvector.o: bits/bitvector.cpp bits/bitvector.h \
 bits/../misc/definitions.h bits/bitbuffer.h
deltavector.o: bits/deltavector.cpp bits/deltavector.h bits/bitvector.h \
 bits/../misc/definitions.h bits/bitbuffer.h
rlevector.o: bits/rlevector.cpp bits/rlevector.h bits/bitvector.h \
 bits/../misc/definitions.h bits/bitbuffer.h bits/../misc/utils.h \
 bits/../misc/definitions.h
vectors.o: bits/vectors.cpp bits/vectors.h bits/deltavector.h \
 bits/bitvector.h bits/../misc/definitions.h bits/bitbuffer.h \
 bits/rlevector.h bits/../misc/utils.h bits/../misc/definitions.h
parameters.o: misc/parameters.cpp misc/parameters.h misc/definitions.h
utils.o: misc/utils.cpp misc/utils.h misc/definitions.h
qsufsort.o: qsufsort/qsufsort.c qsufsort/qsufsort.h \
 qsufsort/../misc/definitions.h
