BitRank.o: BitRank.cpp BitRank.h BlockArray.h Tools.h
EditDistance.o: EditDistance.cpp EditDistance.h Tools.h
TCImplementation.o: TCImplementation.cpp TCImplementation.h BitRank.h \
 BlockArray.h Tools.h TextCollection.h ArrayDoc.h TextStorage.h
TextCollection.o: TextCollection.cpp TextCollection.h Tools.h \
 TCImplementation.h BitRank.h BlockArray.h ArrayDoc.h TextStorage.h
TextCollectionBuilder.o: TextCollectionBuilder.cpp incbwt/rlcsa_builder.h \
 incbwt/rlcsa.h incbwt/bits/vectors.h incbwt/bits/deltavector.h \
 incbwt/bits/bitvector.h incbwt/bits/../misc/definitions.h \
 incbwt/bits/bitbuffer.h incbwt/bits/rlevector.h incbwt/sasamples.h \
 incbwt/misc/definitions.h incbwt/bits/bitbuffer.h \
 incbwt/bits/deltavector.h incbwt/misc/parameters.h \
 incbwt/misc/definitions.h TextCollectionBuilder.h TextCollection.h \
 Tools.h TextStorage.h TCImplementation.h BitRank.h BlockArray.h \
 ArrayDoc.h
TextStorage.o: TextStorage.cpp TextStorage.h TextCollection.h Tools.h
Tools.o: Tools.cpp Tools.h
builder.o: builder.cpp TextCollectionBuilder.h TextCollection.h Tools.h \
 TextStorage.h
maxoverlaps.o: maxoverlaps.cpp
sfoverlap.o: sfoverlap.cpp TextCollectionBuilder.h TextCollection.h \
 Tools.h TextStorage.h EditDistanceOverlap.h EditDistance.h \
 OverlapQuery.h TCImplementation.h BitRank.h BlockArray.h ArrayDoc.h \
 MismatchOverlap.h SuffixFilterOverlap.h SuffixFilterOverlapEd.h
