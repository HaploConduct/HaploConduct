#ifndef DELTAVECTOR_H
#define DELTAVECTOR_H

#include <fstream>

#include "bitvector.h"


namespace CSA
{


/*
  This class is used to construct a DeltaVector.
*/

class DeltaEncoder : public VectorEncoder
{
  public:
    DeltaEncoder(usint block_bytes, usint superblock_size = VectorEncoder::SUPERBLOCK_SIZE);
    ~DeltaEncoder();

    void setBit(usint value);
};


/*
  This is a gap-encoded bit vector using delta coding.
*/

class DeltaVector : public BitVector
{
  public:
    DeltaVector(std::ifstream& file);
    DeltaVector(std::FILE * file);
    DeltaVector(DeltaEncoder& encoder, usint universe_size);
    ~DeltaVector();

//--------------------------------------------------------------------------

    usint rank(usint value, bool at_least = false);

    usint select(usint index);
    usint selectNext();

    pair_type valueAfter(usint value);
    pair_type nextValue();

    pair_type selectRun(usint index, usint max_length);
    pair_type selectNextRun(usint max_length);

    bool isSet(usint value);

//--------------------------------------------------------------------------

    usint reportSize();
};


} // namespace CSA


#endif // DELTAVECTOR_H
