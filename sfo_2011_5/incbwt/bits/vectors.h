#ifndef VECTORS_H
#define VECTORS_H

#include "deltavector.h"
#include "rlevector.h"


namespace CSA
{


/*
  These functions merge two vectors using marked positions.
  The original vectors are deleted.
*/

RLEVector* mergeVectors(RLEVector* first, RLEVector* second, usint* positions, usint n, usint size, usint block_size);

DeltaVector* mergeVectors(DeltaVector* first, DeltaVector* second, usint* positions, usint n, usint size, usint block_size);


} // namespace CSA


#endif // VECTORS_H
