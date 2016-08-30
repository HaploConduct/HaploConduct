#ifndef DEFINITIONS_H
#define DEFINITIONS_H

#include <algorithm>
#include <climits>


namespace CSA
{


#ifdef MASSIVE_DATA_RLCSA
typedef unsigned long usint;
typedef signed long   sint;
#else
typedef unsigned int  usint;
typedef signed int    sint;
#endif


#ifndef uchar
typedef unsigned char           uchar;
#endif
typedef std::pair<usint, usint> pair_type;


inline usint length(usint n)
{
  usint b = 0;
  while(n > 0) { b++; n >>= 1; }
  return b;
}

inline bool isEmpty(const pair_type& data)
{
  return (data.first > data.second);
}

inline usint length(const pair_type& data)
{
  return data.second + 1 - data.first;
}


const usint CHARS = ((usint)1 << CHAR_BIT);
const usint MEGABYTE = 1048576;
const usint MILLION  = 1000000;
const usint WORD_BITS = CHAR_BIT * sizeof(usint);
const usint WORD_MAX = ~((usint)0);

const pair_type EMPTY_PAIR = pair_type(1, 0);


// Previous GET was broken when BITS == WORD_BITS
// Current version works for usints and less
//#define GET(FIELD, BITS) ((FIELD) & ((1 << (BITS)) - 1))
#define GET(FIELD, BITS) ((FIELD) & (WORD_MAX >> (WORD_BITS - (BITS))))
#define LOWER(FIELD, N)  ((FIELD) >> (N))
#define HIGHER(FIELD, N) ((FIELD) << (N))

#define BITS_TO_BYTES(BITS) (((BITS) + CHAR_BIT - 1) / CHAR_BIT)
#define BYTES_TO_WORDS(BYTES) (((BYTES) + sizeof(usint) - 1) / sizeof(usint))
#define BITS_TO_WORDS(BITS) (((BITS) + WORD_BITS - 1) / WORD_BITS)


} // namespace CSA


#endif
