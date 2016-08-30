#include <algorithm>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <iostream>

#include "vectors.h"
#include "../misc/utils.h"


namespace CSA
{


inline void
handleOne(RLEEncoder& encoder, pair_type& run, usint position)
{
  if(run.second == 0)
  {
    run.first = position;
    run.second = 1;
    return;
  }
  if(position == run.first + run.second)
  {
    run.second++;
    return;
  }
  encoder.setRun(run.first, run.second);
  run.first = position;
  run.second = 1;
}

inline void
handleRun(RLEEncoder& encoder, pair_type& run, pair_type& next, usint limit)
{
  if(run.second == 0)
  {
    run.first = next.first;
    run.second = std::min(limit - run.first, next.second);
    next.first += run.second;
    next.second -= run.second;
    return;
  }

  if(next.first == run.first + run.second)
  {
    usint cont = std::min(limit - next.first, next.second);
    run.second += cont;
    next.first += cont;
    next.second -= cont;
    return;
  }

  encoder.setRun(run.first, run.second);
  run.first = next.first;
  run.second = std::min(limit - run.first, next.second);;
  next.first += run.second;
  next.second -= run.second;
}

RLEVector*
mergeVectors(RLEVector* first, RLEVector* second, usint* positions, usint n, usint size, usint block_size)
{
  if((first == 0 && second == 0) || positions == 0) { return 0; }

  pair_type first_run;
  bool first_finished;
  if(first == 0)
  {
    first_run = pair_type(size, 0);
    first_finished = true;
  }
  else
  {
    first_run = first->selectRun(0, size);
    first_run.second++;
    first_finished = false;
  }

  usint second_bit;
  if(second == 0)
  {
    second_bit = n;
  }
  else
  {
    second_bit = second->select(0);
  }

  RLEEncoder encoder(block_size);
  pair_type run = pair_type(size, 0);
  for(usint i = 0; i < n; i++, first_run.first++)
  {
    while(!first_finished && first_run.first < positions[i])
    {
      handleRun(encoder, run, first_run, positions[i]);
      if(first_run.second == 0)
      {
        if(first->hasNext())
        {
          first_run = first->selectNextRun(size);
          first_run.first += i;
          first_run.second++;
        }
        else
        {
          first_finished = true;
        }
      }
    }

    if(i == second_bit) // positions[i] is one
    {
      handleOne(encoder, run, positions[i]);
      second_bit = second->selectNext();
    }
    else  // positions[i] is zero
    {
      if(run.second != 0)
      {
        encoder.setRun(run.first, run.second);
        run.second = 0;
      }
    }
  }

  while(!first_finished)
  {
    handleRun(encoder, run, first_run, size);
    if(first->hasNext())
    {
      first_run = first->selectNextRun(size);
      first_run.first += n;
      first_run.second++;
    }
    else { break; }
  }

  if(run.second != 0)
  {
    encoder.setRun(run.first, run.second);
  }

  delete first; delete second;
  return new RLEVector(encoder, size);
}

//--------------------------------------------------------------------------

DeltaVector*
mergeVectors(DeltaVector* first, DeltaVector* second, usint* positions, usint n, usint size, usint block_size)
{
  if((first == 0 && second == 0) || positions == 0) { return 0; }

  usint first_bit;
  bool first_finished;
  if(first == 0)
  {
    first_bit = 0;
    first_finished = true;
  }
  else
  {
    first_bit = first->select(0);
    first_finished = false;
  }

  usint second_bit;
  if(second == 0)
  {
    second_bit = n;
  }
  else
  {
    second_bit = second->select(0);
  }

  DeltaEncoder encoder(block_size);
  for(usint i = 0; i < n; i++, first_bit++)
  {
    while(!first_finished && first_bit < positions[i])
    {
      encoder.setBit(first_bit);
      if(first->hasNext())
      {
        first_bit = first->selectNext() + i;
      }
      else
      {
        first_finished = true;
      }
    }

    if(i == second_bit) // positions[i] is one
    {
      encoder.setBit(positions[i]);
      second_bit = second->selectNext();
    }
  }

  while(!first_finished)
  {
    encoder.setBit(first_bit);
    if(!first->hasNext()) { break; }
    first_bit = first->selectNext() + n;
  }

  delete first; delete second;
  return new DeltaVector(encoder, size);
}


} // namespace CSA
