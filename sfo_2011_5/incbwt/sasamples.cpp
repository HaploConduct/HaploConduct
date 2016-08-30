#include <algorithm>
#include <fstream>
#include <iostream>

#include "sasamples.h"
#include "misc/utils.h"


namespace CSA
{


SASamples::SASamples(std::ifstream& sample_file, usint sample_rate) :
  rate(sample_rate)
{
  this->indexes = new DeltaVector(sample_file);

  this->size = indexes->getSize();
  this->items = indexes->getNumberOfItems();
  this->integer_bits = length(this->items - 1);

  this->samples = new FastBitBuffer(sample_file, this->items, this->integer_bits);
  this->buildInverseSamples();
}

SASamples::SASamples(usint* array, usint data_size, usint sample_rate) :
  rate(sample_rate),
  size(data_size)
{
  this->items = (this->size - 1) / this->rate + 1;
  this->integer_bits = length(this->items - 1);

  this->samples = new FastBitBuffer(this->items, this->integer_bits);
  DeltaEncoder encoder(SASamples::BLOCK_SIZE);
  for(usint i = 0; i < this->size; i++)
  {
    if(array[i] % rate == 0)
    {
      encoder.setBit(i);
      this->samples->writeItem(array[i] / sample_rate);
    }
  }
  this->indexes = new DeltaVector(encoder, this->size);

  this->buildInverseSamples();
}

SASamples::SASamples(SASamples& index, SASamples& increment, usint* positions, usint number_of_positions) :
  rate(index.rate),
  size(index.size + increment.size),
  items(index.items + increment.items)
{
  this->mergeSamples(index, increment, positions, number_of_positions);

  index.indexes = 0;
  index.samples = 0;
  index.inverse_samples = 0;
  increment.indexes = 0;
  increment.samples = 0;
  increment.inverse_samples = 0;

  this->buildInverseSamples();
}

SASamples::~SASamples()
{
  delete this->indexes;
  delete this->samples;
  delete this->inverse_samples;
}

void
SASamples::writeTo(std::ofstream& sample_file)
{
  this->indexes->writeTo(sample_file);
  this->samples->writeBuffer(sample_file, false);
}

//--------------------------------------------------------------------------

usint
SASamples::inverseSA(usint value)
{
  if(value >= this->size) { return this->size; }
  usint i = this->inverse_samples->readItem(value / this->rate);
  return this->indexes->select(i);
}

pair_type
SASamples::getFirstSampleAfter(usint index)
{
  if(index >= this->size) { return pair_type(this->size, this->size); }

  return this->indexes->valueAfter(index);
}

//--------------------------------------------------------------------------

usint
SASamples::reportSize()
{
  usint bytes = sizeof(*this) - sizeof(this->indexes);
  bytes += this->indexes->reportSize();
  bytes += this->samples->reportSize();
  bytes += this->inverse_samples->reportSize();
  return bytes;
}

//--------------------------------------------------------------------------

void
SASamples::buildInverseSamples()
{
  this->inverse_samples = new FastBitBuffer(this->items, this->integer_bits);
  this->samples->goToItem(0);
  for(usint i = 0; i < this->items; i++)
  {
    this->inverse_samples->goToItem(this->samples->readItem());
    this->inverse_samples->writeItem(i);
  }
}


//--------------------------------------------------------------------------


void
SASamples::mergeSamples(SASamples& index, SASamples& increment, usint* positions, usint n)
{
  DeltaVector* first = index.indexes;
  DeltaVector* second = increment.indexes;
  FastBitBuffer* first_samples = index.samples;
  FastBitBuffer* second_samples = increment.samples;

  usint first_bit = first->select(0);
  bool first_finished = false;
  usint second_bit = second->select(0);
  usint sum = index.items;
  first_samples->goToItem(0);
  second_samples->goToItem(0);

  DeltaEncoder encoder(SASamples::BLOCK_SIZE);
  this->integer_bits = length(this->items - 1);
  this->samples = new FastBitBuffer(this->items, this->integer_bits);
  for(usint i = 0; i < n; i++, first_bit++)
  {
    while(!first_finished && first_bit < positions[i])
    {
      encoder.setBit(first_bit);
      this->samples->writeItem(first_samples->readItem());
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
      this->samples->writeItem(second_samples->readItem() + sum);
      second_bit = second->selectNext();
    }
  }

  while(!first_finished)
  {
    encoder.setBit(first_bit);
    this->samples->writeItem(first_samples->readItem());
    if(!first->hasNext()) { break; }
    first_bit = first->selectNext() + n;
  }

  delete index.indexes;
  delete index.samples;
  delete index.inverse_samples;
  delete increment.indexes;
  delete increment.samples;
  delete increment.inverse_samples;

  this->indexes = new DeltaVector(encoder, size);
}


} // namespace CSA
