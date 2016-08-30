#include <cstdlib>

#include "deltavector.h"


namespace CSA
{


DeltaVector::DeltaVector(std::ifstream& file) :
  BitVector(file)
{
}
DeltaVector::DeltaVector(std::FILE * file) :
  BitVector(file)
{
}


DeltaVector::DeltaVector(DeltaEncoder& encoder, usint universe_size) :
  BitVector(encoder, universe_size)
{
}

DeltaVector::~DeltaVector()
{
}

//--------------------------------------------------------------------------

usint
DeltaVector::rank(usint value, bool at_least)
{
  if(value >= this->size) { return this->items; }
  this->getSample(this->sampleForValue(value));

  while(this->cur < this->block_items && this->val < value)
  {
    this->val += this->buffer->readDeltaCode();
    this->cur++;
  }

  usint idx = this->sample.first + this->cur + 1;
  if(!at_least && this->val > value) { idx--; }
  if(at_least && this->val < value)  { this->getSample(this->block + 1); }
  return idx;
}

usint
DeltaVector::select(usint index)
{
  if(index >= this->items) { return this->size; }
  this->getSample(this->sampleForIndex(index));

  usint lim = index - this->sample.first;
  for(; this->cur < lim; this->cur++)
  {
    this->val += this->buffer->readDeltaCode();
  }

  return this->val;
}

usint
DeltaVector::selectNext()
{
  if(this->cur >= this->block_items)
  {
    this->getSample(this->block + 1);
    return this->val;
  }

  this->cur++;
  this->val += this->buffer->readDeltaCode();
  return this->val;
}

pair_type
DeltaVector::valueAfter(usint value)
{
  if(value >= this->size) { return pair_type(this->size, this->items); }
  this->getSample(this->sampleForValue(value));

  while(this->cur < this->block_items && this->val < value)
  {
    this->val += this->buffer->readDeltaCode();
    this->cur++;
  }
  if(this->val < value)
  {
    this->getSample(this->block + 1);
  }

  return pair_type(this->val, this->sample.first + this->cur);
}

pair_type
DeltaVector::nextValue()
{
  if(this->cur >= this->block_items)
  {
    this->getSample(this->block + 1);
    return pair_type(this->val, this->sample.first);
  }

  this->cur++;
  this->val += this->buffer->readDeltaCode();
  return pair_type(this->val, this->sample.first + this->cur);
}

pair_type
DeltaVector::selectRun(usint index, usint max_length)
{
  return pair_type(this->select(index), 0);
}

pair_type
DeltaVector::selectNextRun(usint max_length)
{
  return pair_type(this->selectNext(), 0);
}

bool
DeltaVector::isSet(usint value)
{
  if(value >= this->size) { return false; }
  this->getSample(this->sampleForValue(value));

  while(this->cur < this->block_items && this->val < value)
  {
    this->val += this->buffer->readDeltaCode();
    this->cur++;
  }

  return (this->val == value);
}

//--------------------------------------------------------------------------

usint
DeltaVector::reportSize()
{
  usint bytes = sizeof(*this);
  bytes += BitVector::reportSize();
  return bytes;
}

//--------------------------------------------------------------------------

DeltaEncoder::DeltaEncoder(usint block_bytes, usint superblock_size) :
  VectorEncoder(block_bytes, superblock_size)
{
}

DeltaEncoder::~DeltaEncoder()
{
}

void
DeltaEncoder::setBit(usint value)
{
  if(this->items == 0)
  {
    this->setFirstBit(value);
    return;
  }
  if(value < this->size) { return; }

  usint diff = value + 1 - this->size;
  this->size = value + 1;
  this->items++;
  if(this->buffer->writeDeltaCode(diff)) { return; }

  // Didn't fit into the block. A new sample & block required.
  this->addNewBlock();
}


} // namespace CSA
