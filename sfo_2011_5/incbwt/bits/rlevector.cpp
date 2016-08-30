#include <cstdlib>

#include "rlevector.h"
#include "../misc/utils.h"


namespace CSA
{


RLEVector::RLEVector(std::ifstream& file) :
  BitVector(file)
{
}

RLEVector::RLEVector(RLEEncoder& encoder, usint universe_size) :
  BitVector(encoder, universe_size)
{
}

RLEVector::~RLEVector()
{
}

//--------------------------------------------------------------------------

usint
RLEVector::rank(usint value, bool at_least)
{
  if(value >= this->size) { return this->items; }

  this->valueLoop(value);

  usint idx = this->sample.first + this->cur + 1;
  if(!at_least && this->val > value)
  {
    idx--;
  }
  if(at_least && this->val < value)
  {
    this->getSample(this->block + 1);
    this->run = 0;
    idx = this->sample.first + this->cur + 1;
  }
  return idx;
}

usint
RLEVector::select(usint index)
{
  if(index >= this->items) { return this->size; }
  this->getSample(this->sampleForIndex(index));
  this->run = 0;

  usint lim = index - this->sample.first;
  while(this->cur < lim)
  {
    this->val += this->buffer->readDeltaCode();
    usint temp = this->buffer->readDeltaCode();
    this->val += temp - 1;
    this->cur += temp;
  }
  if(this->cur > lim)
  {
    this->run = this->cur - lim;
    this->cur -= this->run;
    this->val -= this->run;
  }

  return this->val;
}

usint
RLEVector::selectNext()
{
  if(this->cur >= this->block_items)
  {
    this->getSample(this->block + 1);
    this->run = 0;
    return this->val;
  }

  this->cur++;
  if(this->run > 0)
  {
    this->val++;
    this->run--;
  }
  else
  {
    this->val += this->buffer->readDeltaCode();
    this->run = this->buffer->readDeltaCode() - 1;
  }

  return this->val;
}

pair_type
RLEVector::valueAfter(usint value)
{
  if(value >= this->size) { return pair_type(this->size, this->items); }

  this->valueLoop(value);

  if(this->val < value)
  {
    this->getSample(this->block + 1);
    this->run = 0;
  }

  return pair_type(this->val, this->sample.first + this->cur);
}

pair_type
RLEVector::nextValue()
{
  if(this->cur >= this->block_items)
  {
    this->getSample(this->block + 1);
    this->run = 0;
    return pair_type(this->val, this->sample.first);
  }

  this->cur++;
  if(this->run > 0)
  {
    this->val++;
    this->run--;
  }
  else
  {
    this->val += this->buffer->readDeltaCode();
    this->run = this->buffer->readDeltaCode() - 1;
  }

  return pair_type(this->val, this->sample.first + this->cur);
}

pair_type
RLEVector::selectRun(usint index, usint max_length)
{
  usint value = this->select(index);

  usint len = std::min(max_length, this->run);
  this->run -= len; this->cur += len; this->val += len;

  return pair_type(value, len);
}

pair_type
RLEVector::selectNextRun(usint max_length)
{
  usint value = this->selectNext();

  usint len = std::min(max_length, this->run);
  this->run -= len; this->cur += len; this->val += len;

  return pair_type(value, len);
}

bool
RLEVector::isSet(usint value)
{
  if(value >= this->size) { return false; }

  this->valueLoop(value);

  return (this->val == value);
}

//--------------------------------------------------------------------------

usint
RLEVector::reportSize()
{
  usint bytes = sizeof(*this);
  bytes += BitVector::reportSize();
  return bytes;
}

//--------------------------------------------------------------------------

RLEEncoder::RLEEncoder(usint block_bytes, usint superblock_size) :
  VectorEncoder(block_bytes, superblock_size)
{
}

RLEEncoder::~RLEEncoder()
{
}

void
RLEEncoder::setRun(usint start, usint len)
{
  if(this->items == 0)
  {
    this->setFirstBit(start);
    if(len > 1)
    {
      this->RLEncode(1, len - 1);
    }
    return;
  }
  if(start < this->size || len == 0) { return; }

  // Write as much into the buffer as possible.
  usint diff = start + 1 - this->size;
  usint free_bits = this->buffer->bitsLeft();
  usint code_bits = this->buffer->deltaCodeLength(diff);
  if(free_bits > code_bits) // At least a part of the run fits into the block.
  {
    free_bits -= code_bits;
    usint run_bits = this->buffer->deltaCodeLength(len);
    if(run_bits <= free_bits)
    {
      this->RLEncode(diff, len);
      return;
    }

    // Encode as much as possible and let the rest spill.
    usint llen = 1;
    while(llen + 2 * length(llen + 1) - 1 <= free_bits) { llen++; }
    llen = ((usint)1 << llen) - 1;

    this->RLEncode(diff, llen);
    len -= llen;

    // A new sample will be added.
    this->size++;
    this->items++;
  }
  else
  {
    this->size = start + 1;
    this->items++;
  }

  // Didn't fit into the block. A new sample & block required.
  this->addNewBlock();
  if(len > 1)
  {
    this->RLEncode(1, len - 1);
  }
}


} // namespace CSA
