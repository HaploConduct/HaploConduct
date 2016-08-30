#include <iostream>

#include "rlcsa_builder.h"


namespace CSA
{


RLCSABuilder::RLCSABuilder(usint _block_size, usint _sample_rate, usint _buffer_size) :
  block_size(_block_size), sample_rate(_sample_rate), buffer_size(_buffer_size),
  buffer(0)
{
  this->reset();
}

RLCSABuilder::~RLCSABuilder()
{
  delete this->index;
  delete[] this->buffer;
}

//--------------------------------------------------------------------------

void
RLCSABuilder::insertSequence(char* sequence, usint length, bool delete_sequence)
{
  if(sequence == 0 || length == 0 || !this->ok)
  {
    if(delete_sequence) { delete[] sequence; }
    return;
  }

  if(this->buffer == 0)
  {
    clock_t start = clock();
    RLCSA* temp = new RLCSA((uchar*)sequence, length, this->block_size, this->sample_rate, false, false);
    this->build_time += clock() - start;
    this->addRLCSA(temp, (uchar*)sequence, length + 1, delete_sequence);
    return;
  }

  if(this->buffer_size - this->chars > length)
  {
    memcpy(this->buffer + this->chars, sequence, length);
    if(delete_sequence) { delete[] sequence; }
    this->chars += length;
    this->buffer[this->chars] = 0;
    this->chars++;
  }
  else
  {
    this->flush();
    this->buffer = new uchar[this->buffer_size];
    if(length >= this->buffer_size - 1)
    {
      clock_t start = clock();
      RLCSA* temp = new RLCSA((uchar*)sequence, length, this->block_size, this->sample_rate, false, false);
      this->build_time += clock() - start;
      this->addRLCSA(temp, (uchar*)sequence, length + 1, delete_sequence);
    }
    else
    {
      memcpy(this->buffer + this->chars, sequence, length);
      if(delete_sequence) { delete[] sequence; }
      this->chars += length;
      this->buffer[this->chars] = 0;
      this->chars++;
    }
  }
}

RLCSA*
RLCSABuilder::getRLCSA()
{
  if(this->chars > 0) { this->flush(); }

  RLCSA* temp = this->index;
  this->reset();

  return temp;
}

char*
RLCSABuilder::getBWT(usint& length)
{
  if(this->chars > 0)
  {
    this->flush();
    if(this->buffer_size > 0) { this->buffer = new uchar[this->buffer_size]; }
  }

  if(this->index == 0 || !(this->ok))
  {
    length = 0;
    return 0;
  }

  length = this->index->getSize() + this->index->getNumberOfSequences();
  return (char*)(this->index->readBWT());
}

bool
RLCSABuilder::isOk()
{
  return this->ok;
}

double
RLCSABuilder::getBuildTime()
{
  return this->build_time / (double)CLOCKS_PER_SEC;
}

double
RLCSABuilder::getSearchTime()
{
  return this->search_time / (double)CLOCKS_PER_SEC;
}

double
RLCSABuilder::getMergeTime()
{
  return this->merge_time / (double)CLOCKS_PER_SEC;
}

//--------------------------------------------------------------------------

void
RLCSABuilder::flush()
{
  clock_t start = clock();
  RLCSA* temp = new RLCSA(this->buffer, this->chars, this->block_size, this->sample_rate, true, (this->index == 0));
  this->build_time += clock() - start;
  this->addRLCSA(temp, this->buffer, this->chars, true);
  this->buffer = 0; this->chars = 0;
}

void
RLCSABuilder::addRLCSA(RLCSA* increment, uchar* sequence, usint length, bool delete_sequence)
{
  if(this->index != 0)
  {
    clock_t start = clock();

    usint* positions = new usint[length];
    usint begin = 0;
    for(usint i = 0; i < length - 1; i++)
    {
      if(sequence[i] == 0)
      {
        this->index->reportPositions(&(sequence[begin]), i - begin, &(positions[begin]));
        begin = i + 1;
      }
    }
    this->index->reportPositions(&(sequence[begin]), length - 1 - begin, &(positions[begin]));

    std::sort(positions, positions + length);
    for(usint i = 0; i < length; i++)
    {
      positions[i] += i + 1;  // +1 because the insertion will be after positions[i]
    }
    if(delete_sequence) { delete[] sequence; }

    clock_t mark = clock();
    this->search_time += mark - start;

    RLCSA* merged = new RLCSA(*(this->index), *increment, positions, this->block_size);
    delete[] positions;
    delete this->index;
    delete increment;
    this->index = merged;

    this->merge_time += clock() - mark;  
  }
  else
  {
    this->index = increment;
  }

  this->ok &= this->index->isOk();
}

void
RLCSABuilder::reset()
{
  this->index = 0;

  if(this->buffer_size != 0)
  {
    this->buffer = new uchar[this->buffer_size];
  }
  this->chars = 0;

  this->ok = true;
}


} // namespace CSA
