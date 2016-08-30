#include <cstdlib>

#include "bitvector.h"


namespace CSA
{


BitVector::BitVector(std::ifstream& file) :
  rank_index(0), select_index(0)
{
  file.read((char*)&(this->size), sizeof(this->size));
  file.read((char*)&(this->items), sizeof(this->items));
  file.read((char*)&(this->number_of_blocks), sizeof(this->number_of_blocks));
  file.read((char*)&(this->block_size), sizeof(this->block_size));

  this->array = new usint[this->block_size * this->number_of_blocks];
  file.read((char*)(this->array), this->block_size * this->number_of_blocks * sizeof(usint));
  this->buffer = new FastBitBuffer(this->array, this->block_size);

  this->integer_bits = length(this->size);
  this->samples = new FastBitBuffer(file, 2 * (this->number_of_blocks + 1), this->integer_bits);

  this->indexForRank();
  this->indexForSelect();
}

BitVector::BitVector(std::FILE * file) :
  rank_index(0), select_index(0)
{
    std::fread(&(this->size), sizeof(this->size), 1, file);
    std::fread(&(this->items), sizeof(this->items), 1, file);
    std::fread(&(this->number_of_blocks), sizeof(this->number_of_blocks), 1, file);
    std::fread(&(this->block_size), sizeof(this->block_size), 1, file);

  this->array = new usint[this->block_size * this->number_of_blocks];
  std::fread(this->array, sizeof(usint), this->block_size * this->number_of_blocks, file);
  this->buffer = new FastBitBuffer(this->array, this->block_size);

  this->integer_bits = length(this->size);
  this->samples = new FastBitBuffer(file, 2 * (this->number_of_blocks + 1), this->integer_bits);

  this->indexForRank();
  this->indexForSelect();
}

BitVector::BitVector(VectorEncoder& encoder, usint universe_size) :
  size(universe_size), items(encoder.items),
  block_size(encoder.block_size),
  number_of_blocks(encoder.blocks),
  rank_index(0), select_index(0)
{
  this->array = new usint[this->block_size * this->number_of_blocks];
  this->buffer = new FastBitBuffer(this->array, this->block_size);

  this->integer_bits = length(this->size);
  this->samples = new FastBitBuffer(2 * (this->number_of_blocks + 1), this->integer_bits);

  // Copy & linearize the array.
  usint pos = 0;
  for(std::list<usint*>::iterator iter = encoder.array_blocks.begin(); iter != encoder.array_blocks.end(); iter++)
  {
    memcpy(this->array + pos, *iter, encoder.superblock_bytes);
    pos += encoder.block_size * encoder.blocks_in_superblock;
  }
  memcpy(this->array + pos, encoder.array, encoder.current_blocks * encoder.block_size * sizeof(usint));

  // Compress the samples.
  for(std::list<usint*>::iterator iter = encoder.sample_blocks.begin(); iter != encoder.sample_blocks.end(); iter++)
  {
    usint* buf = *iter;
    for(usint i = 0; i < 2 * encoder.samples_in_superblock; i++)
    {
      this->samples->writeItem(buf[i]);
    }
  }
  for(usint i = 0; i < 2 * encoder.current_samples; i++)
  {
    this->samples->writeItem(encoder.samples[i]);
  }
  this->samples->writeItem(this->items);
  this->samples->writeItem(this->size);

  this->indexForRank();
  this->indexForSelect();
}

BitVector::~BitVector()
{
  delete[] this->array;
  delete   this->buffer;
  delete   this->samples;
  delete   this->rank_index;
  delete   this->select_index;
}

void
BitVector::writeTo(std::ofstream& file)
{
  file.write((char*)&(this->size), sizeof(this->size));
  file.write((char*)&(this->items), sizeof(this->items));
  file.write((char*)&(this->number_of_blocks), sizeof(this->number_of_blocks));
  file.write((char*)&(this->block_size), sizeof(this->block_size));
  file.write((char*)(this->array), this->block_size * this->number_of_blocks * sizeof(usint));
  this->samples->writeBuffer(file, false);
}

void
BitVector::writeTo(FILE* file)
{
    std::fwrite(&(this->size), sizeof(this->size), 1, file);
    std::fwrite(&(this->items), sizeof(this->items), 1, file);
    std::fwrite(&(this->number_of_blocks), sizeof(this->number_of_blocks), 1, file);
    std::fwrite(&(this->block_size), sizeof(this->block_size), 1, file);
    std::fwrite(this->array, sizeof(usint), this->block_size * this->number_of_blocks, file);
    this->samples->writeBuffer(file, false);
}


//--------------------------------------------------------------------------

usint
BitVector::reportSize()
{
  usint bytes = this->buffer->reportSize();
  bytes += this->block_size * this->number_of_blocks * sizeof(usint);
  bytes += this->samples->reportSize();
  bytes += this->rank_index->reportSize();
  bytes += this->select_index->reportSize();
  return bytes;
}

//--------------------------------------------------------------------------

usint
BitVector::sampleForIndex(usint index)
{
  usint low = this->select_index->readItem(index / this->select_rate);
  usint high = this->select_index->readItem(index / this->select_rate + 1);

  this->samples->goToItem(2 * low + 2);
  for(; low < high; low++)
  {
    if(this->samples->readItem() > index) { return low; }
    this->samples->skipItem();
  }

  return low;
}

usint
BitVector::sampleForValue(usint value)
{
  usint low = this->rank_index->readItem(value / this->rank_rate);
  usint high = this->rank_index->readItem(value / this->rank_rate + 1);

  this->samples->goToItem(2 * low + 3);
  for(; low < high; low++)
  {
    if(this->samples->readItem() > value) { return low; }
    this->samples->skipItem();
  }

  return low;
}

//--------------------------------------------------------------------------

void
BitVector::indexForRank()
{
  delete this->rank_index;

  usint value_samples = (this->number_of_blocks + BitVector::INDEX_RATE - 1) / BitVector::INDEX_RATE;
  this->rank_rate = (this->size + value_samples - 1) / value_samples;
  value_samples = (this->size + this->rank_rate - 1) / this->rank_rate + 1;
  this->rank_index = new FastBitBuffer(value_samples, length(this->number_of_blocks - 1));

  usint current = 0, pointer = 0;
  this->samples->goToItem(2);
  while(this->samples->hasNextItem())
  {
    this->samples->skipItem();
    usint limit = this->samples->readItem();
    while(current < limit)
    {
      this->rank_index->writeItem(pointer);
      current += this->rank_rate;
    }
    pointer++;
  }
  this->rank_index->writeItem(this->number_of_blocks - 1);
}

void
BitVector::indexForSelect()
{
  delete this->select_index;

  usint index_samples = (this->number_of_blocks + BitVector::INDEX_RATE - 1) / BitVector::INDEX_RATE;
  this->select_rate = (this->items + index_samples - 1) / index_samples;
  index_samples = (this->items + this->select_rate - 1) / this->select_rate + 1;
  this->select_index = new FastBitBuffer(index_samples, length(this->number_of_blocks - 1));

  usint current = 0, pointer = 0;
  this->samples->goToItem(2);
  while(this->samples->hasNextItem())
  {
    usint limit = this->samples->readItem();
    this->samples->skipItem();
    while(current < limit)
    {
      this->select_index->writeItem(pointer);
      current += this->select_rate;
    }
    pointer++;
  }
  this->select_index->writeItem(this->number_of_blocks - 1);
}

//--------------------------------------------------------------------------

VectorEncoder::VectorEncoder(usint block_bytes, usint superblock_size) :
  size(0), items(0), blocks(0),
  block_size(BYTES_TO_WORDS(block_bytes)),
  superblock_bytes(superblock_size)
{
  this->array = new usint[this->superblock_bytes / sizeof(usint)];
  memset(this->array, 0, this->superblock_bytes);
  this->blocks_in_superblock = this->superblock_bytes / (sizeof(usint) * this->block_size);
  this->current_blocks = 0;

  this->samples = new usint[this->superblock_bytes / sizeof(usint)];
  this->samples_in_superblock = this->superblock_bytes / (2 * sizeof(usint));
  this->current_samples = 0;

  this->buffer = new FastBitBuffer(this->array, this->block_size);
}

VectorEncoder::~VectorEncoder()
{
  delete[] this->array;

  delete this->buffer;
  for(std::list<usint*>::iterator iter = this->array_blocks.begin(); iter != this->array_blocks.end(); iter++)
  {
    delete[] *iter;
  }

  delete[] this->samples;
  for(std::list<usint*>::iterator iter = this->sample_blocks.begin(); iter != this->sample_blocks.end(); iter++)
  {
    delete[] *iter;
  }
}

void
VectorEncoder::addNewBlock()
{
  this->blocks++;
  this->current_blocks++;
  this->current_samples++;

  // Do we need a new superblock for the block?
  if(this->current_blocks > this->blocks_in_superblock)
  {
    this->array_blocks.push_back(this->array);
    this->array = new usint[this->superblock_bytes / sizeof(usint)];
    memset(this->array, 0, this->superblock_bytes);
    this->current_blocks = 1;
  }
  this->buffer->moveBuffer(this->array + (this->block_size * (this->current_blocks - 1)));

  // Do we need a new superblock for the sample?
  if(this->current_samples > this->samples_in_superblock)
  {
    this->sample_blocks.push_back(this->samples);
    this->samples = new usint[this->superblock_bytes / sizeof(usint)];
    this->current_samples = 1;
  }
  this->samples[2 * this->current_samples - 2] = this->items - 1;
  this->samples[2 * this->current_samples - 1] = this->size - 1;
}

void
VectorEncoder::setFirstBit(usint value)
{
  this->samples[0] = 0;
  this->samples[1] = value;

  this->size = value + 1;
  this->items = 1;
  this->blocks = 1;

  this->current_blocks = 1;
  this->current_samples = 1;
}


} // namespace CSA
