#ifndef BITBUFFER_H
#define BITBUFFER_H

#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <cstring>  // defines std::memset, added by Kim

#include "../misc/definitions.h"


namespace CSA
{


template<class Data>
class GenericBitBuffer
{
  public:
    GenericBitBuffer(usint words);
    GenericBitBuffer(usint _items, usint item_size);
    GenericBitBuffer(std::ifstream& file, usint words);
    GenericBitBuffer(std::ifstream& file, usint _items, usint item_size);
    GenericBitBuffer(std::FILE* file, usint _items, usint item_size);
    GenericBitBuffer(Data* buffer, usint words);
    GenericBitBuffer(Data* buffer, usint _items, usint item_size);
    ~GenericBitBuffer();

    void writeBuffer(std::ofstream& file, bool erase = true);
    void writeBuffer(std::FILE * file, bool erase = true);
    void readBuffer(std::ifstream& file, usint words, bool erase = true);
    void setBuffer(Data* buffer, usint words);

    // We assume we are already using a buffer not owned by this.
    inline void moveBuffer(Data* buffer)
    {
      this->data = buffer;
      this->pos = 0;
      this->bits = DATA_BITS;
      this->current = 0;
    }

//--------------------------------------------------------------------------

    inline void reset(bool erase)
    {
      this->pos = 0;
      this->bits = DATA_BITS;
      this->current = 0;
      if(erase)
      {
        memset(this->data, 0, this->size * sizeof(Data));
      }
    }

    inline void skipBits(usint count)
    {
      if(count < this->bits)
      {
        this->bits -= count;
        return;
      }

      count -= this->bits;
      this->pos += 1 + count / DATA_BITS;
      this->bits = DATA_BITS - count % DATA_BITS;
    }

    inline void rewind(usint count)
    {
      this->bits += count;
      if(this->bits > DATA_BITS)
      {
        usint back = (this->bits - 1) / DATA_BITS;
        this->pos -= back;
        this->bits -= back * DATA_BITS;
      }
    }

//--------------------------------------------------------------------------

    inline usint bitsLeft()
    {
      return this->bits + DATA_BITS * (this->size - this->pos - 1);
    }

    inline void writeBits(usint value, usint count)
    {
      while(count >= this->bits)
      {
        count -= this->bits;
        this->data[this->pos] |= GET(LOWER(value, count), this->bits);
        this->pos++; this->bits = DATA_BITS;
      }
      if(count > 0)
      {
        this->bits -= count;
        this->data[this->pos] |= HIGHER(GET(value, count), this->bits);
      }
    }

    // Returns nonzero if bit is 1
    inline usint readBit()
    {
      this->bits--;
      usint bit = this->data[this->pos] & ((usint)1 << this->bits);

      if(this->bits == 0) { this->pos++; this->bits = DATA_BITS; }

      return bit;
    }

    inline usint readBits(usint count)
    {
      usint value = 0;

      while(count >= this->bits)
      {
        count -= this->bits;
        value |= HIGHER(GET(this->data[this->pos], this->bits), count);
        this->pos++; this->bits = DATA_BITS;
      }
      if(count > 0)
      {
        this->bits -= count;
        value |= GET(LOWER(this->data[this->pos], this->bits), count);
      }

      return value;
    }

//--------------------------------------------------------------------------

    /*
      These operations work on fixed-size items.
    */

    inline usint getItemSize()
    {
      return this->item_bits;
    }

/*
    inline void setItemSize(usint new_size)
    {
      if(new_size == 0) { return; }
      this->item_bits = new_size;
    }
*/

    inline void goToItem(usint item)
    {
      usint b = item * this->item_bits;
      this->pos = b / DATA_BITS;
      this->bits = DATA_BITS - b % DATA_BITS;
      this->current = item;
    }

    inline usint readItem()
    {
      this->current++;
      return this->readBits(this->item_bits);
    }

    inline usint readItem(usint item)
    {
      this->goToItem(item);
      return this->readItem();
    }

    inline usint readFirstItem()
    {
      return this->readItem(0);
    }

    inline bool hasNextItem()
    {
      return (this->current < this->items);
    }

    inline void writeItem(usint item)
    {
      this->writeBits(item, this->item_bits);
      this->current++;
    }

    inline void skipItem()
    {
      this->skipBits(this->item_bits);
      this->current++;
    }

//--------------------------------------------------------------------------

    /*
      Delta coding for positive integers
    */

    inline bool canDeltaCode(usint value)
    {
      return this->deltaCodeLength(value) <= this->bitsLeft();
    }

    inline usint deltaCodeLength(usint value)
    {
      usint len = length(value);
      usint llen = length(len);
      return (len + llen + llen - 2);
    }

    // This version returns false if there is no space left for the encoding.
    inline bool writeDeltaCode(usint value)
    {
      usint len = length(value);
      usint llen = length(len);

      if(len + llen + llen - 2 > this->bitsLeft()) { return false; }

      // this->writeBits(0, llen - 1); // Now included in the next writeBits()
      this->writeBits(len, llen + llen - 1);
      this->writeBits(value, len - 1);
      return true;
    }

    // This version assumes the code fits into the buffer.
    inline void writeDeltaCodeDirect(usint value)
    {
      usint len = length(value);
      usint llen = length(len);

      // this->writeBits(0, llen - 1); // Now included in the next writeBits()
      this->writeBits(len, llen + llen - 1);
      this->writeBits(value, len - 1);
    }

    // We assume the code fits into usint:
    //  32-bit:  value < 2^24
    //  64-bit:  value < 2^54
    inline void writeDeltaCodeFast(usint value)
    {
      usint len = length(value);

      value ^= ((usint)1 << (len - 1));
      this->writeBits((len << (len - 1)) | value, len + 2 * length(len) - 2);
    }

    inline usint readDeltaCode()
    {
      usint len = 0;
      while(this->readBit() == 0) { len++; }

      usint temp = (((usint)1 << len) | this->readBits(len)) - 1;
      temp = ((usint)1 << temp) | this->readBits(temp);
      return temp;
    }

//--------------------------------------------------------------------------

    /*
      Gamma coding for positive integers
    */

    inline bool canGammaCode(usint value)
    {
      return this->gammaCodeLength(value) <= this->bitsLeft();
    }

    inline usint gammaCodeLength(usint value)
    {
      return 2 * length(value) - 1;
    }

    // This version returns false if there is no space left for the encoding.
    inline bool writeGammaCode(usint value)
    {
      usint len = length(value);

      if(len > this->bitsLeft()) { return false; }

      this->writeBits(0, len - 1);
      this->writeBits(value, len);
      return true;
    }

    // This version assumes the code fits into the buffer.
    inline void writeGammaCodeDirect(usint value)
    {
      usint len = length(value);

      this->writeBits(0, len - 1);
      this->writeBits(value, len);
    }

    // We assume the code fits into usint:
    //  32-bit:  value < 2^16
    //  64-bit:  value < 2^32
    inline void writeGammaCodeFast(usint value)
    {
      this->writeBits(value, this->gammaCodeLength(value));
    }

    inline usint readGammaCode()
    {
      usint len = 1;
      while(this->readBit() == 0) { len++; }
      return ((usint)1 << len) | this->readBits(len);
    }

//--------------------------------------------------------------------------

    usint reportSize();

//--------------------------------------------------------------------------

  private:
    Data* data;
    usint size, pos, bits;
    usint item_bits, items, current;
    bool  free_buffer;

    const static usint DATA_BITS = sizeof(Data) * CHAR_BIT;

    inline static usint bitsToData(usint _bits) { return (_bits + DATA_BITS - 1) / DATA_BITS; }
};


typedef GenericBitBuffer<uchar> BitBuffer;
typedef GenericBitBuffer<usint> FastBitBuffer;


//--------------------------------------------------------------------------


template<class Data>
GenericBitBuffer<Data>::GenericBitBuffer(usint words) :
  size(words),
  pos(0),
  bits(DATA_BITS),
  item_bits(1),
  items(0),
  current(0),
  free_buffer(true)
{
  this->data = new Data[words];
  memset(this->data, 0, words * sizeof(Data));
}

template<class Data>
GenericBitBuffer<Data>::GenericBitBuffer(usint _items, usint item_size) :
  pos(0),
  bits(DATA_BITS),
  item_bits(item_size),
  items(_items),
  current(0),
  free_buffer(true)
{
  this->size = bitsToData(items * item_size);
  this->data = new Data[this->size];
  memset(this->data, 0, this->size * sizeof(Data));
}

template<class Data>
GenericBitBuffer<Data>::GenericBitBuffer(std::ifstream& file, usint words) :
  size(words),
  pos(0),
  bits(DATA_BITS),
  item_bits(1),
  items(0),
  current(0),
  free_buffer(true)
{
  this->data = new Data[words];
  memset(this->data, 0, words * sizeof(Data));
  file.read((char*)this->data, words * sizeof(Data));
}

template<class Data>
GenericBitBuffer<Data>::GenericBitBuffer(std::ifstream& file, usint _items, usint item_size) :
  size(BITS_TO_WORDS(_items * item_size)),
  pos(0),
  bits(DATA_BITS),
  item_bits(item_size),
  items(_items),
  current(0),
  free_buffer(true)
{
  this->data = new Data[this->size];
  memset(this->data, 0, this->size * sizeof(Data));
  file.read((char*)this->data, this->size * sizeof(Data));
}

template<class Data>
GenericBitBuffer<Data>::GenericBitBuffer(std::FILE* file, usint _items, usint item_size) :
  size(BITS_TO_WORDS(_items * item_size)),
  pos(0),
  bits(DATA_BITS),
  item_bits(item_size),
  items(_items),
  current(0),
  free_buffer(true)
{
  this->data = new Data[this->size];
  memset(this->data, 0, this->size * sizeof(Data));
  std::fread(this->data, sizeof(Data), this->size, file);
}

template<class Data>
GenericBitBuffer<Data>::GenericBitBuffer(Data* buffer, usint words) :
  size(words),
  pos(0),
  bits(DATA_BITS),
  item_bits(1),
  items(0),
  current(0),
  free_buffer(false)
{
  this->data = buffer;
}

template<class Data>
GenericBitBuffer<Data>::GenericBitBuffer(Data* buffer, usint _items, usint item_size) :
  size(BITS_TO_WORDS(_items * item_size)),
  pos(0),
  bits(DATA_BITS),
  item_bits(item_size),
  items(_items),
  current(0),
  free_buffer(false)
{
  this->data = buffer;
}

template<class Data>
GenericBitBuffer<Data>::~GenericBitBuffer()
{
  if(this->free_buffer)
  {
    delete[] this->data;
  }
}

//--------------------------------------------------------------------------

template<class Data>
void
GenericBitBuffer<Data>::writeBuffer(std::ofstream& file, bool erase)
{
  file.write((char*)this->data, this->size * sizeof(Data));
  this->reset(erase);
}

template<class Data>
void
GenericBitBuffer<Data>::writeBuffer(std::FILE* file, bool erase)
{
    std::fwrite(this->data, sizeof(Data), this->size, file);
    this->reset(erase);
}

template<class Data>
void
GenericBitBuffer<Data>::readBuffer(std::ifstream& file, usint words, bool erase)
{
  if(words > this->size || !(this->free_buffer))
  {
    if(this->free_buffer)
    {
      delete[] this->data;
    }
    this->size = words;
    this->data = new Data[words];
    this->free_buffer = true;
  }

  this->reset(erase);
  file.read((char*)this->data, words * sizeof(Data));
}

template<class Data>
void
GenericBitBuffer<Data>::setBuffer(Data* buffer, usint words)
{
  if(this->free_buffer)
  {
    delete[] this->data;
    this->free_buffer = false;
  }

  this->data = buffer;
  this->size = words;
  this->reset(false);
}

//--------------------------------------------------------------------------

template<class Data>
usint
GenericBitBuffer<Data>::reportSize()
{
  usint bytes = sizeof(*this);
  if(this->free_buffer) { bytes += this->size * sizeof(Data); }
  return bytes;
}

//--------------------------------------------------------------------------


} // namespace CSA


#endif // BITBUFFER_H
