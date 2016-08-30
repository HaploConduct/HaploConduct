#ifndef SASAMPLES_H
#define SASAMPLES_H

#include <fstream>

#include "misc/definitions.h"
#include "bits/bitbuffer.h"
#include "bits/deltavector.h"


namespace CSA
{


class SASamples
{
  public:
    const static usint BLOCK_SIZE = 16;

    SASamples(std::ifstream& sample_file, usint sample_rate);
    SASamples(usint* array, usint data_size, usint sample_rate);
    ~SASamples();

    // Destroys contents of index and increment.
    // We assume index and increment have same sample rate.
    SASamples(SASamples& index, SASamples& increment, usint* positions, usint number_of_positions);

    void writeTo(std::ofstream& sample_file);

    // Returns i such that SA[i] = value.
    // If SA[i] is not sampled, returns the next sampled value. (Don't try!)
    // Value is actual 0-based suffix array value.
    // Returns size if value is too large.
    usint inverseSA(usint value);

    // Returns the value of ith sample in suffix array order.
    inline usint getSample(usint i)
    {
      return std::min(this->samples->readItem(i) * this->rate, this->size - 1);
    }

    // Returns (ind, sample number) where ind >= index or (size, ???).
    pair_type getFirstSampleAfter(usint index);

    inline usint getSampleRate() { return this->rate; }
    inline usint getNumberOfSamples() { return this->items; }

    usint reportSize();

  private:
    usint integer_bits;
    usint rate, size, items;

    DeltaVector* indexes;

    FastBitBuffer* samples;
    FastBitBuffer* inverse_samples;

    void buildInverseSamples();

    // Note: contents of original samples are deleted.
    void mergeSamples(SASamples& index, SASamples& increment, usint* positions, usint n);
};


} // namespace CSA


#endif // SASAMPLES_H
