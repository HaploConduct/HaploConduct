#ifndef RLCSA_H
#define RLCSA_H

#include <fstream>
#include <cstring> // defines std::memset, added by Kim

#include "bits/vectors.h"
#include "sasamples.h"
#include "misc/parameters.h"


namespace CSA
{


const std::string SA_EXTENSION = ".sa";
const std::string PSI_EXTENSION = ".psi";
const std::string ARRAY_EXTENSION = ".rlcsa.array";
const std::string SA_SAMPLES_EXTENSION = ".rlcsa.sa_samples";
const std::string PARAMETERS_EXTENSION = ".rlcsa.parameters";


const parameter_type RLCSA_BLOCK_SIZE  = parameter_type("RLCSA_BLOCK_SIZE", 32);
const parameter_type SAMPLE_RATE       = parameter_type("SAMPLE_RATE", 512);
const parameter_type SUPPORT_LOCATE    = parameter_type("SUPPORT_LOCATE", 0);
const parameter_type SUPPORT_DISPLAY   = parameter_type("SUPPORT_DISPLAY", 0);


struct LocateItem
{
  usint value;
  usint offset;
  bool found;
};



class RLCSA
{
  public:
    const static usint ENDPOINT_BLOCK_SIZE = 16;

    RLCSA(const std::string& base_name, bool print = false);

    /*
      If multiple_sequences is true, each 0 is treated as a end of sequence marker.
      There must be nonzero characters between the 0s. data[bytes - 1] must also be 0.
    */ 
    RLCSA(uchar* data, usint bytes, usint block_size, usint sa_sample_rate, bool multiple_sequences, bool delete_data);

    // Destroys contents of index and increment.
    RLCSA(RLCSA& index, RLCSA& increment, usint* positions, usint block_size);
    ~RLCSA();

    void writeTo(const std::string& base_name);

    bool isOk();

    // Returns the closed interval of suffix array containing the matches.
    pair_type count(const std::string& pattern);

    void reportPositions(uchar* data, usint length, usint* positions);

    // Returns SA[range]. User must free the buffer. Latter version uses buffer provided by the user.
    LocateItem* locate(pair_type range);
    LocateItem* locate(pair_type range, LocateItem* data);

    // Returns SA[index].
    usint locate(usint index);

    // Returns Ti[range]. User must free the buffer. Latter version uses buffer provided by the user.
    uchar* display(usint sequence, pair_type range);
    uchar* display(usint sequence, pair_type range, uchar* data);

    // Writes the Psi array into given file. End of sequence markers are not written.
    void decompressInto(std::ofstream& psi_file);

    // Returns the BWT of the collection including end of sequence markers.
    uchar* readBWT();

    // Return value includes the implicit end of sequence markers. To get suffix array indexes,
    // subtract getNumberOfSequences() from the value.
    usint psi(usint index);
    pair_type psi(usint index, usint max_length); // This version returns a run.

    inline bool supportsLocate()       { return this->support_locate; }
    inline bool supportsDisplay()      { return this->support_display; }
    inline usint getSize()              { return this->data_size; }
    inline usint getNumberOfSequences() { return this->number_of_sequences; }

    pair_type getSequence(usint number);

    // Returns the size of the data structure.
    usint reportSize(bool print = false);

  private:
    bool ok;

    RLEVector* array[CHARS];
    SASamples* sa_samples;

    pair_type index_ranges[CHARS];
    usint data_size;

    usint text_chars[CHARS];  // which characters are present in the text
    usint chars;

    usint index_pointers[CHARS]; // which of the above is at i * index_rate
    usint index_rate;

    bool support_locate, support_display;
    usint sample_rate;

    usint number_of_sequences;
    DeltaVector* end_points;

    void locateUnsafe(pair_type range, LocateItem* data);
    bool processRun(pair_type run, LocateItem* data);
    void displayUnsafe(pair_type range, uchar* data);

    inline usint getCharacter(usint pos)
    {
      usint* curr = &(this->text_chars[this->index_pointers[pos / this->index_rate]]);
      while(pos > this->index_ranges[*curr].second) { curr++; }
      return *curr;
    }

    void buildCharIndexes(usint* distribution);
};


// Returns the total number of characters.
usint buildRanges(usint* distribution, pair_type* index_ranges);


} // namespace CSA


#endif // RLCSA_H
