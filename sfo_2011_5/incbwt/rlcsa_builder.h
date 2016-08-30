#ifndef RLCSA_BUILDER_H
#define RLCSA_BUILDER_H

#include <cstdlib>
#include "rlcsa.h"


namespace CSA
{


class RLCSABuilder
{
  public:
    RLCSABuilder(usint _block_size, usint _sample_rate, usint _buffer_size);
    ~RLCSABuilder();

    void insertSequence(char* sequence, usint length, bool delete_sequence);

    // User must free the index. Builder no longer contains it.
    RLCSA* getRLCSA();

    // User must free the BWT. length becomes the length of BWT.
    char* getBWT(usint& length);

    bool isOk();

    // These times are not reset with the rest of the builder.
    double getBuildTime();
    double getSearchTime();
    double getMergeTime();

  private:
    RLCSA* index;

    usint block_size;
    usint sample_rate;
    usint buffer_size;

    uchar* buffer;
    usint chars;

    bool ok;

    clock_t build_time;
    clock_t search_time;
    clock_t merge_time;

    void flush();
    void addRLCSA(RLCSA* increment, uchar* sequence, usint length, bool delete_sequence);
    void reset();
};


} // namespace CSA


#endif // RLCSA_BUILDER_H
