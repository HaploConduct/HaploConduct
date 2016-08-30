#include <algorithm>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <iostream>

#include "rlcsa.h"
#include "misc/utils.h"
#include "qsufsort/qsufsort.h"


namespace CSA
{


RLCSA::RLCSA(const std::string& base_name, bool print) :
  ok(false),
  sa_samples(0),
  end_points(0)
{
  for(usint c = 0; c < CHARS; c++) { this->array[c] = 0; }

  std::string array_name = base_name + ARRAY_EXTENSION;
  std::ifstream array_file(array_name.c_str(), std::ios_base::binary);
  if(!array_file)
  {
    std::cerr << "RLCSA: Error opening Psi array file!" << std::endl;
    return;
  }

  usint distribution[CHARS];
  array_file.read((char*)distribution, CHARS * sizeof(usint));
  this->buildCharIndexes(distribution);

  Parameters parameters;
  parameters.read(base_name + PARAMETERS_EXTENSION);
  for(usint c = 0; c < CHARS; c++)
  {
    if(!isEmpty(this->index_ranges[c])) { this->array[c] = new RLEVector(array_file); }
  }

  this->end_points = new DeltaVector(array_file);
  this->number_of_sequences = this->end_points->getNumberOfItems();

  array_file.read((char*)&(this->sample_rate), sizeof(this->sample_rate));
  array_file.close();

  this->support_locate = parameters.get(SUPPORT_LOCATE);
  this->support_display = parameters.get(SUPPORT_DISPLAY);

  if(this->support_locate || this->support_display)
  {
    std::string sa_sample_name = base_name + SA_SAMPLES_EXTENSION;
    std::ifstream sa_sample_file(sa_sample_name.c_str(), std::ios_base::binary);
    if(!sa_sample_file)
    {
      std::cerr << "RLCSA: Error opening suffix array sample file!" << std::endl;
      return;
    }
    this->sa_samples = new SASamples(sa_sample_file, this->sample_rate);
    sa_sample_file.close();
  }

  if(print) { parameters.print(); }

  this->ok = true;
}

RLCSA::RLCSA(uchar* data, usint bytes, usint block_size, usint sa_sample_rate, bool multiple_sequences, bool delete_data) :
  ok(false),
  sa_samples(0),
  sample_rate(sa_sample_rate),
  end_points(0)
{
  for(usint c = 0; c < CHARS; c++) { this->array[c] = 0; }

  if(!data || bytes == 0)
  {
    std::cerr << "RLCSA: No input data given!" << std::endl;
    if(delete_data) { delete[] data; }
    return;
  }
  if(block_size < 2 * sizeof(usint) || block_size % sizeof(usint) != 0)
  {
    std::cerr << "RLCSA: Block size must be a multiple of " << sizeof(usint) << " bytes!" << std::endl;
    if(delete_data) { delete[] data; }
    return;
  }


  // Do we store SA samples?
  if(this->sample_rate > 0)
  {
    this->support_locate = this->support_display = true;
  }
  else
  {
    this->support_locate = this->support_display = false;
    this->sample_rate = 1;
  }


  // Determine the number of sequences and mark their end points.
  if(multiple_sequences)
  {
    DeltaEncoder endings(RLCSA::ENDPOINT_BLOCK_SIZE);
    this->number_of_sequences = 0;
    usint marker = 0;
    usint padding = 0, chars_encountered = 0;

    for(usint i = 0; i < bytes; i++)
    {
      if(data[i] == 0)
      {
        if(i == marker) { break; }  // Empty sequence.
        this->number_of_sequences++;
        marker = i + 1;

        usint pos = chars_encountered + padding - 1;
        endings.setBit(pos);
        padding = ((pos + this->sample_rate) / this->sample_rate) * this->sample_rate - chars_encountered;
      }
      else
      {
        chars_encountered++;
      }
    }

    if(this->number_of_sequences == 0 || marker != bytes)
    {
      std::cerr << "RLCSA: Collection must consist of 0-terminated nonempty sequences!" << std::endl;
      if(delete_data) { delete[] data; }
      return;
    }
    this->end_points = new DeltaVector(endings, chars_encountered + padding);
  }
  else
  {
    this->number_of_sequences = 1;
    DeltaEncoder endings(RLCSA::ENDPOINT_BLOCK_SIZE, RLCSA::ENDPOINT_BLOCK_SIZE);
    usint pos = bytes - 1;
    endings.setBit(pos);
    usint padding = ((pos + this->sample_rate) / this->sample_rate) * this->sample_rate - bytes;
    this->end_points = new DeltaVector(endings, bytes + padding);
  }


  // Build character tables etc.
  usint distribution[CHARS];
  sint low = CHARS, high = 0;
  for(usint c = 0; c < CHARS; c++) { distribution[c] = 0; }
  for(usint i = 0; i < bytes; i++)
  {
    if(data[i] < low) { low = data[i]; }
    if(data[i] > high) { high = data[i]; }
    distribution[(usint)data[i]]++;
  }
  if(multiple_sequences)
  {
    distribution[0] = 0;
    low = 0;
    high = high + this->number_of_sequences;
  }
  this->buildCharIndexes(distribution);


  // Build suffix array.
  sint* inverse = new sint[bytes + 1];
  if(multiple_sequences)
  {
    sint zeros = 0;
    for(usint i = 0; i < bytes; i++)
    {
      if(data[i] == 0)
      {
        inverse[i] = zeros;
        zeros++;
      }
      else
      {
        inverse[i] = (sint)data[i] + this->number_of_sequences;
      }
    }
  }
  else
  {
    for(usint i = 0; i < bytes; i++) { inverse[i] = (sint)data[i]; }
  }
  if(delete_data) { delete[] data; }
  sint* sa = new sint[bytes + 1];
  suffixsort(inverse, sa, bytes, high + 1, low);


  // Sample SA.
  usint incr = (multiple_sequences ? this->number_of_sequences + 1 : 1);  // No e_of_s markers in SA.
  if(this->support_locate || this->support_display)
  {
    if(multiple_sequences)
    {
      std::cout << "We shouldn't be here!" << std::endl;
      // Real SA starts at sa + incr.
      // for each sequence
      //   find starting position
      //   sample relative to starting position
      // sort samples to SA order
    }
    else
    {
      this->sa_samples = new SASamples((usint*)(sa + incr), this->data_size, this->sample_rate);
    }
  }


  // Build Psi.
  for(usint i = 0; i <= bytes; i++)
  {
    sa[i] = inverse[(sa[i] + 1) % (bytes + 1)];
  }
  delete[] inverse;


  // Build RLCSA.
  usint decr = (multiple_sequences ? 1 : 0);  // Ignore the last e_of_s marker if multiple sequences.
  for(usint c = 0; c < CHARS; c++)
  {
    if(distribution[c] == 0) { this->array[c] = 0; continue; }

    usint* curr = (usint*)(sa + index_ranges[c].first + incr);
    usint* limit = (usint*)(sa + index_ranges[c].second + incr + 1);
    RLEEncoder encoder(block_size);
    pair_type run(*curr, 1);
    curr++;

    for(; curr < limit; curr++)
    {
      if(*curr == run.first + run.second) { run.second++; }
      else
      {
        encoder.setRun(run.first - decr, run.second);
        run = pair_type(*curr, 1);
      }
    }
    encoder.setRun(run.first - decr, run.second);

    this->array[c] = new RLEVector(encoder, this->data_size + incr - decr);
  }
  delete[] sa;


  this->ok = true;
}

RLCSA::RLCSA(RLCSA& index, RLCSA& increment, usint* positions, usint block_size) :
  ok(false),
  sa_samples(0),
  end_points(0)
{
  for(usint c = 0; c < CHARS; c++) { this->array[c] = 0; }

  if(!index.isOk() || !increment.isOk())
  {
    return; // Fail silently. Actual error has already been reported.
  }
  if(positions == 0)
  {
    std::cerr << "RLCSA: Positions for insertions not available!" << std::endl;
    return;
  }
  if(index.sample_rate != increment.sample_rate)
  {
    std::cerr << "RLCSA: Cannot combine indexes with different sample rates!" << std::endl;
    std::cout << "Index: " << index.sample_rate << ", increment: " << increment.sample_rate << std::endl;
    return;
  }

  if(index.sa_samples == 0 || increment.sa_samples == 0)
  {
    this->support_locate = this->support_display = false;
  }
  else
  {
    this->support_locate = this->support_display = true;
  }


  // Build character tables etc.
  usint distribution[CHARS];
  for(usint c = 0; c < CHARS; c++)
  {
    distribution[c] = length(index.index_ranges[c]) + length(increment.index_ranges[c]);
  }
  this->buildCharIndexes(distribution);
  this->sample_rate = index.sample_rate;


  // Combine end points of sequences.
  this->number_of_sequences = index.number_of_sequences + increment.number_of_sequences;
  DeltaEncoder* endings = new DeltaEncoder(RLCSA::ENDPOINT_BLOCK_SIZE);

  endings->setBit(index.end_points->select(0));
  for(usint i = 1; i < index.number_of_sequences; i++)
  {
    endings->setBit(index.end_points->selectNext());
  }
  usint sum = index.end_points->getSize();
  delete index.end_points; index.end_points = 0;

  endings->setBit(sum + increment.end_points->select(0));
  for(usint i = 1; i < increment.number_of_sequences; i++)
  {
    endings->setBit(sum + increment.end_points->selectNext());
  }
  sum += increment.end_points->getSize();
  delete increment.end_points; increment.end_points = 0;

  this->end_points = new DeltaVector(*endings, sum);
  delete endings;


  // Combine Psi.
  usint psi_size = this->data_size + this->number_of_sequences;
  for(usint c = 0; c < CHARS; c++)
  {
    if(distribution[c] == 0) { this->array[c] = 0; continue; }
    this->array[c] = mergeVectors(index.array[c], increment.array[c], positions, increment.data_size + increment.number_of_sequences, psi_size, block_size);
    index.array[c] = 0;
    increment.array[c] = 0;

    if(this->array[c] == 0)
    {
      std::cerr << "RLCSA: Merge failed for vectors " << c << "!" << std::endl;
      return;
    }
  }


  // Combine suffix array samples.
  if(this->support_locate || this->support_display)
  {
    positions += increment.number_of_sequences;
    for(usint i = 0; i < increment.data_size; i++)
    {
      positions[i] -= this->number_of_sequences;
    }
    this->sa_samples = new SASamples(*(index.sa_samples), *(increment.sa_samples), positions, increment.data_size);
  }

  this->ok = true;
}

RLCSA::~RLCSA()
{
  for(usint c = 0; c < CHARS; c++) { delete this->array[c]; }
  delete this->sa_samples;
  delete this->end_points;
}

void
RLCSA::writeTo(const std::string& base_name)
{
  std::string array_name = base_name + ARRAY_EXTENSION;
  std::ofstream array_file(array_name.c_str(), std::ios_base::binary);
  if(!array_file)
  {
    std::cerr << "RLCSA: Error creating Psi array file!" << std::endl;
    return;
  }

  usint distribution[CHARS];
  for(usint c = 0; c < CHARS; c++)
  {
    distribution[c] = length(this->index_ranges[c]);
  }
  array_file.write((char*)distribution, CHARS * sizeof(usint));
  for(usint c = 0; c < CHARS; c++)
  {
    if(this->array[c] != 0)
    {
      this->array[c]->writeTo(array_file);
    }
  }

  this->end_points->writeTo(array_file);
  array_file.write((char*)&(this->sample_rate), sizeof(this->sample_rate));
  array_file.close();

  if(this->support_locate || this->support_display)
  {
    std::string sa_sample_name = base_name + SA_SAMPLES_EXTENSION;
    std::ofstream sa_sample_file(sa_sample_name.c_str(), std::ios_base::binary);
    if(!sa_sample_file)
    {
      std::cerr << "RLCSA: Error creating suffix array sample file!" << std::endl;
      return;
    }

    this->sa_samples->writeTo(sa_sample_file);
    sa_sample_file.close();
  }
}

bool
RLCSA::isOk()
{
  return this->ok;
}

//--------------------------------------------------------------------------

pair_type
RLCSA::count(const std::string& pattern)
{
  if(pattern.length() == 0) { return pair_type(0, this->data_size - 1); }

  pair_type index_range = this->index_ranges[(usint)*(pattern.rbegin())];
  if(isEmpty(index_range)) { return index_range; }
  index_range.first += this->number_of_sequences;
  index_range.second += this->number_of_sequences;

  for(std::string::const_reverse_iterator iter = ++pattern.rbegin(); iter != pattern.rend(); iter++)
  {
    RLEVector* vector = this->array[(usint)*iter];
    usint start = this->index_ranges[(usint)*iter].first;

    index_range.first = start + vector->rank(index_range.first, true) - 1 + this->number_of_sequences;
    index_range.second = start + vector->rank(index_range.second) - 1 + this->number_of_sequences;

    if(isEmpty(index_range)) { return index_range; }
  }

  // Suffix array indexes are 0-based.
  index_range.first -= this->number_of_sequences;
  index_range.second -= this->number_of_sequences;

  return index_range;
}

void
RLCSA::reportPositions(uchar* data, usint length, usint* positions)
{
  if(data == 0 || length == 0 || positions == 0) { return; }

  usint current = this->number_of_sequences - 1;
  positions[length] = current; // "immediately after current"
  for(sint i = (sint)(length - 1); i >= 0; i--)
  {
//     positions[i] = current; // "immediately after current"
    usint c = (usint)data[i];
    if(array[c])
    {
      current = this->index_ranges[c].first + this->array[c]->rank(current) - 1 + this->number_of_sequences;
    }
    else
    {
      if(c < this->text_chars[0]) // No previous characters either.
      {
        current = this->number_of_sequences - 1;
      }
      else
      {
        current = this->index_ranges[c].first - 1 + this->number_of_sequences;
      }
    }
    positions[i] = current; // "immediately after current"
  }
}

//--------------------------------------------------------------------------

LocateItem*
RLCSA::locate(pair_type range)
{
  if(!(this->support_locate) || isEmpty(range) || range.second >= this->data_size) { return 0; }

  LocateItem* data = new LocateItem[range.second + 1 - range.first];
  if(!data) { return 0; }
  this->locateUnsafe(range, data);

  return data;
}

LocateItem*
RLCSA::locate(pair_type range, LocateItem* data)
{
  if(!(this->support_locate) || isEmpty(range) || range.second >= this->data_size || data == 0) { return 0; }
  this->locateUnsafe(range, data);
  return data;
}

void
RLCSA::locateUnsafe(pair_type range, LocateItem* data)
{
  usint items = range.second + 1 - range.first;
  for(usint i = 0, j = range.first; i < items; i++, j++)
  {
    data[i].value = j + this->number_of_sequences;
    data[i].offset = 0;
    data[i].found = false;
  }

  bool found = false;
  while(!found)
  {
    found = true;
    pair_type run = EMPTY_PAIR;
    for(usint i = 0; i < items; i++)
    {
      if(data[i].found)
      {
        continue; // The run might continue after this.
      }
      else if(isEmpty(run))
      {
        run = pair_type(i, i);
      }
      else if(data[i].value - data[run.first].value == i - run.first)
      {
        run.second = i;
      }
      else
      {
        found &= this->processRun(run, data);
        run = pair_type(i, i);
      }
    }
    if(!isEmpty(run)) { found &= this->processRun(run, data); }
  }
}

bool
RLCSA::processRun(pair_type run, LocateItem* data)
{
  bool found = true;
  usint run_start = 0, run_left = 0;
  pair_type next_sample = pair_type(0, 0);

  for(usint i = run.first; i <= run.second; i++)
  {
    if(data[i].found)
    {
      if(run_left > 0) { run_left--; }
      continue;
    }
    if(data[i].value < this->number_of_sequences) // Implicit sample here.
    {
      data[i].value = this->end_points->select(data[i].value) + 1 - data[i].offset;
      data[i].found = true;
      if(run_left > 0) { run_left--; }
      continue;
    }
    if(next_sample.first < data[i].value) // Need another sample.
    {
      next_sample = this->sa_samples->getFirstSampleAfter(data[i].value - this->number_of_sequences);
      next_sample.first += this->number_of_sequences;
    }
    if(data[i].value < next_sample.first) // No sample found for current position.
    {
      if(run_left > 0)
      {
        data[i].value = data[run_start].value + i - run_start;
        run_left--;
      }
      else
      {
        pair_type value = this->psi(data[i].value - this->number_of_sequences, run.second - i);
        data[i].value = value.first;
        run_left = value.second;
        run_start = i;
      }
      data[i].offset++;
      found = false;
    }
    else  // Sampled position found.
    {
      data[i].value = this->sa_samples->getSample(next_sample.second) - data[i].offset;
      data[i].found = true;
      if(run_left > 0) { run_left--; }
    }
  }
  return found;
}

usint
RLCSA::locate(usint index)
{
  if(!(this->support_locate) || index >= this->data_size) { return this->data_size; }

  usint offset = 0;
  index += this->number_of_sequences;
  while(true)
  {
    if(index < this->number_of_sequences) // Implicit sample here
    {
      return this->end_points->select(index) + 1 - offset;
    }
    pair_type next_sample = this->sa_samples->getFirstSampleAfter(index - this->number_of_sequences);
    next_sample.first += this->number_of_sequences;
    if(next_sample.first == index)
    {
      return this->sa_samples->getSample(next_sample.second) - offset;
    }
    index = this->psi(index - this->number_of_sequences);
    offset++;
  }
}

//--------------------------------------------------------------------------

uchar*
RLCSA::display(usint sequence, pair_type range)
{
  if(!(this->support_display) || isEmpty(range)) { return 0; }

  pair_type seq_range = this->getSequence(sequence);
  if(isEmpty(seq_range)) { return 0; }

  range.first += seq_range.first; range.second += seq_range.first;
  if(range.second > seq_range.second) { return 0; }

  uchar* data = new uchar[range.second + 1 - range.first];
  if(!data) { return 0; }
  this->displayUnsafe(range, data);

  return data;
}

uchar*
RLCSA::display(usint sequence, pair_type range, uchar* data)
{
  if(!(this->support_display) || isEmpty(range) || data == 0) { return 0; }

  pair_type seq_range = this->getSequence(sequence);
  if(isEmpty(seq_range)) { return 0; }

  range.first += seq_range.first; range.second += seq_range.first;
  if(range.second > seq_range.second) { return 0; }

  this->displayUnsafe(range, data);
  return data;
}

void
RLCSA::displayUnsafe(pair_type range, uchar* data)
{
  usint i = range.first - range.first % this->sa_samples->getSampleRate();

  usint pos = this->sa_samples->inverseSA(i);

  for(; i < range.first; i++)
  {
    pos = this->psi(pos) - this->number_of_sequences;
  }
  for(; i <= range.second; i++)
  {
    data[i - range.first] = this->getCharacter(pos);
    pos = this->psi(pos) - this->number_of_sequences;
  }
}

//--------------------------------------------------------------------------

void
RLCSA::decompressInto(std::ofstream& psi_file)
{
  for(usint c = 0; c < CHARS; c++)
  {
    if(!(this->array[c])) { continue; }
    usint value = this->array[c]->select(0);
    psi_file.write((char*)&(value), sizeof(value));
    while(this->array[c]->hasNext())
    {
      value = this->array[c]->selectNext();
      psi_file.write((char*)&value, sizeof(value));
    }
  }
}

uchar*
RLCSA::readBWT()
{
  usint n = this->data_size + this->number_of_sequences;

  uchar* bwt = new uchar[n];
  memset(bwt, 0, n);

  for(usint c = 0; c < CHARS; c++)
  {
    RLEVector* vector = this->array[c];
    if(vector != 0)
    {
      bwt[vector->select(0)] = c;
      while(vector->hasNext()) { bwt[vector->selectNext()] = c; }
    }
  }

  return bwt;
}

//--------------------------------------------------------------------------

usint
RLCSA::psi(usint index)
{
  if(index >= this->data_size)
  {
    return this->data_size + this->number_of_sequences;
  }

  usint c = this->getCharacter(index);
  return this->array[c]->select(index - this->index_ranges[c].first);
}

pair_type
RLCSA::psi(usint index, usint max_length)
{
  if(index >= this->data_size)
  {
    return pair_type(this->data_size + this->number_of_sequences, 0);
  }

  usint c = this->getCharacter(index);
  return this->array[c]->selectRun(index - this->index_ranges[c].first, max_length);
}

//--------------------------------------------------------------------------

pair_type
RLCSA::getSequence(usint number)
{
  if(number >= this->number_of_sequences) { return EMPTY_PAIR; }

  pair_type result;
  if(number == 0)
  {
    result.first = 0;
    result.second = this->end_points->select(number);
  }
  else
  {
    if(this->sa_samples == 0) { return EMPTY_PAIR; }
    usint d = this->sa_samples->getSampleRate();
    result.first = d * ((this->end_points->select(number - 1) / d) + 1);
    result.second = this->end_points->selectNext();
  }

  return result;
}

usint
RLCSA::reportSize(bool print)
{
  usint bytes = 0, temp = 0;

  for(usint c = 0; c < CHARS; c++)
  {
    if(this->array[c]) { temp += this->array[c]->reportSize(); }
  }
  if(print) { std::cout << "Psi array:       " << temp << std::endl; }
  bytes += temp;

  if(this->support_locate || this->support_display)
  {
    temp = this->sa_samples->reportSize();
    if(print) { std::cout << "SA samples:      " << temp << std::endl; }
    bytes += temp;
  }

  temp = sizeof(*this) + this->end_points->reportSize();
  if(print) { std::cout << "RLCSA overhead:  " << temp << std::endl; }
  bytes += temp;

  if(print)
  {
    std::cout << "Total size:      " << bytes << std::endl;
    std::cout << std::endl;
  }

  return bytes;
}

//--------------------------------------------------------------------------

void
RLCSA::buildCharIndexes(usint* distribution)
{
  this->data_size = buildRanges(distribution, this->index_ranges);

  usint i = 0, c = 0;
  for(; c < CHARS; c++)
  {
    if(!isEmpty(this->index_ranges[c]))
    {
      this->text_chars[i] = c;
      i++;
    }
  }
  this->chars = i;

  this->index_rate = std::max((this->data_size + CHARS - 1) / CHARS, (usint)1);
  usint current = 0;

  for(c = 0, i = 0; c < this->chars; c++)
  {
    pair_type range = this->index_ranges[this->text_chars[c]];
    while(current <= range.second)
    {
      this->index_pointers[i] = c;
      current += this->index_rate;
      i++;
    }
  }
}

usint
buildRanges(usint* distribution, pair_type* index_ranges)
{
  if(distribution == 0 || index_ranges == 0) { return 0; }

  usint sum = 0;
  for(usint c = 0; c < CHARS; c++)
  {
    if(distribution[c] == 0)
    {
      if(sum == 0) { index_ranges[c] = EMPTY_PAIR; }
      else         { index_ranges[c] = pair_type(sum, sum - 1); }
    }
    else
    {
      index_ranges[c].first = sum;
      sum += distribution[c];
      index_ranges[c].second = sum - 1;
    }
  }

  return sum;
}


} // namespace CSA
