#ifndef UTILS_H
#define UTILS_H

#include <fstream>
#include <list>

#include "definitions.h"


namespace CSA
{


std::streamoff fileSize(std::ifstream& file);
std::streamoff fileSize(std::ofstream& file);


std::ostream& operator<<(std::ostream& stream, pair_type data);


void readRows(std::ifstream& file, std::list<std::string>& rows, bool skipEmptyRows);


} // namespace CSA


#endif // UTILS_H
