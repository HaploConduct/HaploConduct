#include "utils.h"


namespace CSA
{


std::streamoff
fileSize(std::ifstream& file)
{
  std::streamoff curr = file.tellg();

  file.seekg(0, std::ios::end);
  std::streamoff size = file.tellg();
  file.seekg(0, std::ios::beg);
  size -= file.tellg();

  file.seekg(curr, std::ios::beg);
  return size;
}

std::streamoff
fileSize(std::ofstream& file)
{
  std::streamoff curr = file.tellp();

  file.seekp(0, std::ios::end);
  std::streamoff size = file.tellp();
  file.seekp(0, std::ios::beg);
  size -= file.tellp();

  file.seekp(curr, std::ios::beg);
  return size;
}

std::ostream&
operator<<(std::ostream& stream, pair_type data)
{
  return stream << "(" << data.first << ", " << data.second << ")";
}

void
readRows(std::ifstream& file, std::list<std::string>& rows, bool skipEmptyRows)
{
  while(file)
  {
    std::string buf;
    std::getline(file, buf);
    if(skipEmptyRows && buf.length() == 0) { continue; }
    rows.push_back(buf);
  }
}


} // namespace CSA
