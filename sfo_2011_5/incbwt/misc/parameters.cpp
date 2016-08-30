#include <cstdlib>
#include <iostream>
#include <map>

#include "parameters.h"


namespace CSA
{


Parameters::Parameters()
{
}

Parameters::~Parameters()
{
}

bool
Parameters::contains(const std::string& key)
{
  return (this->parameters.find(key) != this->parameters.end());
}

usint
Parameters::get(const std::string& key)
{
  std::map<std::string, usint>::iterator iter = this->parameters.find(key);
  if(iter == this->parameters.end()) { return 0; }
  return iter->second;
}

usint
Parameters::get(const parameter_type& param)
{
  return this->get(param.first);
}

void
Parameters::set(const std::string& key, usint value)
{
  this->parameters[key] = value;
}

void
Parameters::set(const parameter_type& param)
{
  this->set(param.first, param.second);
}

void
Parameters::read(std::ifstream& file)
{
  while(file)
  {
    std::string key;
    std::string c;
    usint value;

    file >> key >> c >> value;
    if(c == "=") { this->parameters[key] = value; }
  }
}

void Parameters::read(const std::string& file_name)
{
  std::ifstream file(file_name.c_str(), std::ios_base::binary);
  if(!file)
  {
    std::cerr << "Cannot open parameter file " << file_name << " for reading!" << std::endl;
    return;
  }
  this->read(file);
  file.close();
}

void
Parameters::print()
{
  this->write(std::cout);
  std::cout << std::endl;
}

void
Parameters::write(std::ostream& stream)
{
  for(std::map<std::string, usint>::iterator iter = this->parameters.begin(); iter != this->parameters.end(); iter++)
  {
    stream << iter->first << " = " << iter->second << std::endl;
  }
}

void
Parameters::write(const std::string& file_name)
{
  std::ofstream file(file_name.c_str(), std::ios_base::binary);
  if(!file)
  {
    std::cerr << "Cannot open parameter file " << file_name << " for writing!" << std::endl;
    return;
  }
  this->write(file);
  file.close();
}


} // namespace CSA
