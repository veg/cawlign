// Copyright (c) 2018 Daniel Zilles
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

#ifndef CONFIGPARSER_HPP
#define CONFIGPARSER_HPP

#include <string>
#include <vector>
#include <map>
#include <iostream>
#include <sstream>

typedef std::map<std::string, std::vector<std::string>> configList; 

class ConfigParser {

  public:
  ConfigParser(std::ifstream& configFile);

    template<typename T>
    T aConfig(std::string section, std::string name, size_t pos = 0);
    template<typename T>
    std::vector<T> aConfigVec(std::string section, std::string name);
    
    
  private:
    static void       handleMissingKey (std::string);
    configList mConfigurations;

};

template <>
bool ConfigParser::aConfig<bool>(std::string section, std::string name, size_t pos);

template <typename T>
T ConfigParser::aConfig(std::string section, std::string configName, size_t pos) {

  T tmp;
    
  const auto& mConfigRef = mConfigurations;
  auto search = mConfigRef.find(section + " - " + configName);
    
  if (search == mConfigRef.end()) {
      handleMissingKey (std::string("Could not find required configuration section ") + section + std::string(" key ") + configName);
  }
    
  std::istringstream iss(search->second[0]);

  if (search->second[0].find( "0x" ) != std::string::npos)
    iss >> std::hex >> tmp;
  else
    iss >> std::dec >> tmp;

  return tmp;
}

template <>
std::vector<bool> ConfigParser::aConfigVec<bool>(std::string section, std::string name);

template <typename T>
std::vector<T> ConfigParser::aConfigVec(std::string section, std::string configName) {

    
  const auto& mConfigRef = mConfigurations;
  auto search = mConfigRef.find(section + " - " + configName);
      
  if (search == mConfigRef.end()) {
    handleMissingKey (std::string("Could not find required configuration section ") + section + std::string(" key ") + configName);
  }


  std::vector<T>  tmp(search->second.size());
  for (unsigned i = 0; i < search->second.size(); i++) {

    std::istringstream iss(search->second[i]);

    if (search->second[i].find( "0x" ) != std::string::npos)
      iss >> std::hex >> tmp[i];
    else
      iss >> std::dec >> tmp[i];
  }
  return tmp;
}
#endif
