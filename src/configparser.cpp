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

#include "configparser.hpp"
#include <fstream>
#include <sstream>
#include <vector>
#include <iostream>
#include <algorithm>
#include <stdexcept>

#include "argparse.hpp"
#include "alignment.h"

using namespace std;

/**
 * Constructor for `ConfigParser`.
 *
 * This function reads and parses a configuration file. It processes lines to strip whitespace, 
 * ignores comments, and splits key-value pairs using the `=` symbol. If the key is within a section 
 * (indicated by `[...]`), the key is prefixed with the section name. Multiple values for a key 
 * are stored in a vector. Parsed configurations are stored in a map with section and key names 
 * combined as the key.
 *
 * @param configFile A reference to an `ifstream` representing the open configuration file.
 * @throw std::runtime_error if a parsing error occurs in the configuration file.
 */
ConfigParser::ConfigParser(ifstream& configFile) {

  if (!configFile.is_open()) return;
  string line;
  string sectionName;
  string configName;
  string rawConfigContent;
  vector<string> configContent;

  size_t lineNbr = 1;
 
  while( getline(configFile, line) ) {

    line.erase(std::remove_if( line.begin(), line.end(),
    [l = std::locale {}](auto ch) { return std::isspace(ch, l); }), line.end());

    if (line.empty()) {
      lineNbr++;
      continue;
    }

    if (line.at(0) == ';') {
      lineNbr++;
      continue;
    }

    if(line.at(0) == '[' && line.at(line.length()-1) == ']') {
      sectionName = line.substr(1, line.length()-2);
      continue;
    }

    size_t equalSignPos = line.find('=');
    if(equalSignPos != string::npos) {
      configName = line.substr(0, equalSignPos);
      rawConfigContent = line.substr(equalSignPos+1, line.length()-1);
     
      stringstream ss(rawConfigContent);
      while (ss.good()) {
        string tmp;
        getline(ss, tmp, ',');
        configContent.push_back(tmp);
      }
    }

    else {
      string errorMessage = ":" + to_string(lineNbr) + ": parsing error \n" + line;
      throw std::runtime_error(errorMessage);
    }

    mConfigurations.insert(std::make_pair(sectionName + " - " + configName, configContent));
    configContent.clear();
    lineNbr++;
  }
  configFile.close();
}

/**
 * Template specialization for retrieving boolean configuration values.
 *
 * This function retrieves a specific boolean configuration value from a section and configuration 
 * name. The value is parsed as either `true` (for "true", "TRUE", or "1") or `false` (for "false", 
 * "FALSE", or "0"). If the value does not match any of these, the function defaults to `false`.
 *
 * @param section The section name in the configuration file.
 * @param configName The key within the section.
 * @param pos The position of the value in the vector (if there are multiple values).
 * @return `true` if the configuration value is recognized as true, otherwise `false`.
 */
template <>
bool ConfigParser::aConfig<bool>(std::string section, std::string configName, size_t pos) {

  bool tmp;
  std::string config = mConfigurations[section + " - " + configName][pos];
  std::istringstream iss(config);

  if (config == "true" ||
      config == "TRUE" ||
      config == "1") {
    return true;
  }

  else if (config == "false" ||
           config == "FALSE" ||
      config == "0") {
    return false;
  }
  else
    return false;
}

/**
 * Template specialization for retrieving a vector of boolean configuration values.
 *
 * This function retrieves a vector of boolean values associated with a configuration key in a section. 
 * Each value is parsed as either `true` (for "true", "TRUE", or "1") or `false` (for "false", 
 * "FALSE", or "0"). If a value does not match any of these, it defaults to `false`.
 *
 * @param section The section name in the configuration file.
 * @param configName The key within the section.
 * @return A vector of boolean values parsed from the configuration.
 */
template <>
std::vector<bool> ConfigParser::aConfigVec<bool>(std::string section, std::string configName) {

  std::vector<std::string> config;

  config = mConfigurations[section + " - " + configName];

  std::vector<bool>  tmp(config.size());
  for (unsigned i = 0; i < config.size(); i++) {

    if (config[i] == "true" ||
        config[i] == "TRUE" ||
        config[i] == "1") {
      tmp[i] = true;
    }

    else if (config[i] == "false" ||
             config[i] == "FALSE" ||
             config[i] == "0") {
      tmp[i] = false;
    }
    else
      tmp[i] = false;

  }
  return tmp;
}

/**
 * Handles missing configuration keys by triggering an error.
 *
 */
void ConfigParser::handleMissingKey (std::string error) {
    argparse::ERROR_NO_USAGE(error.c_str());
}













