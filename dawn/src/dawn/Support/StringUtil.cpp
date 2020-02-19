//===--------------------------------------------------------------------------------*- C++ -*-===//
//                          _
//                         | |
//                       __| | __ ___      ___ ___
//                      / _` |/ _` \ \ /\ / / '_  |
//                     | (_| | (_| |\ V  V /| | | |
//                      \__,_|\__,_| \_/\_/ |_| |_| - Compiler Toolchain
//
//
//  This file is distributed under the MIT License (MIT).
//  See LICENSE.txt for details.
//
//===------------------------------------------------------------------------------------------===//

#include "dawn/Support/StringUtil.h"
#include <sstream>

namespace dawn {

extern std::string decimalToOrdinal(int dec) {
  std::string decimal = std::to_string(dec);
  std::string suffix;

  switch(decimal.back()) {
  case '1':
    suffix = "st";
    break;
  case '2':
    suffix = "nd";
    break;
  case '3':
    suffix = "rd";
    break;
  default:
    suffix = "th";
  }

  if(dec > 10)
    suffix = "th";

  return decimal + suffix;
}

std::string indent(const std::string& string, int amount) {
  // This could probably be done faster (it's not really speed-critical though)
  std::istringstream iss(string);
  std::ostringstream oss;
  std::string spacer(amount, ' ');
  bool firstLine = true;
  for(std::string line; std::getline(iss, line);) {
    if(!firstLine)
      oss << spacer;
    oss << line;
    if(!iss.eof())
      oss << "\n";
    firstLine = false;
  }
  return oss.str();
}

void tokenize(const std::string& str, char delim, std::vector<std::string>& elems) {
  std::string elem;
  std::stringstream ss(str);
  while(std::getline(ss, elem, delim)) {
    elems.emplace_back(std::move(elem));
  }
}

void tokenize(const std::string& str, char delim, std::vector<double>& elems) {
  std::string elem;
  std::stringstream ss(str);
  while(std::getline(ss, elem, delim)) {
    elems.emplace_back(std::stod(elem));
  }
}

} // namespace dawn
