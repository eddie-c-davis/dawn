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

#include "dawn/Support/FileUtil.h"

namespace dawn {

StringRef getFilename(StringRef path) { return path.substr(path.find_last_of('/') + 1); }

StringRef getFilenameWithoutExtension(StringRef path) {
  auto filename = getFilename(path);
  return filename.substr(0, filename.find_last_of("."));
}

StringRef getExtension(StringRef filename) {
  return filename.substr(filename.find_last_of(".") - 1);
}

std::string readFile(const std::string& file) {
  std::ifstream is(file);
  return std::string((std::istreambuf_iterator<char>(is)), std::istreambuf_iterator<char>());
}

std::string readPipe(const std::string& path) {
  std::string command = path + " 2>&1";
  std::array<char, 1024> buffer;
  std::string result;

  FILE* pipe = popen(command.c_str(), "r");
  if(pipe) {
    while(fgets(buffer.data(), 1024, pipe)) {
      result += buffer.data();
    }
    pclose(pipe);
  }

  return result;
}

} // namespace dawn
