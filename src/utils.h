/**
 * Copyright (c) 2016-present, Facebook, Inc.
 * All rights reserved.
 *
 * This source code is licensed under the BSD-style license found in the
 * LICENSE file in the root directory of this source tree. An additional grant
 * of patent rights can be found in the PATENTS file in the same directory.
 */

#ifndef FASTTEXT_UTILS_H
#define FASTTEXT_UTILS_H

#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include <ios>
#include <vector>

namespace fasttext
{

namespace utils
{

int64_t size(std::ifstream &);
void seek(std::ifstream &, int64_t);

std::vector<std::string> splitString(std::string str, char delimiter);
std::vector<std::string> splitStringDelims(std::string str, std::string delims);

} // namespace utils

} // namespace fasttext

#endif
