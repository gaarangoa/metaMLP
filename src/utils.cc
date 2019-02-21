/**
 * Copyright (c) 2016-present, Facebook, Inc.
 * All rights reserved.
 *
 * This source code is licensed under the BSD-style license found in the
 * LICENSE file in the root directory of this source tree. An additional grant
 * of patent rights can be found in the PATENTS file in the same directory.
 */

#include "utils.h"

namespace fasttext
{

namespace utils
{

int64_t size(std::ifstream &ifs)
{
    ifs.seekg(std::streamoff(0), std::ios::end);
    return ifs.tellg();
}

void seek(std::ifstream &ifs, int64_t pos)
{
    ifs.clear();
    ifs.seekg(std::streampos(pos));
}

std::vector<std::string> splitStringDelims(std::string str, std::string delimiter)
{
    std::vector<std::string> internal;

    auto start = 0U;
    auto end = str.find(delimiter);
    std::string token;

    while (end != std::string::npos)
    {
        token = str.substr(start, end - start);
        internal.push_back(token);
        start = end + delimiter.length();
        end = str.find(delimiter, start);
    }

    token = str.substr(start, end - start);
    internal.push_back(token);
    return internal;
}

std::vector<std::string> splitString(std::string str, char delimiter)
{
    std::vector<std::string> internal;
    std::stringstream ss(str); // Turn the string into a stream.
    std::string tok;

    while (getline(ss, tok, delimiter))
    {
        internal.push_back(tok);
    }

    return internal;
}

} // namespace utils

} // namespace fasttext
