//============================================================================
// Name        : Hashing.cpp
// Author      : Rukasz
// Version     :
// Copyright   : Your copyright notice
// Description : Ansi-style
//============================================================================

#include "Hashing.h"

template <typename HashMethod, typename InputContainer, typename KeyContainer>
std::string createHmac(const InputContainer& input, const KeyContainer& key)
{
  std::string hash = hmac<HashMethod>(&input[0], input.size(), &key[0], key.size());
  return hash;
}

//use HMAC on strings
string useHMAC(string inputData,string secretKey) {

  string mdString = createHmac<SHA256>(inputData, secretKey);

  return mdString;
} //end useHMAC
