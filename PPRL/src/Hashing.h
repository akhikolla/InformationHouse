//============================================================================
// Name        : Hashing.h
// Author      : Rukasz
// Version     :
// Copyright   : Your copyright notice
// Description : Ansi-style
//============================================================================

#ifndef HASHING_H
#define HASHING_H

#include <iostream>
#include <stdio.h>
#include <cstring>
#include "sha256.h"
#include "hmac.h"

using namespace std;

/**
 * Hashes data with the secret key.
 *
 * @param key
 */
string useHMAC(string inputData, string secretKey);

#endif
