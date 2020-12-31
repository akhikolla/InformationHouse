/*
 * MTB_Metaphone.cpp
 *
 *  Created on: 06.07.2017
 *      Author: schnell-42
 */

#include "MTB_Metaphone.h"

/**
* Five values in the English language
*/
string vowels = "AEIOU";
/**
* Variable used in Metaphone algorithm
*/
string frontv = "EIY";
/**
* Variable used in Metaphone algorithm
*/
string varson = "CSPTG";
/**
* The max code length for metaphone is 4
*/
int maxCodeLen = 4;

string MTB_Metaphone::getName() {
  return "Metaphone";
}

double MTB_Metaphone::getRelativeValue(string var1, string var2) {
  if (metaphone(var1).compare(metaphone(var2)) == 0) {
    return 1.0;
  }
  return 0.0;
}

double MTB_Metaphone::getAbsoluteValue(string var1, string var2) {

  return getRelativeValue(var1,  var2);
}

/**
* A class to generate phonetic code. The initial Java implementation, William B. Brogden. December, 1997 Permission given by wbrogden for code to be used
* anywhere. "Hanging on the Metaphone" by Lawrence Philips <i>Computer Language </i> of Dec. 1990, p 39
*
* @author wbrogden@bga.com
* @author bayard@generationjava.com
* @author Tim O'Brien
* @version $Id: Metaphone.java,v 1.1 2005/11/23 17:37:12 Administrator Exp $
*/

void setMaxCodeLen(int maxCodeLenIN) {
  maxCodeLen = maxCodeLenIN;
}


int getMaxCodeLen() {
  return maxCodeLen;
}

/**
* Find the metaphone value of a String. This is similar to the soundex algorithm, but better at finding similar sounding words. All input is converted to
* upper case. Limitations: Input format is expected to be a single ASCII word with only characters in the A - Z range, no punctuation or numbers.
*
* @param txt
*            String to find the metaphone code for
* @return A metaphone code corresponding to the String supplied
*/
string metaphone(string txt) {
// #pragma clang diagnostic push
// #pragma clang diagnostic ignored "-Wtautological-compare"
  int mtsz = 0;
  bool hard = false;
  if ((txt.length() == 0)) {
    return "";
  }
  // single character is itself
  if (txt.length() == 1) {
    toUpper(txt);
    return txt;
  }
  char *inwd = new char[txt.size() + 1];
  toUpper(txt);
  strcpy(inwd, txt.c_str());
  string tmpS;
  string local; // manipulate
  string code; // output
  // handle initial 2 characters exceptions
  switch (inwd[0]) {
  case 'K':
  case 'G':
  case 'P': /* looking for KN, etc */
if (inwd[1] == 'N') {
  local.append(txt.substr(1, txt.length() - 1));
} else {
  local.append(inwd);
}
break;
  case 'A': /* looking for AE */
if (inwd[1] == 'E') {
  local.append(txt.substr(1, txt.length() - 1));
} else {
  local.append(inwd);
}
break;
  case 'W': /* looking for WR or WH */
if (inwd[1] == 'R') { // WR -> R
  local.append(txt.substr(1, txt.length() - 1));
  break;
}
if (inwd[1] == 'H') {
  local.append(txt.substr(1, txt.length() - 1));
  local[0] = 'W'; // WH -> W
} else {
  local.append(inwd);
}
break;
  case 'X': /* initial X becomes S */
inwd[0] = 'S';
    local.append(inwd);
    break;
  default:
    local.append(inwd);
  } // now local has working string with initials fixed
  unsigned wdsz = local.length();
  unsigned n = 0;
  while ((mtsz < maxCodeLen) // max code size of 4 works well
           && (n < wdsz)) {
    char symb = local[n];
    // remove duplicate letters except C
    if ((symb != 'C') && (n > 0) && (local[n - 1] == symb)) {
      n++;
    } else { // not dup
      switch (symb) {
      case 'A':
      case 'E':
      case 'I':
      case 'O':
      case 'U':
        if (n == 0) {
          code.push_back(symb);
          mtsz++;
        }
        break; // only use vowel if leading char
      case 'B':
        if ((n > 0) && !((n + 1) == wdsz) // not MB at end of word
              && (local[n - 1] == 'M')) {
          code.push_back(symb);
        } else {
          code.push_back(symb);
        }
        mtsz++;
        break;
      case 'C': // lots of C special cases
        /* discard if SCI, SCE or SCY */
        if ((n > 0) && (local[n - 1] == 'S') && ((n + 1) < wdsz)
              && (frontv.find(local[n + 1]) != string::npos)) {
          break;
        }
        tmpS = local;
        if (tmpS.find("CIA", n) == n) { // "CIA" -> X
          code.push_back('X');
          mtsz++;
          break;
        }
        if (((n + 1) < wdsz)
              && (frontv.find(local[n + 1])!= string::npos)) {
          code.push_back('S');
          mtsz++;
          break; // CI,CE,CY -> S
        }
        if ((n > 0) && (tmpS.find("SCH", n - 1) == (n - 1))) { // SCH->sk
          code.push_back('K');
          mtsz++;
          break;
        }
        if (tmpS.find("CH", n) == n) { // detect CH
          if ((n == 0) && (wdsz >= 3) // CH consonant -> K consonant
                && (vowels.find(local[2]) == string::npos)) {
            code.push_back('K');
          } else {
            code.push_back('X'); // CHvowel -> X
          }
          mtsz++;
        } else {
          code.push_back('K');
          mtsz++;
        }
        break;
      case 'D':
        if (((n + 2) < wdsz) // DGE DGI DGY -> J
              && (local[n + 1] == 'G')
              && (frontv.find(local[n + 2]) != string::npos)) {
              code.push_back('J');
          n += 2;
        } else {
          code.push_back('T');
        }
        mtsz++;
        break;
      case 'G': // GH silent at end or before consonant
        if (((n + 2) == wdsz) && (local[n + 1] == 'H')) {
          break;
        }
        if (((n + 2) < wdsz) && (local[n + 1] == 'H')
              && (vowels.find(local[n + 2]) == string::npos)) {
          break;
        }
        tmpS = local;
        if (((n > 0) && (tmpS.find("GN", n) == n))
              || (tmpS.find("GNED", n) == n)) {
          break; // silent G
        }
        if ((n > 0) && (local[n - 1] == 'G')) {
          hard = true;
        } else {
          hard = false;
        }
        if (((n + 1) < wdsz)
              && (frontv.find(local[n + 1]) != string::npos)
              && (!hard)) {
              code.push_back('J');
        } else {
          code.push_back('K');
        }
        mtsz++;
        break;
      case 'H':
        if ((n + 1) == wdsz) {
          break; // terminal H
        }
        if ((n > 0)
              && (varson.find(local[n - 1]) != string::npos)) {
          break;
        }
        /*if (vowels.find(local[n + 1]) >= 0) { */ // unsigned int is always greater_equal zero
          code.push_back('H');
          mtsz++; // Hvowel
        /*}*/
        break;
      case 'F':
      case 'J':
      case 'L':
      case 'M':
      case 'N':
      case 'R':
        code.push_back(symb);
        mtsz++;
        break;
      case 'K':
        if (n > 0) { // not initial
          if (local[n - 1] != 'C') {
            code.push_back(symb);
          }
        } else {
          code.push_back(symb); // initial K
        }
        mtsz++;
        break;
      case 'P':
        if (((n + 1) < wdsz) && (local[n + 1] == 'H')) {
          // PH -> F
          code.push_back('F');
        } else {
          code.push_back(symb);
        }
        mtsz++;
        break;
      case 'Q':
        code.push_back('K');
        mtsz++;
        break;
      case 'S':
        tmpS = local;
        if ((tmpS.find("SH", n) == n)
              || (tmpS.find("SIO", n) == n)
              || (tmpS.find("SIA", n) == n)) {
              code.push_back('X');
        } else {
          code.push_back('S');
        }
        mtsz++;
        break;
      case 'T':
        tmpS = local; // TIA TIO -> X
        if ((tmpS.find("TIA", n) == n)
              || (tmpS.find("TIO", n) == n)) {
          code.push_back('X');
          mtsz++;
          break;
        }
        if (tmpS.find("TCH", n) == n) {
          break;
        }
        // substitute numeral 0 for TH (resembles theta after all)
        if (tmpS.find("TH", n) == n) {
          code.push_back('0');
        } else {
          code.push_back('T');
        }
        mtsz++;
        break;
      case 'V':
        code.push_back('F');
        mtsz++;
        break;
      case 'W':
      case 'Y': // silent if not followed by vowel
        if (((n + 1) < wdsz)) {
          code.push_back(symb);
          mtsz++;
        }
        break;
      case 'X':
        code.push_back('K');
        code.push_back('S');
        mtsz += 2;
        break;
      case 'Z':
        code.push_back('S');
        mtsz++;
        break;
      } // end switch
      n++;
    } // end else from symb != 'C'
    if (mtsz > maxCodeLen) {
      code = code.substr(0,maxCodeLen);
    }
  }
  delete[] inwd;
// #pragma clang diagnostic pop
  return code;
}

