//============================================================================
// Name        : Bloomencoder.cpp
// Author      : Rukasz
// Version     :
// Copyright   : Your copyright notice
// Description : Ansi-style
//============================================================================

#include "Standardisation.h"

void ReplaceAllSubstr( string& source, const string& from, const string& to )
{
  string newString;
  newString.reserve( source.length() );  // avoids a few memory allocations

  string::size_type lastPos = 0;
  string::size_type findPos;

  while( string::npos != ( findPos = source.find( from, lastPos )))
  {
    newString.append( source, lastPos, findPos - lastPos );
    newString += to;
    lastPos = findPos + from.length();
  }

  // Care for the rest after last occurrence
  newString += source.substr( lastPos );

  source.swap( newString );
}


void ReplaceAllChar( string& source, char from, const string& to )
{
  string newString;
  newString.reserve( source.length() );  // avoids a few memory allocations

  string::size_type lastPos = 0;
  string::size_type findPos;

  while( string::npos != ( findPos = source.find( from, lastPos )))
  {
    newString.append( source, lastPos, findPos - lastPos );
    newString += to;
    lastPos = findPos+1;
  }

  // Care for the rest after last occurrence
  newString += source.substr( lastPos );

  source.swap( newString );
}


void toUpper(string &s){
  transform(s.begin(), s.end(),s.begin(), ::toupper);
}

void standardisation(string &var) {

  delFuell(var);
  delTitel(var);
  replaceNonAscii(var);
} //end standardisation

void preprocess(string& var) {
  toUpper(var);


} //end preprocess

/**
* checks  char
*/
bool invalidChar (char c)
{
  return !((c>=32 &&c <=126));//||(c==32)||(c==39)||(c==58)||(c==59)||(c==95));
}

bool invalidCharPunct (char c)
{
  return !((c <48)||(c>=57 && c <=65)||(c>90));
}

void replaceUmlaut(char c1, char c2){

}

void replaceNonAscii(string& var) {
  preprocess(var);
  ReplaceAllSubstr(var, "\"", "");

  ReplaceAllChar(var, 192, "A");
  ReplaceAllChar(var, 193, "A");
  ReplaceAllChar(var, 194, "A");
  ReplaceAllChar(var, 195, "A");
  ReplaceAllChar(var, 196, "AE");
  ReplaceAllChar(var, 197, "A");
  ReplaceAllChar(var, 198, "AE");
  ReplaceAllChar(var, 224, "A");
  ReplaceAllChar(var, 225, "A");
  ReplaceAllChar(var, 226, "A");
  ReplaceAllChar(var, 227, "A");
  ReplaceAllChar(var, 228, "AE");
  ReplaceAllChar(var, 229, "A");
  ReplaceAllChar(var, 230, "AE");

  ReplaceAllChar(var, 199, "C");
  ReplaceAllChar(var, 231, "C");

  ReplaceAllChar(var, 200, "E");
  ReplaceAllChar(var, 201, "E");
  ReplaceAllChar(var, 202, "E");
  ReplaceAllChar(var, 203, "E");
  ReplaceAllChar(var, 232, "E");
  ReplaceAllChar(var, 233, "E");
  ReplaceAllChar(var, 234, "E");
  ReplaceAllChar(var, 235, "E");

  ReplaceAllChar(var, 204, "I");
  ReplaceAllChar(var, 205, "I");
  ReplaceAllChar(var, 206, "I");
  ReplaceAllChar(var, 207, "I");
  ReplaceAllChar(var, 236, "I");
  ReplaceAllChar(var, 237, "I");
  ReplaceAllChar(var, 238, "I");
  ReplaceAllChar(var, 239, "I");

  ReplaceAllChar(var, 210, "O");
  ReplaceAllChar(var, 211, "O");
  ReplaceAllChar(var, 212, "O");
  ReplaceAllChar(var, 213, "O");
  ReplaceAllChar(var, 214, "OE");
  ReplaceAllChar(var, 242, "O");
  ReplaceAllChar(var, 243, "O");
  ReplaceAllChar(var, 244, "O");
  ReplaceAllChar(var, 245, "O");
  ReplaceAllChar(var, 246, "OE");

  ReplaceAllChar(var, 223, "SS");
  ReplaceAllSubstr(var, "Ş", "S");
  ReplaceAllSubstr(var, "ş", "S");

  ReplaceAllChar(var, 252, "UE");
  ReplaceAllChar(var, 220, "UE");
  ReplaceAllChar(var, 217, "U");
  ReplaceAllChar(var, 218, "U");
  ReplaceAllChar(var, 219, "U");
  ReplaceAllChar(var, 249, "U");
  ReplaceAllChar(var, 250, "U");
  ReplaceAllChar(var, 251, "U");
  ReplaceAllChar(var, 252, "U");
  //remove everything left, that is non ascii
  var.erase(remove_if(var.begin(),var.end(), invalidChar), var.end());
} //end replaceNonAscii

void delTitel(string& var) {
  preprocess(var);
  ReplaceAllSubstr(var, "BARONESSE ", " ");
  ReplaceAllSubstr(var, "BARONIN ", " ");
  ReplaceAllSubstr(var, "BARON ", " ");
  ReplaceAllSubstr(var, "BR ", " ");
  ReplaceAllSubstr(var, "DOKTOR ", " ");
  ReplaceAllSubstr(var, "DR ", " ");
  ReplaceAllSubstr(var, "DRIN ", " ");
  ReplaceAllSubstr(var, "FREIFRAU ", " ");
  ReplaceAllSubstr(var, "FREIHERRIN ", " ");
  ReplaceAllSubstr(var, "FREIHERR ", " ");
  ReplaceAllSubstr(var, "FREIIN ", " ");
  ReplaceAllSubstr(var, "FUERSTIN ", " ");
  ReplaceAllSubstr(var, "FUERST ", " ");
  ReplaceAllSubstr(var, "GRAEFIN ", " ");
  ReplaceAllSubstr(var, "GROSSHERZOEGIN ", " ");
  ReplaceAllSubstr(var, "GROSSHERZOG ", " ");
  ReplaceAllSubstr(var, "GROSSHERZOGIN ", " ");
  ReplaceAllSubstr(var, "HERZOEGIN ", " ");
  ReplaceAllSubstr(var, "HERZOG ", " ");
  ReplaceAllSubstr(var, "HERZOGIN ", " ");
  ReplaceAllSubstr(var, "KRONRINZESSIN ", " ");
  ReplaceAllSubstr(var, "KRONPRINZ ", " ");
  ReplaceAllSubstr(var, "KURPRINZESSIN ", " ");
  ReplaceAllSubstr(var, "KURPRINZ ", " ");
  ReplaceAllSubstr(var, "LANDGRAEFIN ", " ");
  ReplaceAllSubstr(var, "LANDGRAF ", " ");
  ReplaceAllSubstr(var, "MARKGRAEFIN ", " ");
  ReplaceAllSubstr(var, "MARKGRAF ", " ");
  ReplaceAllSubstr(var, "MED ", " ");
  ReplaceAllSubstr(var, "PD ", " ");
  ReplaceAllSubstr(var, "PRINZESSIN ", " ");
  ReplaceAllSubstr(var, "PRINZ ", " ");
  ReplaceAllSubstr(var, "PROFESSORIN ", " ");
  ReplaceAllSubstr(var, "PROFESSOR ", " ");
  ReplaceAllSubstr(var, "PROFIN ", " ");
  ReplaceAllSubstr(var, "PROF ", " ");
  ReplaceAllSubstr(var, "REICHGRAEFIN ", " ");
  ReplaceAllSubstr(var, "REICHGRAF ", " ");
  ReplaceAllSubstr(var, "REICHSFREIHERRIN ", " ");
  ReplaceAllSubstr(var, "REICHSFREIHERR ", " ");
  ReplaceAllSubstr(var, "REICHSVIKAR ", " ");
  ReplaceAllSubstr(var, "RITTER ", " ");
  ReplaceAllSubstr(var, "SR ", " ");
  ReplaceAllSubstr(var, "GRAF ", " ");
} //end delTitel

void delFuell(string & var) {
  preprocess(var);

  ReplaceAllSubstr(var, " AL ", " ");
  ReplaceAllSubstr(var, " AM ", " ");
  ReplaceAllSubstr(var, " AN ", " ");
  ReplaceAllSubstr(var, " AUF ", " ");
  ReplaceAllSubstr(var, " BEN ", " ");
  ReplaceAllSubstr(var, " DA ", " ");
  ReplaceAllSubstr(var, " DEL ", " ");
  ReplaceAllSubstr(var, " DEM ", " ");
  ReplaceAllSubstr(var, " DEN ", " ");
  ReplaceAllSubstr(var, " DER ", " ");
  ReplaceAllSubstr(var, " DES ", " ");
  ReplaceAllSubstr(var, " DI ", " ");
  ReplaceAllSubstr(var, " DOS ", " ");
  ReplaceAllSubstr(var, " DU ", " ");
  ReplaceAllSubstr(var, " EL ", " ");
  ReplaceAllSubstr(var, " EN ", " ");
  ReplaceAllSubstr(var, " ET ", " ");
  ReplaceAllSubstr(var, " LA ", " ");
  ReplaceAllSubstr(var, " LE ", " ");
  ReplaceAllSubstr(var, " L ", " ");
  ReplaceAllSubstr(var, " MC ", " ");
  ReplaceAllSubstr(var, " MAC ", " ");
  ReplaceAllSubstr(var, " MED ", " ");
  ReplaceAllSubstr(var, " O ", " ");
  ReplaceAllSubstr(var, " TER ", " ");
  ReplaceAllSubstr(var, " UND ", " ");
  ReplaceAllSubstr(var, " VAN ", " ");
  ReplaceAllSubstr(var, " VON ", " ");
  ReplaceAllSubstr(var, " VOR ", " ");
  ReplaceAllSubstr(var, " VOM ", " ");
  ReplaceAllSubstr(var, " V ", " ");
  ReplaceAllSubstr(var, " Y ", " ");
  ReplaceAllSubstr(var, " ZUM ", " ");
  ReplaceAllSubstr(var, " ZUR ", " ");
  ReplaceAllSubstr(var, " ZU ", " ");
  ReplaceAllSubstr(var, " DE ", " ");
  ReplaceAllSubstr(var, " D ", " ");
  ReplaceAllSubstr(var, " A ", " ");

} //end delFuell


