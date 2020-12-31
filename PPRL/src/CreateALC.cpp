#include "CreateALC.h"

using namespace std;


string createALC(vector<string> vars, vector<bool> soundex, string password) {
	string alc;
	if (vars.size()!=soundex.size()){
		Rcpp::Rcerr << "The length of the input vector of the variables must have the same size as the soundex vector!" << endl;
		alc = "";
		return alc;
	}

	for (int i = 0; i < (int) vars.size(); i++) {
		createALCHelper(vars[i], soundex[i]);
		alc = alc + vars[i];
	}
	return useHMAC(alc, password);
}

void createALCHelper(string &var, bool soundex) {
  if (soundex) {
    soundexC(var);
  } else {
    replaceNonAscii(var);
    //var = deletePunct(var); TODO überprüfen
  }
}

void soundexC(string &var){
  standardisation(var);
  //delete vowels
  char first = var[0];
  deleteVowels(var);
  deleteYWH(var);
  string rest;
  codeConsonants(var);
  replaceDuplicates(rest);
  cutToThree(rest);
  fillZero(rest);
  var = first + rest;

  //TODO Doppelte zu Beginn streichen?
}

void codeConsonants(string& var) {

	if (var.size()>= 1) {
		var.replace(0, 1, "");
	}
	ReplaceAllSubstr(var, "B", "1");
	ReplaceAllSubstr(var, "F", "1");
	ReplaceAllSubstr(var, "P", "1");
	ReplaceAllSubstr(var, "V", "1");

	ReplaceAllSubstr(var, "C", "2");
	ReplaceAllSubstr(var, "G", "2");
	ReplaceAllSubstr(var, "J", "2");
	ReplaceAllSubstr(var, "K", "2");
	ReplaceAllSubstr(var, "Q", "2");
	ReplaceAllSubstr(var, "S", "2");
	ReplaceAllSubstr(var, "X", "2");
	ReplaceAllSubstr(var, "Z", "2");

	ReplaceAllSubstr(var, "D", "3");
	ReplaceAllSubstr(var, "T", "3");

	ReplaceAllSubstr(var, "L", "4");

	ReplaceAllSubstr(var, "M", "5");
	ReplaceAllSubstr(var, "N", "5");

	ReplaceAllSubstr(var, "R", "6");


}

void replaceDuplicates(string& var) {

  while (var.find("11") < var.size()) {
    ReplaceAllSubstr(var, "11", "1");
  }
  while (var.find("22") < var.size()) {
    ReplaceAllSubstr(var, "22", "2");
  }
  while (var.find("33") < var.size()) {
    ReplaceAllSubstr(var, "33", "3");
  }
  while (var.find("44") < var.size()) {
    ReplaceAllSubstr(var, "44", "4");
  }
  while (var.find("55") < var.size()) {
    ReplaceAllSubstr(var, "55", "5");
  }
  while (var.find("66") < var.size()) {
    ReplaceAllSubstr(var, "66", "6");
  }
}

void deleteVowels(string &var) {
	//delete vowels, if not at the beginning
	char first = var[0];
	string keepFirst = "";
	if (first=='A'|| first=='E'||first=='I'|| first=='O'||first=='U'){
		keepFirst = first;
	}
	ReplaceAllSubstr(var, "A", "");
	ReplaceAllSubstr(var, "E", "");
	ReplaceAllSubstr(var, "I", "");
	ReplaceAllSubstr(var, "O", "");
	ReplaceAllSubstr(var, "U", "");
	var =  keepFirst+var;
}


void deleteYWH(string &var) {
	//delete Y, W and H if not at the beginning
	char first = var[0];
	string keepFirst = "";
	if (first=='H'|| first=='Y'||first=='W'){
		keepFirst = first;
	}
	ReplaceAllSubstr(var, "Y", "");
	ReplaceAllSubstr(var, "H", "");
	ReplaceAllSubstr(var, "W", "");
	var =  keepFirst+var;
}



void cutToThree(string &var) {
	if (var.size() > 3)
		var = var.substr(0,3);
}

void fillZero(string& var) {
	if (var.size() == 0)
		var = "000";
	if (var.size() == 1)
		var = var + "00";
	if (var.size() == 2)
		var = var + "0";

}

