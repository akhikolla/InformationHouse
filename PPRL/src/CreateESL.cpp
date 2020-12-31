#include "CreateESL.h"


using namespace std;

string createESL(vector<string> vars, vector<vector<int>> code, string password) {
	string esl;
	for (int i = 0; i < (int) vars.size(); i++) {
		vars[i] = createESLHelper(vars[i], code[i]);
		esl  = esl  + vars[i];
	}
	return useHMAC(esl, password);
}

string createESLHelper(string var, vector<int> code) {
	string esl;
	replaceNonAscii(var);
	//var = deletePunct(var);
	//fill up date with zeros if it consists of only one digit
	string checkVarDigit = "0123456789";
	if (var.size() == 1 && strstr(checkVarDigit.c_str(),var.c_str()))
		{
		   var = "0" + var;
	}
	if(!(code.size() ==1 && code[0] == 0)&& code.size() >0 ){
		for (unsigned i = 0 ; i < code.size(); i++){
			if (code[i]>0 && (unsigned)code[i] <=var.size())
			esl = esl + var.at(code[i]-1);
		}
	}
	else
		esl = var;


	return esl;
}

