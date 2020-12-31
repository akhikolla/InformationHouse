/*
 * MTB_Reth.cpp
 *
 *  Created on: 05.07.2017
 *      Author: schnell-42
 */
/**
 *  Reth, Hans-Peter von and Schek, Hans-Jörg (1977) Eine Zugriffsmethode für die phonetische Ähnlichkeitssuche,  Heidelberg Scientific Center technical reports
by IBM Deutschland GmbH
 * @return
 */

#include "MTB_Reth.h"

string MTB_Reth::getName() {
	return "Reth";
}
double MTB_Reth::getRelativeValue(string o1, string o2) {
	if (encodeReth(o1).compare(encodeReth(o2))==0)
			{
				return 1.0;
			}
return 0.0;
}

double MTB_Reth::getAbsoluteValue(string o1, string o2) {

return getRelativeValue( o1,  o2);
}

string encodeReth(string input)
	{
		replaceNonAscii(input);
		toUpper(input);

		ReplaceAllSubstr(input,"AA", "A");
		ReplaceAllSubstr(input,"AH", "A");
		ReplaceAllSubstr(input,"P", "B");
		ReplaceAllSubstr(input,"BB", "B");
		ReplaceAllSubstr(input,"PP", "B");
		ReplaceAllSubstr(input,"BP", "B");
		ReplaceAllSubstr(input,"PB", "B");
		ReplaceAllSubstr(input,"T", "D");
		ReplaceAllSubstr(input,"DD", "D");
		ReplaceAllSubstr(input,"TT", "D");
		ReplaceAllSubstr(input,"TH", "D");
		ReplaceAllSubstr(input,"DT", "D");
		ReplaceAllSubstr(input,"TD", "D");
		ReplaceAllSubstr(input,"EE", "E");
		ReplaceAllSubstr(input,"EH", "E");
		ReplaceAllSubstr(input,"AE", "E");
		ReplaceAllSubstr(input,"AEH", "E");
		ReplaceAllSubstr(input,"V", "F");
		ReplaceAllSubstr(input,"W", "F");
		ReplaceAllSubstr(input,"FF", "F");
		ReplaceAllSubstr(input,"PH", "F");
		ReplaceAllSubstr(input,"C", "G");
		ReplaceAllSubstr(input,"K", "G");
		ReplaceAllSubstr(input,"GG", "G");
		ReplaceAllSubstr(input,"KK", "G");
		ReplaceAllSubstr(input,"CK", "G");
		ReplaceAllSubstr(input,"CK", "G");
		ReplaceAllSubstr(input,"GK", "G");
		ReplaceAllSubstr(input,"KG", "G");
		ReplaceAllSubstr(input,"Y", "I");
		ReplaceAllSubstr(input,"IE", "I");
		ReplaceAllSubstr(input,"IE", "I");
		ReplaceAllSubstr(input,"IEH", "I");
		ReplaceAllSubstr(input,"LL", "L");
		ReplaceAllSubstr(input,"MM", "M");
		ReplaceAllSubstr(input,"NN", "N");
		ReplaceAllSubstr(input,"OO", "O");
		ReplaceAllSubstr(input,"OH", "O");
		ReplaceAllSubstr(input,"RR", "R");
		ReplaceAllSubstr(input,"SZ", "S");
		ReplaceAllSubstr(input,"SS", "S");
		ReplaceAllSubstr(input,"UH", "U");
		ReplaceAllSubstr(input,"GS", "X");
		ReplaceAllSubstr(input,"KS", "X");
		ReplaceAllSubstr(input,"CHS", "X");
		ReplaceAllSubstr(input,"CKS", "X");
		ReplaceAllSubstr(input,"C", "Z");
		ReplaceAllSubstr(input,"TZ", "Z");
		ReplaceAllSubstr(input,"OEH", "OE");
		ReplaceAllSubstr(input,"UEH", "UE");
		ReplaceAllSubstr(input,"AI", "AI");
		ReplaceAllSubstr(input,"AY", "AI");
		ReplaceAllSubstr(input,"EI", "AI");
		ReplaceAllSubstr(input,"EY", "AI");
		ReplaceAllSubstr(input,"EU", "OI");
		ReplaceAllSubstr(input,"AEU", "OI");
		ReplaceAllSubstr(input,"KW", "QU");
		ReplaceAllSubstr(input,"CH", "SCH");
		ReplaceAllSubstr(input,"ZIO", "TIO");
		ReplaceAllSubstr(input,"TUI", "ZUI");

		input = replaceEnding(input, "ER",'R');
		input = replaceEnding(input, "EL",'R');
		input = replaceEnding(input, "H",' ');
		ReplaceAllSubstr(input," ", "");
		return input;
	}


