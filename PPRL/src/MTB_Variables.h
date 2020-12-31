/*
 * MTB_Variables.h
 *
 *  Created on: 23.05.2017
 *      Author: schnell-42
 */

#ifndef MTB_VARIABLES_H_
#define MTB_VARIABLES_H_

using namespace std;

/**
 * Holds Variables for the data set.
 */

class MTB_Variables{
public:
	/**
	 * constructor
	 */
	MTB_Variables();
	~MTB_Variables();

	/**
	 * Return the list of variables
	 *
	 * @return
	 */
	vector<string> getVariableList();
	/**
	 * Sets the list of variables
	 * @param
	 */
	void setVariableList(vector<string>);
	/**
	 * Return the number of variables
	 *
	 * @return
	 */
	int getNumberOfVariables();
	/**
	* Sets the number of variables
	* @param
	*/
	void setNumberOfVariables(int);

private:
	vector<string> variableList;
	int number;

};



#endif /* MTB_VARIABLES_H_ */
