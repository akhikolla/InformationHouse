/*
 * MTB_Data.h
 *
 *  Created on: 24.05.2017
 *      Author: schnell-42
 */

#ifndef MTB_DATA_H_
#define MTB_DATA_H_

#include <string>
#include <vector>

using namespace std;

/**
 * Holds data (one field) plus its variable name and id
 */

class MTB_Data {
public:
	/**
	 * constructor
	 */
	MTB_Data();
	~MTB_Data();

private:
	string varibleName = "NA";
	string id = "0";

};

/**
 * Holds one string plus its variable name and id
 */
class MTB_StringData {
public:
	/**
	 *
	 * @param variablename
	 * @param data
	 * @param id
	 */
	MTB_StringData(string variablename, string data, string id);
	/**
	 * empty constructor
	 */
	MTB_StringData();
	~MTB_StringData();

	/**
	 * Return data
	 *
	 * @param
	 * @param
	 */
	string getData();
	/**
	 * Sets the data
	 *
	 * @param
	 * @param
	 */
	void setData(string data);

	/**
	 * Return ID
	 *
	 * @param
	 * @param
	 */
	string getID();
	/**
	 * Sets the ID
	 *
	 * @param
	 * @param
	 */
	void setID(string id);
	/**
	 * Return data
	 *
	 * @param
	 * @param
	 */
	string getVariablename();
	/**
	 * Sets the data at variablename and id
	 *
	 * @param
	 * @param
	 */
	void setVariablename(string variablename);

private:
	string id;
	string data;
	string variablename;
	//int freq =1;
};

/**
 * Holds one string plus its variable name and id
 */
class MTB_StringVectorData {
public:
	/**
	 *
	 * @param variablename
	 * @param data
	 * @param id
	 */
	MTB_StringVectorData(string variablename, vector<string> id,
			vector<string> data, vector<string> blockingData);
	/**
	 * empty constructor
	 */
	MTB_StringVectorData();
	~MTB_StringVectorData();

	/**
	 * Returns data
	 *
	 * @param
	 * @param
	 */
	vector<string> getData();
	/**
	 * Sets data
	 *
	 * @param
	 * @param
	 */
	void setData(vector<string> data);

	/**
	 * Returns IDs
	 *
	 * @param
	 * @param
	 */
	vector<string> getID();
	/**
	 * Sets IDs
	 *
	 * @param
	 * @param
	 */
	void setID(vector<string> id);
	/**
	 * Returns data
	 *
	 * @param
	 * @param
	 */
	vector<string> getBlockingData();
	/**
	 * Sets blocking data
	 *
	 * @param
	 * @param
	 */
	void setBlockingData(vector<string> data);
	/**
	 * Returns variable name
	 *
	 * @param
	 * @param
	 */
	string getVariablename();
	/**
	 * Sets the data at variablename and id
	 *
	 * @param
	 * @param
	 */
	void setVariablename(string variablename);

private:
	vector<string> id;
	vector<string> data;
	vector<string> blockingData;
	string variablename;
};

#endif /* MTB_DATA_H_ */
