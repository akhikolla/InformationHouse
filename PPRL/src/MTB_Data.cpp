/*
 * MTB_Data.cpp
 *
 *  Created on: 24.05.2017
 *      Author: schnell-42
 */

#include "MTB_Data.h"
MTB_Data::MTB_Data() {
}

MTB_Data::~MTB_Data() {
}


/**
 * MTB_StringData
 */
MTB_StringData::MTB_StringData() {
}

MTB_StringData::MTB_StringData(string variablename, string data, string id) {
	this->data = data;
	this->id = id;
	this->variablename = variablename;
}

MTB_StringData::~MTB_StringData() {
}


string MTB_StringData::getData() {
	return this->data;
}


void MTB_StringData::setData(string data) {
	this->data = data;
}


string MTB_StringData::getID() {
	return this->id;
}


void MTB_StringData::setID(string id) {
	this->id = id;
}


string MTB_StringData::getVariablename() {
	return this->variablename;
}


void MTB_StringData::setVariablename(string variablename) {
	this->variablename = variablename;
}


/**
 *  MTB_StringVectorData
 */
MTB_StringVectorData::MTB_StringVectorData() {
}

MTB_StringVectorData::MTB_StringVectorData(string variablename,
		vector<string> id, vector<string> data,  vector<string> blockingData) {
	this->data = data;
	this->id = id;
	this->blockingData = blockingData;
	this->variablename = variablename;
}

MTB_StringVectorData::~MTB_StringVectorData() {
}


vector<string> MTB_StringVectorData::getData() {
	return this->data;
}


void MTB_StringVectorData::setData(vector<string> data) {
	this->data = data;
}


vector<string> MTB_StringVectorData::getBlockingData() {
	return this->blockingData;
}


void MTB_StringVectorData::setBlockingData(vector<string> blockingData) {
	this->blockingData = blockingData;
}


vector<string> MTB_StringVectorData::getID() {
	return this->id;
}


void MTB_StringVectorData::setID(vector<string> id) {
	this->id = id;
}


string MTB_StringVectorData::getVariablename() {
	return this->variablename;
}


void MTB_StringVectorData::setVariablename(string variablename) {
	this->variablename = variablename;
}

