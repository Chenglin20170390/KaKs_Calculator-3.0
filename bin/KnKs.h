/************************************************************
* Copyright (c) CNCB-NGDC, BIG, CAS
* All rights reserved.

* Filename: KnKs.h
* Abstract: Declaration of KNKS class for noncoding sequencess.

* Version: 1.0
* Author: Zhang Zhang  (zhanghzhang@big.ac.cn)
* Date: Jan.2, 2021
*************************************************************/
#if !defined(KNKS_H)
#define  KNKS_H

#include "base.h"
#include "GY94.h"

using namespace std;

/* KAKS class */
class KNKS : public Base {
public:
	KNKS();
	~KNKS();

	/* Main function to call kaks's methods */
	bool Run(int argc, const char* argv[]);
	/* Read and Calculate seq, called in "Run" main function */
	bool ReadCalculateSeq(string filename);

	/* Initialize class, ready for running */
	int Initialize();
	/* Unitialize class, for unloading */
	int Uninitialize();
	/* Check whether str is numeric */
	bool isNum(string str);

	/* Get results for Windows */
	string getResult4Win();

protected:
	/* Use several methods to calculate ka/ks */
	bool calculateKnKs(string seq1, string seq2);
	/* Get corrected distance based on JC69 */
	double getDistanceJC69(double d);
	/* Get corrected distance based on JC69 */
	double getDistanceK2P(double ts, double tv);
	/* Get corrected distance based on JC69 */
	double getDistanceHKY(double ts, double tv, double pa, double pt, double pg, double pc);
	
	/* Show help information */
	void helpInfo();
	/* Program information */
	void programInfo();

	/* Check the sequence whether is valid or not */
	bool checkNCValid(string name, string str);
	bool checkValid(string name, string str);
	/* Parse the input parameters */
	bool parseParameter(int argc, const char* argv[]);
	/* Show input parameters' information on screen */
	void showParaInfo();
	/* Get title information for writing into file */
	string getNCTitleInfo();
	string getCDSTitleInfo();

public:
	double Kn;
	double nc_len, nc_ts, nc_tv, nc_kappa, nc_GC;
	/* Methods' name and reference */
	vector<string> method_name;
	vector<string> method_ref;

	/* Parameters' title in outputing file */
	vector<string> titleInfo;

	/* Results for windows parser that shows results on ListCtrl */
	string result4Win;
	/* Flag of results for Windoes */
	bool isOK4Win;

	/* Noncoding sequence file name */
	string input_nc_filename;
	/* Coding sequence file name */
	string input_coding_filename;
	/* Input mutation rate */
	double mutation_rate;
	/* File name for noncoding output */
	string output_nc_filename;
	/* File name for coding output */
	string output_coding_filename;
	/* Sequence count */
	int number;
	/* Running time */
	int hh, mm, ss;

private:
	/* Noncoding results */
	string nc_result;
	/* Detailed coding results */
	string coding_result;
	/* Output stream */
	ofstream os;
	/* A pair of sequence */
	string seq1, seq2;

	/* CDS method */
	GY94 zz;

};

#endif

