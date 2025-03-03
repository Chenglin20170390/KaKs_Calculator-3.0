/************************************************************
* Copyright (c) CNCB-NGDC, BIG, CAS
* All rights reserved.
 
* Filename: KaKs.h
* Abstract: Declaration of KAKS class including several methods.

* Version: 1.0
* Author: Zhang Zhang  (zhanghzhang@genomics.org.cn)
* Date: Feb.2, 2005

* Version: 3.0
* Author: Zhang Zhang  (zhanghzhang@big.ac.cn)
* Date: Jan.2, 2021
*************************************************************/
#if !defined(KAKS_H)
#define  KAKS_H

#include "base.h"
#include "NG86.h"
#include "LWL85.h"
#include "LPB93.h"
#include "GY94.h"
#include "YN00.h"
#include "MYN.h"
#include "MSMA.h"

using namespace std;

/* KAKS class */
class KAKS: public Base {
	
public:	
	KAKS();
	~KAKS();

	/* Main function to call kaks's methods */
	bool Run(int argc, const char* argv[]);		
	/* Read and Calculate seq, called in "Run" main function */
	bool ReadCalculateSeq(string filename);
	
	/* Get the result for Windows, depending on the bool isOK4Win */
	string getResult4Win();

	/* Initialize class, ready for running */
	int Initialize();
	/* Unitialize class, for unloading */
	int Uninitialize();

protected:		
	/* Use several methods to calculate ka/ks */
	bool calculateKaKs();
	/* Show help information */
	void helpInfo();
	/* Show help information */
	void programInfo();

	/* NONE: an in-house algorithm in BIG, that is NG86 without correction */
	void start_NONE();
	/* NG86 */
	void start_NG86();
	/* LWL85 */
	void start_LWL85();
	/* Modified LWL85 */
	void start_MLWL85();
	/* LPB93 */
	void start_LPB93();
	/* Modified LPB93 */
	void start_MLPB93();
	/* GY94 */
	void start_GY94();	
	/* YN00 */
	void start_YN00();
	/* MYN */
	void start_MYN();	
	/* Model Selection and Model Averaging */
	void start_MSMA();
	


	/* Check the sequence whether is valid or not */
	bool checkValid(string name, string str);
	/* Parse the input parameters */
	bool parseParameter(int argc, const char* argv[]);
	/* Show input parameters' information on screen */
	void showParaInfo();
	/* Get title information for writing into file */
	string getTitleInfo();

public:
	/* Methods' name and reference */
	vector<string> method_name;
	vector<string> method_ref;
	
	/* Parameters' title in outputing file */
	vector<string> titleInfo;

	/* Results for windows parser that shows results on ListCtrl */
	string result4Win;
	/* Flag for outputing in Windows */
	bool isOK4Win;
	
	/* File name for output */
	string output_filename;
	/* Sequence file name */
	string seq_filename;

	/* Flag for whether to run NG86, MLWL85, MLPB93, GY94, YN00, MYN, MS/A=model selection/averaging */
	bool none, ng86, lwl85, lpb93, yn00, mlwl85, mlpb93, gy94, myn06, ms06, ma06;	
	/* Number of compared pairwise sequences */
	unsigned long number;	//Maybe too many
	/* Running time:  */
	int hh, mm, ss;

protected:
	/* File name for detailed results for model selection */
	string detail_filename;
	/* Detailed results */
	string details; 
	
private:
	/* The temporary results for write into file */
	string result;
	/* Output stream */
	ofstream os;
	/* A pair of sequence */
	string seq1, seq2;
}; 

#endif



