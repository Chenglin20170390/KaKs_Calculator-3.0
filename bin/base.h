/******************************************************************
* Copyright (c) CNCB-NGDC, BIG, CAS
* All rights reserved.
 
* Filename: base.h
* Abstract: Declaration of base class for KaKs methods.

* Version: 1.0
* Author: Zhang Zhang  (zhanghzhang@genomics.org.cn)
* Date: Feb.2, 2005

* Version: 1.2 
* Author: Zhang Zhang  (zhanghzhang@genomics.org.cn)
* Date: Apr. 2006

* Version: 3.0
* Author: Zhang Zhang  (zhanghzhang@big.ac.cn)
* Date: Jan.2, 2021
******************************************************************/
#pragma warning(disable:4786)		//Disable warning when using vector
#pragma warning(disable:4018)		//Disable warning when using vector
#pragma warning(disable:4244)		//Disable warning when using vector


#if !defined(BASE_H)
#define  BASE_H

#define VERSION		"Version 3.0, Nov. 2021"
#define AUTHOR_NAME	"Zhang Zhang"
#define AUTHOR_MAIL	"zhangzhang@big.ac.cn"

#define PACKAGE_NAME "KaKs_Calculator"
#define PACKAGE_REF "KaKs_Calculator 3.0: calculating selective pressure on coding and non-coding sequences. Genomics Proteomics Bioinformatics 2021"
#define KAKS_NAME "KaKs"
#define KAKS_DESC "calculating selective pressure on coding sequences."
#define KNKS_NAME "KnKs"
#define KNKS_DESC "calculating selective pressure on non-coding sequences."

#define CODONLENGTH 3			//Length of codon
#define DNASIZE 4			//A C G T
#define XSIZE DNASIZE*DNASIZE   //Size of one group AXX (X=A,C,G,T) 
#define CODON 64				//Codon Size
#define NULL 0					//Zero
#define NA -1					//Not Available
#define NCODE	33				//Number of genetic codes
#define NNCODE  NCODE*2			//Double of the number genetic codes		
#define SMALLVALUE 1e-6			//Value near to zero
#define NUMBER_OF_RATES	6		//Number of substitution rates
#define MODELCOUNT	14			//Number of candidate models

#define gammap(x,alpha) (alpha*(1-pow(x,-1/alpha)))
#define square(a) ((a)*(a))

#define min2(a,b) ((a)<(b)?(a):(b))
#define max2(a,b) ((a)>(b)?(a):(b))
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

/* Stanard lib of C++ */
#include<string>
#include<iostream>
#include<sstream>
#include<fstream>
#include<vector>
#include<map>
#include<stdlib.h>
#include<math.h>
#include<time.h>
#include<string.h>


using namespace std;

/*****Global variables*****/
extern const char* transl_table[];	//Genetic codon table
extern int genetic_code;			//ID of codon table from 1 to 33
extern int GeneticCode[][64];

extern string seq_name;				//Pairwise sequences' name
extern unsigned long length;		//Length of compared sequences
extern double GC[4];				//GC Contents of entire sequences(GC[0]) and three codon positions (GC[1--3])
//End of Global variables


/* Convert one type to any other type */
template<class out_type,class in_value>
	out_type CONVERT(const in_value & t) {
		stringstream stream;
		//Put the value 't' into the stream
		stream<<t;			
		out_type result;
		//Put the stream into the 'result'
		stream>>result;

		return result;
	}

class Codon {
	
public:
	Codon() {}
	
	Codon(string RorD, string q, string aa="") {
		RobDiv = RorD;
		Quarter = q;
	}
	/* Robustness or Diversity */
	string RobDiv;
	/* Quarter Location: 1:low gc, 2:gc p1, 3:gc p2, 4:high gc */
	string Quarter;
	/* Amino acid */
	string aa;
};

typedef map<int, string> IntString;
typedef map<string, int> StringInt;
typedef map<string, Codon> StringCodon;


class Base {

public:
	Base();

	/* Read seqeunces and their names */
	bool readAXTSeq(string filename, vector<string> &vec_names, vector<string> &vec_seqs);
	/* Write the content into file */
	bool writeFile(string output_filename, const char* result);

	/* Check pairwise noncoding sequences valid or not */
	bool checkPairwiseNoncoding(string seq, string &msg); 
	/* Check pairwise coding sequences valid or not */
	bool checkPairwiseCoding(string &seq, string &msg);

	/* Parse results */
	string parseOutput();	
	/* Format string for outputing into file */
	void addString(string &result, string str, string flag="\t");

	/* Generate a radnom integer */
	int getRandom();
	
	/* Convert a char-T,C,A,G into a digit 0,1,2,3, respectively */
	int  convertChar(char ch);
	/* Convert a digit-0,1,2,3 into a char T,C,A,G, respectively */
	char convertInt(int ch);
	/* Convert a string to uppercase */
	string stringtoUpper(string str);
	
	/* Calculate the amino acid of codon by codon */
	char getAminoAcid(string codon);
	/* Calculate the amino acid of codon by codon's id*/
	char getAminoAcid(int id);
	/* Get the number of stop codon in a given genetic code table */
	int getNumNonsense(int genetic_code);

	/* Return the codon's id from codon table */
	int getID(string codon);
	/* Return a codon according to the id */
	string getCodon(int IDcodon);
	/* Get GCC of entire sequences and of three codon positions */
	void getGCContent(string str, int cds=1);

	/* Sum array's elements */
	double sumArray(double x[], int end, int begin=0);
	int sumArray(int x[], int end, int begin=0);

	/* Init value to array */
	int initArray(double x[], int n, double value=0.0);
	int initArray(int x[], int n, int value=0);

	/* Elements in array are mutipled by scale */
	int scaleArray(double scale, double x[], int n);
	/* Sqrt of the sum of the elements' square */
	double norm(double x[], int n);
	/* Copy array's values one by one: to[] = from[] */
	int copyArray(double from[], double to[], int n);
	/* Sum of 'n' products multiplied by two elements x[], y[] */
	double innerp(double x[], double y[], int n);
	/* Set x[i,j]=0 when x!=j and x[i,j]=1 when x=j */
	int initIdentityMatrix(double x[], int n);

	/* Compute p-value by Fisher exact test to justify the validity of ka/ks */
	double fisher(double cs, double us, double cn, double un);
	/* factorial */
	double factorial(double n);


public:
	/* Name of method for calculating ka/ks */
	string name;
	/* Sysnonymous sites(S) and nonsysnonymous sites(N) */
	double S, N;
	/* Number of sysnonymous differences(Sd) and nonsysnonymous differences(Nd), snp=Sd+Nd */
	double Sd, Nd, snp;
	/* Number of sysnonymous substitutions(ks) and nonsysnonymous substitutions(ka) per site */
	double Ka, Ks;
	/* Standard Errors for Ka and Ks */
	double SEKa, SEKs;
	/* Transitional(Si) and transversional(Vi) substitutions */
	double Si[5], Vi[5];
	/* 0-fold, 2-fold, 4-fold */
	double L[5];
	/* Total Numbers of substitutions per i-th type site: K[i]=A[i]+B[i] */
	double K[5];
	/* 	   T  C  A  G
		T  -  6  7  8  
		C  0  -  9  10
		A  1  3  -  11
		G  2  4  5  -				*/	

	/* Transition/transversion rate ratio, tc: between pyrimidines, ag: between purines */
	double kappa, kappatc, kappaag;

	/* Divergence distance, substitutions per site: t = (S*Ks+N*Ka)/(S+N)*/
	double t;
	/* Maximum Likelihood Value */
	double lnL;
	/* Akaike Information Criterion (AICc) */
	double AICc;
	/* Akaike Weight in Model Selection */
	double AkaikeWeight;
	/* Name of mutation model according to AICc */
	string model;
	/* Transition/transversion mutation ratio(s) or substitution rates */
	double KAPPA[NUMBER_OF_RATES];

	/* Store Maximum Likilhood results for Model Selection and Model Averaging */
	struct MLResult {
		string result;	//estimated results
		double AICc;	//the value of a modified AIC
		double freq[CODON];		//Codon frequency
		double rate[NUMBER_OF_RATES];	//Six substitution rates
		double w;	//Ka/Ks
		double t;	//divergence distance
	};
        
	//ID(1~64) to Codon
	IntString ID2Codon;
	//Codon to ID
	StringInt Codon2ID;
        //64 Codons
	StringCodon Codon64;	
};

#endif


