/************************************************************
 * Copyright (C) CNCB-NGDC, BIG, CAS
 * All rights reserved.

 * Filename: KNKS.cpp
 * Abstract: Declaration of KNKS class including several methods.

 * Version: 1.0
 * Author: Zhang Zhang (zhanghzhang@genomics.org.cn)
 * Date: Feb.2, 2005

 * Version: 3.0
 * Author: Zhang Zhang (zhanghzhang@big.ac.cn)
 * Date: Jan.17, 2021
 *************************************************************/

#include "KnKs.h"

KNKS::KNKS() {

	string items_nc[] = { "Sequence", "Kn", "Ks", "Kn/Ks", "Length", "Substitutions", "Kappa", "GC" };
	string items_cds[] = { "CDS", "Ka", "Ka/Ks", "CDS-Length", "CDS-Substitutions", "CDS-Kappa", "GC(1:2:3)" };
	int i;
	for (i = 0; i < sizeof(items_nc) / sizeof(string); i++) titleInfo.push_back(items_nc[i]);
	for (i = 0; i < sizeof(items_cds) / sizeof(string); i++) titleInfo.push_back(items_cds[i]);

	//Load Methods' Names and References for -h in linux, also for windows' tool tip
	//method_name.push_back("ZZ");
	//method_ref.push_back("Zhang Z (2021) in preparation");

	Initialize();

}

KNKS::~KNKS() {

	Uninitialize();

	//Free memory
	titleInfo.clear();
	method_name.clear();
	method_ref.clear();
}

int KNKS::Initialize() {

	GY94 tmp("HKY");
	
	zz = tmp;

	result4Win = nc_result = coding_result = seq_name = seq1 = seq2 = "";
	input_nc_filename = input_coding_filename = output_nc_filename = output_coding_filename = "";
	nc_result = coding_result = "";
	genetic_code = 1;
	number = 0;

	mutation_rate = NA;

	nc_GC = 0.0;
	hh = mm = ss = 0;

	return 1;
}

int KNKS::Uninitialize() {

	if (os.is_open()) {
		os.close();
	}
	hh = mm = ss = 0;
	return 1;
}

string KNKS::getResult4Win() {
	if (isOK4Win) {
		return result4Win;
	}
	return "";
}

bool KNKS::checkValid(string name, string str) {
	bool flag = true;
	long i;

	try {

		//Check whether sequences are equal in length
		if (str.length() % 2 != 0) {
			cout << "[Error. The sequences are not equal in length.]";
			throw 1;
		}

		string str1 = str.substr(0, str.length() / 2);
		string str2 = str.substr(str.length() / 2, str.length() / 2);

		//Check whether (sequence length)/3==0
		if (str1.length() % 3 != 0 || str2.length() % 3 != 0) {
			cout << "[Error. The sequences are not codon-based alignment.]";
			throw 1;
		}

		//Delete gap and stop codon
		bool found;
		int j;
		for (i = 0; i < str1.length(); i += 3) {
			for (found = false, j = 0; j < 3 && !found; j++) {
				if (str1[j + i] == '-' || str2[j + i] == '-') {
					found = true;
				}
				str1[i + j] = toupper(str1[i + j]);
				str2[i + j] = toupper(str2[i + j]);
				if (convertChar(str1[i + j]) == -1 || convertChar(str2[i + j]) == -1) {
					found = true;
				}
			}

			if ((getAminoAcid(str1.substr(i, 3)) == '!') || (getAminoAcid(str2.substr(i, 3)) == '!')) {
				found = true;
			}

			if (found) {
				str1 = str1.replace(i, 3, "");
				str2 = str2.replace(i, 3, "");
				i -= 3;
			}
		}

		//pass value into private variables
		seq1 = str1;
		seq2 = str2;
		//pass value into extern variables
		seq_name = name;
		length = str1.length();
	}
	catch (...) {
		flag = false;
	}

	return flag;
}

/****************************************************
 * Function: ReadCalculateSeq
 * Input Parameter: string
 * Output: Read sequences, check sequences' validity
				  and calculate Ka and Ks.
 * Return Value: True if succeed, otherwise false.

 * Note: Using axt file for low memory
 *****************************************************/
bool KNKS::ReadCalculateSeq(string filename) {

	//Record the start time
	time_t time_start = time(NULL);

	bool flag = true;

	try {
		showParaInfo();	//Show information on display

		vector<string> vec_nc_names, vec_nc_seqs;
		//Read noncoding sequences
		if (!readAXTSeq(input_nc_filename, vec_nc_names, vec_nc_seqs)) {
			throw 1;
		}

		vector<string> vec_coding_names, vec_coding_seqs;
		//mutation_rate is provided by user
		if (mutation_rate!=NA) {
			Ks = mutation_rate;
		}
		else {//mutation rate needs to be inferred from adjacent coding sequences
			if (!readAXTSeq(input_coding_filename, vec_coding_names, vec_coding_seqs) ||
				vec_nc_names.size() != vec_coding_names.size()) {
				throw 1;
			}
		}

		//Output stream
		if (output_nc_filename!= "" && output_nc_filename.length() > 0) {
			os.open(output_nc_filename.c_str());
		}

		nc_result = getNCTitleInfo();
		coding_result = getCDSTitleInfo();

		string temp = "";
		for (int i = 0; i < vec_nc_seqs.size(); i++) {
			//Check str's validility and calculate
			cout << "[" << i + 1 << "] " << vec_nc_names[i] << "\t";
			if (checkNCValid(vec_nc_names[i], vec_nc_seqs[i])) {
				try {

					//Noncoding sequences: Kn 
					calculateKnKs(seq1, seq2);
					//cout << Kn;
					//cout << "[Noncoding Kn]" << "\t";
					
					if (mutation_rate!=NA) {//Ks is provided by user
						/*"Sequence", "Kn", "Ks", "Kn/Ks", "Length", "Substitutions", "Kappa", "GC" */
						addString(nc_result, vec_nc_names[i]);
						addString(nc_result, CONVERT<string>(Kn));
						addString(nc_result, CONVERT<string>(Ks));
						if (Ks < SMALLVALUE || Ks == NA || Kn == NA) {
							addString(nc_result, "NA");
						}
						else {
							addString(nc_result, CONVERT<string>(Kn / Ks));
						}
						addString(nc_result, CONVERT<string>(nc_len));
						addString(nc_result, CONVERT<string>(nc_ts + nc_tv));
						addString(nc_result, CONVERT<string>(nc_kappa));
						addString(nc_result, CONVERT<string>(nc_GC));

						//CDS items
						addString(nc_result, "NA");
						addString(nc_result, "NA");
						addString(nc_result, "NA");
						addString(nc_result, "NA");
						addString(nc_result, "NA");
						addString(nc_result, "NA");
						addString(nc_result, "NA", "\n");
					}
					else {//Ks is inferred from input CDS 
						if (checkValid(vec_coding_names[i], vec_coding_seqs[i])) {
							getGCContent(seq1 + seq2); 
							coding_result += zz.Run(seq1.c_str(), seq2.c_str());
							Ks = zz.Ks;
						}
						/*"Sequence", "Kn", "Ks", "Kn/Ks", "Length", "Substitutions", "Kappa", "GC" */
						addString(nc_result, vec_nc_names[i]);
						addString(nc_result, CONVERT<string>(Kn));
						addString(nc_result, CONVERT<string>(Ks));
						if (Ks < SMALLVALUE || Ks == NA || Kn == NA) {
							addString(nc_result, "NA");
						}
						else {
							addString(nc_result, CONVERT<string>(Kn / Ks));
						}
						addString(nc_result, CONVERT<string>(nc_len));
						addString(nc_result, CONVERT<string>(nc_ts + nc_tv));
						addString(nc_result, CONVERT<string>(nc_kappa));
						addString(nc_result, CONVERT<string>(nc_GC));

						/* "CDS", "Ka", "Ka/Ks", "CDS-Length", "CDS-Substitutions", "CDS-Kappa", "GC(1:2:3)" */
						addString(nc_result, vec_coding_names[i]);
						addString(nc_result, CONVERT<string>(zz.Ka));
						addString(nc_result, CONVERT<string>(zz.Ka / zz.Ks));
						addString(nc_result, CONVERT<string>(seq1.length()));
						addString(nc_result, CONVERT<string>(zz.snp));
						addString(nc_result, CONVERT<string>(zz.KAPPA[0]));
						//GC Content
						string tmp = "";
						tmp = CONVERT<string>(GC[0]);	tmp += "(";
						tmp += CONVERT<string>(GC[1]);	tmp += ":";
						tmp += CONVERT<string>(GC[2]);	tmp += ":";
						tmp += CONVERT<string>(GC[3]);	tmp += ")";
						addString(nc_result, tmp, "\n");
					}
					
					cout << "[OK]";
					number++;
					//cout << nc_result;

					//add a lock "isOK4Win" to avoid the program collapse 
					isOK4Win = false;
					result4Win += nc_result;
					isOK4Win = true;

					//Write into the file
					if (output_nc_filename.length() > 0 && os.is_open()) {
						os<<nc_result.c_str();
						os.flush();
					}
					nc_result = "";

				}catch(...) {
					cout << "[Error in calculating]";
				}
			}//end of if
			cout << endl;
		}//end of for

		cout << nc_result;

	}
	catch (...) {
		flag = false;
	}

	//Time used for running
	time_t t = time(NULL) - time_start;
	hh = t / 3600;
	mm = (t % 3600) / 60;
	ss = t - (t / 60) * 60;

	return flag;
}

/* Get corrected distance based on JC69 */
double KNKS::getDistanceJC69(double p) {
	double d = 1 - (4 * p) / 3;
	if (d < 0.0) {
		d = NA;
	}
	else {
		d = log(d);
		if (d > 0.0)
			d = NA;
		else
			d = (-3.0)*d / 4.0;
	}

	return d;
}

/* Get corrected distance based on K2P */
/* d = -1/2 log(1-2S-V) - 1/4 log(1-2V)  */
double KNKS::getDistanceK2P(double ts, double tv) {
	double d1 = 1 - 2 * ts - tv;
	double d2 = 1 - 2 * tv;
	double d = NA;
	if (d1 > 0.0 && d2 > 0.0) {
		d = -0.5 * log(d1) - 0.25 * log(d2);
	}
	return d;
}

/**************************************************
* Get corrected distance based on HKY 
* Equations 1.27 and 1.28 in Book "Computational Molecular Evolution" by Ziheng Yang 
***************************************************/
double KNKS::getDistanceHKY(double ts, double tv, double pa, double pt, double pg, double pc) {
	double pr = pa + pg;
	double py = pt + pc;
	double p1 = pt * pc / py + pa * pg / pr;
	double p2 = pt * pc * pr / py + pa * pg * py / pr;

	double a = 1 - ts / (2 * p1) - tv * p2 / (2 * (pt*pc*pr + pa*pg*py));
	double b = 1 - tv / (2 * py * pr);
	double d = NA;
	if (a > 0.0 && b > 0.0) {
		a = -log(a);
		b = -log(b);
		d = 2 * p1 * a - 2 * (p2 - py*pr) * b;
		//nc_kappa = a / b - 1.0;
		nc_kappa = 2 * ts / tv;
	}
	return d;
}

/**************************************************
 * Function: checkNCValid
 * Input Parameter: string, string
 * Output: Check validity of pairwise noncoding sequences
 * Return Value: True if succeed, otherwise false.
 ***************************************************/
bool KNKS::checkNCValid(string name, string str) {
	bool flag = true;

	try {

		//Check whether sequences are equal in length
		if (str.length() % 2 != 0) {
			cout << "[Error. The sequences are not equal in length.]";
			throw 1;
		}

		//pass value into private variables
		seq1 = str.substr(0, str.length() / 2);
		seq2 = str.substr(str.length() / 2, str.length() / 2);

		//pass value into extern variables
		seq_name = name;
		length = seq1.length();
	}
	catch (...) {
		flag = false;
	}
	return flag;
}

/**************************************************
 * Function: Run
 * Input Parameter: int, const char* []
 * Output: Calculate Ka/Ks and output.
 * Return Value: void
 ***************************************************/
bool KNKS::Run(int argc, const char* argv[]) {

	bool flag = true;

	try {
		//Judge whether input parameters are legal
		if (!parseParameter(argc, argv)) {
			throw 1;
		}

		//Read sequences and calculate Kn & Ks
		ReadCalculateSeq(input_nc_filename);

		//writeFile(output_nc_filename, nc_result.c_str());
		//Output results
		cout << "Outputing results: ";
		cout << output_nc_filename;
		//Output CDS details
		if ( writeFile(output_coding_filename, coding_result.c_str()) == true) {
			cout << "\t" << output_coding_filename;
		}
		cout << endl;

		//Print on display
		cout << "Mission accomplished. (Time elapsed: ";
		if (hh) cout << hh << ":" << mm << ":" << ss;
		else cout << mm << ":" << ss;
		cout<< ")" << endl;
	}
	catch (...) {
		flag = false;
	}

	return flag;
}

bool KNKS::isNum(string str) {
	stringstream sin(str);
	double d;
	char c;
	if (!(sin >> d))
		return false;
	if (sin >> c)
		return false;
	return true;
}

/**************************************************
 * Function: parseParameter
 * Input Parameter: int, const char* []
 * Output: Parse the input parameters
 * Return Value: bool
 ***************************************************/
bool KNKS::parseParameter(int argc, const char* argv[]) {

	bool flag = true;
	int i;
	string temp;

	try {

		//No parameter
		if (argc == 1) {
			throw 1;
		}
		else if (argc == 2) {//help information
			temp = argv[1];
			if (temp == "-H" || temp == "-h") {
				helpInfo();
				flag = false;
			}
			else {
				throw 1;
			}
		}
		else {
			//cout<<"test"<<endl;
			//parse parameters
			int input_nc_flag = 0, input_coding_flag = 0, output_nc_flag = 0, output_coding_flag = 0, codeflag = 0;
			for (i = 1; i < argc; i++) {

				temp = stringtoUpper(argv[i]);
				//Input coding axt file
				if (temp == "-I") {
					if ((i + 1) < argc && input_nc_flag == 0) {
						input_nc_filename = argv[++i];
						input_nc_flag++;
					}
					else {
						throw 1;
					}
				}
				//Input noncoding axt file or numeric mutation rate
				else if (temp == "-J") {
					if ((i + 1) < argc && input_coding_flag == 0) {
						string tmp = argv[++i];
						if (isNum(tmp)) mutation_rate = CONVERT<double>(tmp);
						else input_coding_filename = tmp;
						input_coding_flag++;
					}
					else {
						throw 1;
					}
				}
				//Output file
				else if (temp == "-O") {
					if ((i + 1) < argc && output_nc_flag == 0) {
						output_nc_filename = argv[++i];
						output_nc_flag++;
					}
					else {
						throw 1;
					}
				}//Genetic Code Table
				else if (temp == "-C") {
					if ((i + 1) < argc && codeflag == 0) {
						genetic_code = CONVERT<int>(argv[++i]);
						if (genetic_code < 1 || genetic_code > NCODE || strlen(transl_table[2 * (genetic_code - 1)]) < 1)
							throw 1;
						codeflag++;
					}
					else {
						throw 1;
					}
				}//Details for coding estimates
				else if (temp == "-D") {
					if ((i + 1) > argc) throw 1;
					output_coding_filename = argv[++i];
					output_coding_flag++;
				}
				else throw 1;
			}

			//If no input or output file, report error
			if (input_nc_flag == 0 || input_coding_flag == 0 || output_nc_flag == 0) {
				throw 1;
			}
			//If Ks is provided but CDS results are required, report error
			if (mutation_rate != NA && output_coding_flag == 1) {
				throw 1;
			}
		}
	}
	catch (...) {
		cout << "Input parameter(s) error." << endl;
		cout << "For help information: " << KNKS_NAME << " -h" << endl;
		flag = false;
	}

	return flag;
}

/*******************************************************
 * Function: Calculate
 * Input Parameter: void
 * Output: Calculate KNKS and output results.
 * Return Value: void
 *
 * Note:
 ********************************************************/
bool KNKS::calculateKnKs(string seq1, string seq2) {

	bool flag = true;

	try {
		//Estimate Ka and Ks
		Kn = nc_len = nc_ts = nc_tv = nc_kappa = 0.0;
		double pi[4];
		initArray(pi, 4);
		for (int i = 0; i < seq1.length(); i++) {
			//t,c,a,g for 0,1,2,3; otherwise -1 
			int c1 = convertChar(seq1[i]);
			int c2 = convertChar(seq2[i]);
			if (c1 > -1 && c2 > -1) {
				nc_len++;
				pi[c1]++;
				pi[c2]++;
				if (c1!=c2) {
					Kn++;
					if (c1+c2==1 || c1+c2==5) nc_ts++;
					else nc_tv++;
				}

				if (c1 == 1 || c1 == 3) nc_GC++;
				if (c2 == 1 || c2 == 3) nc_GC++;
			}
		}
		//cout << ts+tv << " / " << seq1.length() << "  ()"<<len<<"= " << (ts+tv) / len << "\t";
		//cout << " corrected JC69: " << getDistanceJC69((ts + tv)/len) << "\t";

		scaleArray(0.5/ nc_len, pi, 4);
		nc_GC = 0.5*nc_GC/nc_len; 
		Kn = getDistanceHKY(nc_ts/nc_len, nc_tv/nc_len, pi[2], pi[0], pi[3], pi[1]);
		//cout << " corrected HKY: " << Kn << "\t";
	}
	catch (...) {
		flag = false;
	}

	return flag;
}

/************************************************
 * Function: showParaInfo
 * Input Parameter:
 * Output: print on display
 * Return Value:
 *************************************************/
void KNKS::showParaInfo() {

	programInfo();

	//Input parameters
	cout << "Input noncoding file: " << input_nc_filename << endl;
	if (mutation_rate==NA) cout << "Input coding file: " << input_coding_filename << endl;
	else cout << "Input neutral mutation rate or Ks: " << mutation_rate << endl;
	
	//Output file(s)
	if (output_coding_filename.length() > 0) {
		cout << "Output files: " << output_nc_filename << "\t" << output_coding_filename << endl;
	}
	else {
		cout << "Output file: " << output_nc_filename << endl;
	}

	cout << "Genetic code: " << transl_table[2 * (genetic_code - 1) + 1] << endl;
	cout << "Please wait while reading sequences and calculating..." << endl;
}

/************************************************
 * Function: getTitleInfo
 * Input Parameter:
 * Output: get title information of outputing file
 * Return Value: string
 *************************************************/
string KNKS::getNCTitleInfo() {

	string title = "";
	int i;
	//Add "\t" to items except the last one
	for (i=0; i < titleInfo.size() - 1; i++) addString(title, titleInfo[i]);
	//The last item is added by "\n"
	addString(title, titleInfo[i], "\n");

	return title;
}

string KNKS::getCDSTitleInfo() {

	string title = "";
	string items[] = { "Sequence", "Method", "Ka", "Ks", "Ka/Ks",
		"P-Value(Fisher)", "Length", "S-Sites", "N-Sites", "Fold-Sites(0:2:4)",
		"Substitutions", "Syn-Subs", "Nonsyn-Subs", "Fold-Syn-Subs(0:2:4)", "Fold-Nonsyn-Subs(0:2:4)",
		"Divergence-Distance", "Substitution-Rate-Ratio(rTC:rAG:rTA:rCG:rTG:rCA/rCA)", "GC(1:2:3)", "ML-Score", "AICc",
		"Akaike-Weight", "Model" }; //21 outputing items

	//Format of outputing parameters
	int i;
	for (i = 0; i < sizeof(items) / sizeof(string) - 1; i++) addString(title, items[i]);

	addString(title, items[i], "\n");

	return title;
}

/***********************************
 * Function: programInfo
 * Input Parameter:
 * Output: Display program infomation.
 * Return Value: void
 ************************************/
void KNKS::programInfo() {
	cout << "********************************************************************" << endl;
	cout << "  Package: " << PACKAGE_NAME << " (" << VERSION << ")" << endl;
	cout << "  Program: " << KNKS_NAME << endl;
	cout << "  Description: " << KNKS_DESC << endl;
	cout << "  Reference: " << PACKAGE_REF << endl;
	cout << "********************************************************************" << endl;
}

/***********************************
 * Function: helpInfo
 * Input Parameter:
 * Output: Display help infomation.
 * Return Value: void
 ************************************/
void KNKS::helpInfo() {
	int i, j;

	programInfo();
	cout << endl;

	cout << "Usage: " << KNKS_NAME << " <options>" << endl;
	cout << "\t-i\tInput axt file name of non-coding aligned sequences [string]" << endl;
	cout << "\t-j\tInput axt file name of adjacent coding aligned sequences [string] || neutral mutation rate or Ks [numeric]" << endl;
	cout << "\t-o\tOutput file name for noncoding estimates [string]" << endl;
	cout << "\t-d\tWhether to output details of coding estimates [int, 0 or 1; default = 1)" << endl;
	cout << "\t-c\tGenetic code table [int, 1~33; default = 1-Standard Code)." << endl;
	for (i = 0, j = 0; i < NNCODE; i += 2) {
		if (strlen(transl_table[i]) > 0) {
			cout << "\t\t  " << transl_table[i + 1] << "\t";
			j++;
		}
		if (j == 2) {
			cout << endl;
			j = 0;
		}
	}
	cout << endl;
	cout << "\t\t  (More information about the Genetic Codes: http://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi?mode=c)" << endl;
	
	
	cout << "\t-h\tHelp information" << endl;
	cout << endl;

	cout << "Example:" << endl;
	cout << "\t" << KNKS_NAME << " -i test.axt -j adj.axt -o test.axt.knks\t//use 'adj.axt' to deduce neutral mutation rate" << endl;
	cout << "\t" << KNKS_NAME << " -i test.axt -j 0.618   -o test.axt.knks\t//use 0.618 as netural mutation rate" << endl;
	
	cout << endl;

	cout << "Please cite:" << endl;
	cout << "\t" << PACKAGE_REF << endl;
	cout << "See also: " << endl;
	cout << "\t" << "https://ngdc.cncb.ac.cn/biocode/tools/BT000001" << endl;
	cout << "Contact information:" << endl;
	cout << "\t" << "Please send bugs or advice to " << AUTHOR_NAME << " at " << AUTHOR_MAIL << "." << endl;
	cout << endl;
}