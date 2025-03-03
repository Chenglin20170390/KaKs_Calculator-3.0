/************************************************************
 * Copyright (C) CNCB-NGDC, BIG, CAS
 * All rights reserved.
 
 * Filename: KaKs.cpp
 * Abstract: Declaration of KAKS class including several methods.

 * Version: 1.0
 * Author: Zhang Zhang (zhanghzhang@genomics.org.cn)
 * Date: Feb.2, 2005

 * Version: 3.0
 * Author: Zhang Zhang (zhanghzhang@big.ac.cn)
 * Date: Jan.17, 2021
 *************************************************************/

#include "KaKs.h"

KAKS::KAKS() {

    string items[] = {"Sequence", "Method", "Ka", "Ks", "Ka/Ks",
        "P-Value(Fisher)", "Length", "S-Sites", "N-Sites", "Fold-Sites(0:2:4)",
        "Substitutions", "Syn-Subs", "Nonsyn-Subs", "Fold-Syn-Subs(0:2:4)", "Fold-Nonsyn-Subs(0:2:4)",
        "Divergence-Distance", "Substitution-Rate-Ratio(rTC:rAG:rTA:rCG:rTG:rCA/rCA)", "GC(1:2:3)", "ML-Score", "AICc",
        "Akaike-Weight", "Model"}; //21 outputing items

    //Format of outputing parameters
    int i;
    for (i = 0; i<sizeof (items) / sizeof (string); i++) titleInfo.push_back(items[i]);

    //Load Methods' Names and References for -h in linux, also for windows' tool tip
    method_name.push_back("NG");
    method_ref.push_back("Nei M and Gojobori T (1986) Mol Biol Evol");

    method_name.push_back("LWL");
    method_ref.push_back("Li WH, Wu CI and Luo CC (1985) Mol Biol Evol");

    method_name.push_back("LPB");
    method_ref.push_back("Li WH (1993) J Mol Evol & Pamilo P and Bianchi NO (1993) Mol Biol Evol");

    method_name.push_back("MLWL");
    method_ref.push_back("Tzeng YH, Pan R and Li WH (2004) Mol Biol Evol");

    method_name.push_back("MLPB");
    method_ref.push_back("Tzeng YH, Pan R and Li WH (2004) Mol Biol Evol");

    method_name.push_back("GY");
    method_ref.push_back("Goldman N and Yang Z (1994) Mol Biol Evol");

    method_name.push_back("YN");
    method_ref.push_back("Yang Z and Nielsen R (2000) Mol Biol Evol");

    method_name.push_back("MYN");
    method_ref.push_back("Zhang Z, Li J and Yu J (2006) BMC Evolutionary Biology");

    method_name.push_back("MS");
    method_ref.push_back("Zhang Z, et al (2006) Geno Prot Bioinfo");

    method_name.push_back("MA");
    method_ref.push_back("Zhang Z, et al (2006) Geno Prot Bioinfo");
    
    //method_name.push_back("ZZ");
    //method_ref.push_back("Zhang Z (2021) in preparation");

    Initialize();

}

KAKS::~KAKS() {

    Uninitialize();

    //Free memory
    titleInfo.clear();
    method_name.clear();
    method_ref.clear();
}

int KAKS::Initialize() {

    none = ng86 = lpb93 = lwl85 = mlwl85 = mlpb93 = yn00 = gy94 = myn06 = ms06 = ma06 = false;
    result4Win = result = seq_name = seq1 = seq2 = "";
    seq_filename = output_filename = detail_filename = "";
    result = details = "";
    genetic_code = 1;
    number = 0;

    return 1;
}

int KAKS::Uninitialize() {

    if (os.is_open()) {
        os.close();
    }

	hh = mm = ss = 0;

    return 1;
}

string KAKS::getResult4Win() {
	if (isOK4Win) {
		return result4Win;
	}
	return "";
}

/****************************************************
 * Function: ReadCalculateSeq
 * Input Parameter: string
 * Output: Read sequences, check sequences' validity
                  and calculate Ka and Ks.
 * Return Value: True if succeed, otherwise false.
 
 * Note: Using axt file for low memory
 *****************************************************/
bool KAKS::ReadCalculateSeq(string filename) {

	//Record the start time
	time_t time_start = time(NULL);

    bool flag = true;

    try {
		ifstream is(filename.c_str());
		if (!is) {
			cout << "Error in opening file..." << endl;
			throw 1;
		}

		showParaInfo();	//Show information on display

		result = getTitleInfo();

		//Output stream
		if (output_filename != "" && output_filename.length() > 0) {
			os.open(output_filename.c_str());
		}

		string temp = "", name = "", str = "";

		while (getline(is, temp, '\n')) {
			name = temp;
			str = "";

			getline(is, temp, '\n');
			while (temp != "") {
				str += temp;
				getline(is, temp, '\n');
			}

			//string msg = "";
			//if (checkPairwiseCoding(str, msg)) {
			if (checkValid(name, str)) {
				cout << "[" << ++number << "] " << name << "\t";
				bool isOK = calculateKaKs();
				if (isOK == false) {
					cout << "[Error in calculating]";
					throw 1;
				}
				else {
					cout << "[OK]";
				}
			}
			//else {
			//	cout << msg.c_str();
			//}
			cout << endl;
		}
		is.close();
		is.clear();
    } catch (...) {
        flag = false;
    }

	//Time used for running
	time_t t = time(NULL) - time_start;
	hh = t / 3600;
	mm = (t % 3600) / 60;
	ss = t - (t / 60) * 60;

    return flag;
}

/**************************************************
 * Function: checkValid
 * Input Parameter: string, string
 * Output: Check validity of pairwise sequences
 * Return Value: True if succeed, otherwise false. 
 ***************************************************/
bool KAKS::checkValid(string name, string str) {
    bool flag = true;
    long i;

    try {

		//Check whether sequences are equal in length
		if (str.length() % 2 !=0) {
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
    } catch (...) {
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
bool KAKS::Run(int argc, const char* argv[]) {

    bool flag = true;

    try {
        //Judge whether input parameters are legal
        if (!parseParameter(argc, argv)) {
            throw 1;
        }

        //Read sequences and calculate Ka & Ks
		ReadCalculateSeq(seq_filename);
        
		//Output results
		cout << "Outputing results: ";
		cout << output_filename;
        //Write details for model-selected method
		if ( writeFile(detail_filename, (getTitleInfo() + details).c_str())==true) {
			cout << "\t" << detail_filename;
		}
		cout << endl;

        //Print on display
        cout << "Mission accomplished. (Time elapsed: ";
        if (hh) cout << hh << ":" << mm << ":" << ss << ")" << endl;
        else cout << mm << ":" << ss << ")" << endl;
    } catch (...) {
        flag = false;
    }

    return flag;
}

/**************************************************
 * Function: parseParameter
 * Input Parameter: int, const char* []
 * Output: Parse the input parameters
 * Return Value: bool 
 ***************************************************/
bool KAKS::parseParameter(int argc, const char* argv[]) {

    bool flag = true;
    int i;
    string temp;

    try {

        //No parameter
        if (argc == 1) {
            if (seq_filename == "") throw 1;
            //else for only one parameter - seq_filename
        } else if (argc == 2) {//help information
            temp = argv[1];
            if (temp == "-H" || temp == "-h") {
                helpInfo();
                flag = false;
            } else {
                throw 1;
            }
        } else {
            //cout<<"test"<<endl;
            //parse parameters
            int inputflag = 0, outputflag = 0, codeflag = 0;
            for (i = 1; i < argc; i++) {

                temp = stringtoUpper(argv[i]);
                //Input axt file
                if (temp == "-I") {
                    if ((i + 1) < argc && inputflag == 0) {
                        seq_filename = argv[++i];
                        inputflag++;
                    } else {
                        throw 1;
                    }
                }//Output file
                else if (temp == "-O") {
                    if ((i + 1) < argc && outputflag == 0) {
                        output_filename = argv[++i];
                        outputflag++;
                    } else {
                        throw 1;
                    }
                }//Genetic Code Table
                else if (temp == "-C") {
                    if ((i + 1) < argc && codeflag == 0) {
                        genetic_code = CONVERT<int>(argv[++i]);
                        if (genetic_code < 1 || genetic_code > NCODE || strlen(transl_table[2 * (genetic_code - 1)]) < 1)
                            throw 1;
                        codeflag++;
                    } else {
                        throw 1;
                    }
                }//Details for Model Selection
                else if (temp == "-D") {
                    if ((i + 1) > argc) throw 1;
                    detail_filename = argv[++i];

                }//Algorithm(s) selected
                else if (temp == "-M") {
                    if ((i + 1) > argc) throw 1;
                    temp = stringtoUpper(argv[++i]);
                    if (temp == "NONE") none = true;
                    else if (temp == "NG") ng86 = true;
                    else if (temp == "LWL") lwl85 = true;
                    else if (temp == "LPB") lpb93 = true;
                    else if (temp == "MLPB") mlpb93 = true;
                    else if (temp == "MLWL") mlwl85 = true;
                    else if (temp == "GY") gy94 = true;
                    else if (temp == "YN") yn00 = true;
                    else if (temp == "MYN") myn06 = true;
                    else if (temp == "MS") ms06 = true;
                    else if (temp == "MA") ma06 = true;
                    else if (temp == "ALL") {
                        ng86 = lpb93 = lwl85 = mlwl85 = mlpb93 = gy94 = yn00 = myn06 = ms06 = ma06 = true;
                    } else throw 1;
                } else throw 1;
            }

            //If no input or output file, report error
            if (inputflag == 0 || outputflag == 0) throw 1;

            //Default: use ma to to calculate Ka and Ks
            if (!(none + ng86 + lpb93 + lwl85 + mlwl85 + mlpb93 + gy94 + yn00 + myn06 + ms06 + ma06)) {
				ma06 = true;
            }
        }
    } catch (...) {
        cout << "Input parameter(s) error." << endl;
        cout << "For help information: "<< KAKS_NAME <<" -h" << endl;
        flag = false;
    }

    return flag;
}

/*******************************************************
 * Function: Calculate
 * Input Parameter: void
 * Output: Calculate kaks and output results.
 * Return Value: void
 *
 * Note: 
 ********************************************************/
bool KAKS::calculateKaKs() {

    bool flag = true;

    try {

		//Get GCC at three codon positions
		getGCContent(seq1 + seq2);

        //Estimate Ka and Ks
        if (none) start_NONE();
        if (ng86) start_NG86();
        if (lwl85) start_LWL85();
        if (mlwl85) start_MLWL85();
        if (lpb93) start_LPB93();
        if (mlpb93) start_MLPB93();
        if (gy94) start_GY94();
        if (yn00) start_YN00();
        if (myn06) start_MYN();
        if (ms06 || ma06) start_MSMA();
        
		//add a lock "isOK4Win" to avoid the program collapse 
		isOK4Win = false;
        result4Win += result;
		isOK4Win = true;

        //Write into the file
        if (output_filename.length() > 0 && os.is_open()) {
            os << result;
            os.flush();
        }
        result = "";
    } catch (...) {
        flag = false;
    }

    return flag;
}

//NONE: NG without correction for multiple substitution

void KAKS::start_NONE() {

    NONE zz;
    result += zz.Run(seq1, seq2);
}

//NG

void KAKS::start_NG86() {

    NG86 zz;
    result += zz.Run(seq1, seq2);
}

//LWL

void KAKS::start_LWL85() {

    LWL85 zz;
    result += zz.Run(seq1, seq2);
}

//MLWL

void KAKS::start_MLWL85() {

    MLWL85 zz;
    result += zz.Run(seq1, seq2);
}

//LPB

void KAKS::start_LPB93() {

    LPB93 zz;
    result += zz.Run(seq1, seq2);
}

//MLPB

void KAKS::start_MLPB93() {

    MLPB93 zz;
    result += zz.Run(seq1, seq2);
}

//GY

void KAKS::start_GY94() {

    GY94 zz("HKY");
    result += zz.Run(seq1.c_str(), seq2.c_str());
}

//YN

void KAKS::start_YN00() {

    YN00 zz;
    result += zz.Run(seq1, seq2);
}

//MYN

void KAKS::start_MYN() {

    MYN zz;
    result += zz.Run(seq1, seq2);
}

/************************************************
 * Function: start_MSMA
 * Input Parameter: void
 * Output: Calculate Ka and Ks using the method of 
                  model selection or model averaging.
 * Return Value: void
 *************************************************/
void KAKS::start_MSMA() {

    vector<MLResult> result4MA; //generated by MS and used by MA

    //Model Selection
    MS zz1;
    string tmp = zz1.Run(seq1.c_str(), seq2.c_str(), result4MA, details);
    if (ms06) {
        result += tmp;
    }

    //Model Averaging
    if (ma06) {
        MA zz2;
        result += zz2.Run(seq1.c_str(), seq2.c_str(), result4MA);
    }
}


/************************************************
 * Function: showParaInfo
 * Input Parameter: 
 * Output: print on display
 * Return Value: 
 *************************************************/
void KAKS::showParaInfo() {
	
	programInfo();

	cout << "Input file: "<< seq_filename << endl;
	
	//Output file(s)
	if (detail_filename.length() > 0) {
		cout << "Output files: " << output_filename << "\t" << detail_filename <<endl;
	}
	else {
		cout << "Output file: " << output_filename << endl;
	}
    
	//Methods
	cout << "Method(s): ";
    if (none) cout << "NONE" << " ";
    if (ng86) cout << "NG" << " ";
    if (lwl85) cout << "LWL" << " ";
    if (mlwl85) cout << "MLWL" << " ";
    if (lpb93) cout << "LPB" << " ";
    if (mlpb93) cout << "MLPB" << " ";
    if (gy94) cout << "GY" << " ";
    if (yn00) cout << "YN" << " ";
    if (myn06) cout << "MYN" << " ";
    if (ms06) cout << "MS" << " ";
    if (ma06) cout << "MA" << " ";
	cout << endl;

	cout << "Genetic code: " << transl_table[2 * (genetic_code - 1) + 1] << endl;
    cout << "Please wait while reading sequences and calculating..." << endl;
}

/************************************************
 * Function: getTitleInfo
 * Input Parameter: 
 * Output: get title information of outputing file
 * Return Value: string
 *************************************************/
string KAKS::getTitleInfo() {

    string title = "";
    int i = 0;

    if (titleInfo.size() > 0) {
        //Add "\t" to items except the last one
        for (; i < titleInfo.size() - 1; i++) addString(title, titleInfo[i]);
        //The last item is added by "\n"
        addString(title, titleInfo[i], "\n");
    }

    return title;
}

/***********************************
 * Function: programInfo
 * Input Parameter:
 * Output: Display program infomation.
 * Return Value: void
 ************************************/
void KAKS::programInfo() {
	cout << "********************************************************************" << endl;
	cout << "  Package: " << PACKAGE_NAME << " (" << VERSION << ")"<<endl;
	cout << "  Program: " << KAKS_NAME << endl;
	cout << "  Description: " << KAKS_DESC << endl;
	cout << "  Reference: " << PACKAGE_REF << endl;
	cout << "********************************************************************" << endl;
}

/***********************************
 * Function: helpInfo
 * Input Parameter: 
 * Output: Display help infomation.
 * Return Value: void
 ************************************/
void KAKS::helpInfo() {
    int i, j;

	programInfo();
	cout << endl;

    cout << "Usage: "<< KAKS_NAME <<" <options>" << endl;
    cout << "\t-i\tInput axt file name, containing coding aligned sequences [string]" << endl;
    cout << "\t-o\tOutput file name for saving results [string]" << endl;
    cout << "\t-c\tGenetic code table [int, 1~33; default = 1-Standard Code]" << endl;
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

    cout << "\t-m\tMethods for estimating Ka and Ks and theirs references [string, default = MA]" << endl;
    for (i = 0; i < method_name.size(); i++) {
        cout << "\t\t  " << method_name[i].c_str();
        cout << "\t\t" << method_ref[i].c_str() << endl;
    }
    cout << "\t\t  ALL(including all above methods)" << endl;

    cout << "\t-d\tFile name for details about each candidate model when using the method of MS" << endl;
	cout << "\t-h\tHelp information" << endl; 
	cout << endl;

    cout << "Example:" << endl;
    cout << "\t" << KAKS_NAME << " -i test.axt -o test.axt.kaks\t//use MA method based on a more suitable model and standard code" << endl;
    cout << "\t" << KAKS_NAME << " -i test.axt -o test.axt.kaks -c 2\t//use MA method and vertebrate mitochondrial code" << endl;
    cout << "\t" << KAKS_NAME << " -i test.axt -o test.axt.kaks -m LWL -m MYN\t//use LWL and MYN methods, and standard Code" << endl;

    cout << endl;

	cout << "Please cite:" << endl;
	cout << "\t" << PACKAGE_REF << endl;
	cout << "See also: "<< endl;
	cout << "\t" << "https://bigd.big.ac.cn/biocode/tools/BT000001" << endl;
	cout << "Contact information:" << endl;
    cout << "\t" << "Please send bugs or advice to "<<AUTHOR_NAME<<" at "<<AUTHOR_MAIL<<"."<<endl;
    cout << endl;
}

