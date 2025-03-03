/************************************************************
* Copyright (C) 2005, BIG of Chinese Academy of Sciences
* All rights reserved.
 
* Filename: base.cpp
* Abstract: Definition of base class for KaKs methods.

* Version: 1.0
* Author: Zhang Zhang (zhanghzhang@genomics.org.cn)
* Date: Feb.2, 2005

*************************************************************/

#include "base.h"


/******** Global variables ********/
/*						The Genetic Codes 
http://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi?mode=c
	Last update of the Genetic Codes: Jan 28, 2021 */
int genetic_code=1; //from 1 to 33
/* Genetic standard codon table, !=stop codon */
const char* transl_table[] = {
 "FFLLSSSSYY!!CC!WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "1-Standard Code",
 "FFLLSSSSYY!!CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSS!!VVVVAAAADDEEGGGG", "2-Vertebrate Mitochondrial Code",
 "FFLLSSSSYY!!CCWWTTTTPPPPHHQQRRRRIIMMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "3-Yeast Mitochondrial Code",
 "FFLLSSSSYY!!CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "4-Mold, Protozoan, and Coelenterate Mitochondrial Code and the Mycoplasma/Spiroplasma Code",
 "FFLLSSSSYY!!CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSSSVVVVAAAADDEEGGGG", "5-Invertebrate Mitochondrial Code",
 "FFLLSSSSYYQQCC!WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "6-Ciliate, Dasycladacean and Hexamita Nuclear Code",
 "", "7-",
 "", "8-",
 "FFLLSSSSYY!!CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG", "9-Echinoderm and Flatworm Mitochondrial Code",
 "FFLLSSSSYY!!CCCWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "10-Euplotid Nuclear Code",
 "FFLLSSSSYY!!CC!WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "11-Bacterial, Archaeal and Plant Plastid Code",
 "FFLLSSSSYY!!CC!WLLLSPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "12-Alternative Yeast Nuclear Code",
 "FFLLSSSSYY!!CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSGGVVVVAAAADDEEGGGG", "13-Ascidian Mitochondrial Code",
 "FFLLSSSSYYY!CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG", "14-Alternative Flatworm Mitochondrial Code",
 "", "15-",
 "FFLLSSSSYY!LCC!WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "16-Chlorophycean Mitochondrial Code",
 "", "17-",
 "", "18-",
 "", "19-",
 "", "20-",
 "FFLLSSSSYY!!CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNNKSSSSVVVVAAAADDEEGGGG", "21-Trematode Mitochondrial Code",
 "FFLLSS!SYY!LCC!WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "22-Scenedesmus obliquus Mitochondrial Code",
 "FF!LSSSSYY!!CC!WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "23-Thraustochytrium Mitochondrial Code",
 "FFLLSSSSYY!!CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSSKVVVVAAAADDEEGGGG", "24-Rhabdopleuridae Mitochondrial Code",
 "FFLLSSSSYY!!CCGWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "25-Candidate Division SR1 and Gracilibacteria Code",
 "FFLLSSSSYY!!CC!WLLLAPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "26-Pachysolen tannophilus Nuclear Code",
 "FFLLSSSSYYQQCCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "27-Karyorelict Nuclear Code",
 "FFLLSSSSYYQQCCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "28-Condylostoma Nuclear Code",
 "FFLLSSSSYYYYCC!WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "29-Mesodinium Nuclear Code",
 "FFLLSSSSYYEECC!WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "30-Peritrich Nuclear Code",
 "FFLLSSSSYYEECCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "31-Blastocrithidia Nuclear Code",
 "", "32-",
 "FFLLSSSSYYY!CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSSKVVVVAAAADDEEGGGG", "33-Cephalodiscidae Mitochondrial UAA-Tyr Code"
};

string seq_name;		//sequences' name
unsigned long length;	//sequences' length
double GC[4];			//GC Content
//********End of Global variables**********

//Constructor function
Base::Base() {

	//Codons
	int i, id = 0;
    string nuc = "TCAG";
	for (i=0; i<64; i++) {
            ID2Codon[i] = getCodon(i);
        }
            
	IntString::iterator pos;
	for(pos=ID2Codon.begin(); pos!=ID2Codon.end(); pos++) {
		//cout<<pos->first << " "<<pos->second<<endl;
		Codon2ID[pos->second] = pos->first;
	}
        
        //Codons & Genetic Codes
	Codon64["AAA"] = Codon("D", "1", "KKKKKK!!NKKKKNKK!!!!NKK");
	Codon64["AAG"] = Codon("D", "1", "KKKKKK!!KKKKKKKK!!!!KKK");
	Codon64["AAT"] = Codon("D", "1", "NNNNNN!!NNNNNNNN!!!!NNN");
	Codon64["AAC"] = Codon("D", "1", "NNNNNN!!NNNNNNNN!!!!NNN");
	Codon64["TAA"] = Codon("D", "1", "*****Q!!*****Y**!!!!***");
	Codon64["TAG"] = Codon("D", "1", "*****Q!!******QL!!!!*L*");
	Codon64["TAT"] = Codon("D", "1", "YYYYYY!!YYYYYYYY!!!!YYY");
	Codon64["TAC"] = Codon("D", "1", "YYYYYY!!YYYYYYYY!!!!YYY");
	Codon64["ATA"] = Codon("D", "1", "IMMIMI!!IIIIMIII!!!!MII");
	Codon64["ATG"] = Codon("D", "1", "MMMMMM!!MMMMMMMM!!!!MMM");
	Codon64["ATT"] = Codon("D", "1", "IIIIII!!IIIIIIII!!!!III");
	Codon64["ATC"] = Codon("D", "1", "IIIIII!!IIIIIIII!!!!III");
	Codon64["TTA"] = Codon("D", "1", "LLLLLL!!LLLLLLLL!!!!LL*");
	Codon64["TTG"] = Codon("D", "1", "LLLLLL!!LLLLLLLL!!!!LLL");
	Codon64["TTT"] = Codon("D", "1", "FFFFFF!!FFFFFFFF!!!!FFF");
	Codon64["TTC"] = Codon("D", "1", "FFFFFF!!FFFFFFFF!!!!FFF");

	Codon64["GAA"] = Codon("D", "2", "EEEEEE!!EEEEEEEE!!!!EEE");
	Codon64["GAG"] = Codon("D", "2", "EEEEEE!!EEEEEEEE!!!!EEE");
	Codon64["GAT"] = Codon("D", "2", "DDDDDD!!DDDDDDDD!!!!DDD");
	Codon64["GAC"] = Codon("D", "2", "DDDDDD!!DDDDDDDD!!!!DDD");
	Codon64["CAA"] = Codon("D", "2", "QQQQQQ!!QQQQQQQQ!!!!QQQ");
	Codon64["CAG"] = Codon("D", "2", "QQQQQQ!!QQQQQQQQ!!!!QQQ");
	Codon64["CAT"] = Codon("D", "2", "HHHHHH!!HHHHHHHH!!!!HHH");
	Codon64["CAC"] = Codon("D", "2", "HHHHHH!!HHHHHHHH!!!!HHH");
	Codon64["GTA"] = Codon("R", "2", "VVVVVV!!VVVVVVVV!!!!VVV");
	Codon64["GTG"] = Codon("R", "2", "VVVVVV!!VVVVVVVV!!!!VVV");
	Codon64["GTT"] = Codon("R", "2", "VVVVVV!!VVVVVVVV!!!!VVV");
	Codon64["GTC"] = Codon("R", "2", "VVVVVV!!VVVVVVVV!!!!VVV");
	Codon64["CTA"] = Codon("R", "2", "LLTLLL!!LLLLLLLL!!!!LLL");
	Codon64["CTG"] = Codon("R", "2", "LLTLLL!!LLLSLLLL!!!!LLL");
	Codon64["CTT"] = Codon("R", "2", "LLTLLL!!LLLLLLLL!!!!LLL");
	Codon64["CTC"] = Codon("R", "2", "LLTLLL!!LLLLLLLL!!!!LLL");

	Codon64["AGA"] = Codon("D", "3", "R*RRSR!!SRRRGSRR!!!!SRR");
	Codon64["AGG"] = Codon("D", "3", "R*RRSR!!SRRRGSRR!!!!SRR");
	Codon64["AGT"] = Codon("D", "3", "SSSSSS!!SSSSSSSS!!!!SSS");
	Codon64["AGC"] = Codon("D", "3", "SSSSSS!!SSSSSSSS!!!!SSS");
	Codon64["TGA"] = Codon("D", "3", "*WWWW*!!WC**WW**!!!!W**");
	Codon64["TGG"] = Codon("D", "3", "WWWWWW!!WWWWWWWW!!!!WWW");
	Codon64["TGT"] = Codon("D", "3", "CCCCCC!!CCCCCCCC!!!!CCC");
	Codon64["TGC"] = Codon("D", "3", "CCCCCC!!CCCCCCCC!!!!CCC");
	Codon64["ACA"] = Codon("R", "3", "TTTTTT!!TTTTTTTT!!!!TTT");
	Codon64["ACG"] = Codon("R", "3", "TTTTTT!!TTTTTTTT!!!!TTT");
	Codon64["ACT"] = Codon("R", "3", "TTTTTT!!TTTTTTTT!!!!TTT");
	Codon64["ACC"] = Codon("R", "3", "TTTTTT!!TTTTTTTT!!!!TTT");
	Codon64["TCA"] = Codon("R", "3", "SSSSSS!!SSSSSSSS!!!!S*S");
	Codon64["TCG"] = Codon("R", "3", "SSSSSS!!SSSSSSSS!!!!SSS");
	Codon64["TCT"] = Codon("R", "3", "SSSSSS!!SSSSSSSS!!!!SSS");
	Codon64["TCC"] = Codon("R", "3", "SSSSSS!!SSSSSSSS!!!!SSS");

	Codon64["GGA"] = Codon("R", "4", "GGGGGG!!GGGGGGGG!!!!GGG");
	Codon64["GGG"] = Codon("R", "4", "GGGGGG!!GGGGGGGG!!!!GGG");
	Codon64["GGT"] = Codon("R", "4", "GGGGGG!!GGGGGGGG!!!!GGG");
	Codon64["GGC"] = Codon("R", "4", "GGGGGG!!GGGGGGGG!!!!GGG");
	Codon64["CGA"] = Codon("R", "4", "RRRRRR!!RRRRRRRR!!!!RRR");
	Codon64["CGG"] = Codon("R", "4", "RRRRRR!!RRRRRRRR!!!!RRR");
	Codon64["CGT"] = Codon("R", "4", "RRRRRR!!RRRRRRRR!!!!RRR");
	Codon64["CGC"] = Codon("R", "4", "RRRRRR!!RRRRRRRR!!!!RRR");
	Codon64["GCA"] = Codon("R", "4", "AAAAAA!!AAAAAAAA!!!!AAA");
	Codon64["GCG"] = Codon("R", "4", "AAAAAA!!AAAAAAAA!!!!AAA");
	Codon64["GCT"] = Codon("R", "4", "AAAAAA!!AAAAAAAA!!!!AAA");
	Codon64["GCC"] = Codon("R", "4", "AAAAAA!!AAAAAAAA!!!!AAA");
	Codon64["CCA"] = Codon("R", "4", "PPPPPP!!PPPPPPPP!!!!PPP");
	Codon64["CCG"] = Codon("R", "4", "PPPPPP!!PPPPPPPP!!!!PPP");
	Codon64["CCT"] = Codon("R", "4", "PPPPPP!!PPPPPPPP!!!!PPP");
	Codon64["CCC"] = Codon("R", "4", "PPPPPP!!PPPPPPPP!!!!PPP"); 
        

	for(i=0; i<5; i++) {
		Si[i] = Vi[i] = L[i] = NULL;
	}
	for (i=0; i<NUMBER_OF_RATES; i++) {
		KAPPA[i] = 1;
	}
	
	SEKa = SEKs = AICc = lnL = AkaikeWeight = NA;
	Ka = Ks = Sd = Nd = S = N = snp = t = kappa = NULL;

	model = "";
}

/********************************************
* Function: addString
* Input Parameter: string, string, string
* Output: result = result + str + flag
* Return Value: void
* Note: flag = "\t" (default) or "\n"
*********************************************/
void Base::addString(string &result, string str, string flag) {
	result += str;
	result += flag;
}

/**********************************************************************
* Function: getAminoAcid
* Input Parameter: codon or codon's id
* Output: Calculate the amino acid according to codon or codon's id.
* Return Value: char 
***********************************************************************/
char Base::getAminoAcid(string codon) {
	return transl_table[2*(genetic_code-1)][getID(codon)];
}
char Base::getAminoAcid(int id) {
	return transl_table[2*(genetic_code-1)][id];
}

/**********************************
* Function: getNumNonsense
* Input Parameter: int
* Output: get the number of nonsense codons
* Return Value: int
***********************************/
int Base::getNumNonsense(int genetic_code) {

	int num, i;
	for(num=i=0; i<CODON; i++) {
		if(getAminoAcid(i)=='!') num++;
	}

	return num;
}

/********************************************
* Function: getID
* Input Parameter: codon
* Output: Get codon's id in array of codon_table.
* Return Value: int
*********************************************/
int Base::getID(string codon) {
	return (convertChar(codon[0])*XSIZE + convertChar(codon[1])*DNASIZE + convertChar(codon[2]));
}

/********************************************
* Function: getCodon
* Input Parameter: int
* Output: Get the codon according to id;
		  a reverse funtion of getID.
* Return Value: string
*********************************************/
string Base::getCodon(int IDcodon) {
	
	string codon = "TTT";

	if (IDcodon>=0 && IDcodon<64) {
		codon[0]=convertInt(IDcodon/16); 
		codon[1]=convertInt((IDcodon%16)/4);
		codon[2]=convertInt(IDcodon%4);
	}

	return codon;
}

/*********************************************
* Function: convertChar
* Input Parameter: ch as char
* Output: Convert a char-T,C,A,G into a digit
*		  0,1,2,3, respectively.
* Return Value: int.
**********************************************/
int Base::convertChar(char ch) {
	int ret = -1;
	switch(ch) {
		case 'T':case 'U':
			ret = 0;
			break;
		case 'C':
			ret = 1;
			break;
		case 'A':
			ret = 2;
			break;
		case 'G':
			ret = 3;
			break;
	}
	return ret;
}

/********************************************
* Function: convertInt
* Input Parameter: int
* Output: Convert a digit- 0,1,2,3 into a 
*		  char-T,C,A,G, respectively.
* Return Value: char
*********************************************/
char Base::convertInt(int i) {
	char ch = '-';
	switch(i) {
		case 0:
			ch = 'T';
			break;
		case 1:
			ch = 'C';
			break;
		case 2:
			ch = 'A';
			break;
		case 3:
			ch = 'G';
			break;
	}
	return ch;
}

/********************************************
* Function: stringtoUpper
* Input Parameter: string
* Output: upper string
* Return Value: string
*********************************************/
string Base::stringtoUpper(string str) {
	int j;
	for(j=0; j<str.length(); j++) str[j] = toupper(str[j]);
	return str;
}

/********************************************
* Function: getRandom
* Input Parameter: void
* Output: Generate a radnom integer
* Return Value: int
*********************************************/
int Base::getRandom() {	
	srand((unsigned)time(NULL));
	return rand();
}

/********************************************
* Function: initArray
* Input Parameter: array of int/double, int, int/double(default=0)
* Output: Init the array x[0...n-1]=value
* Return Value: int
*********************************************/
int Base::initArray(double x[], int n, double value) {
	int i; 
	for(i=0; i<n; i++) x[i] = value;
	return 0;
}

int Base::initArray(int x[], int n, int value) {
	int i; 
	for(i=0; i<n; i++) x[i] = value;
	return 0;
}

/********************************************
* Function: sumArray
* Input Parameter: double/int, int, int(default=0)
* Output: Sum of array x[]
* Return Value: double/int
*********************************************/
double Base::sumArray(double x[], int end, int begin) { 
	int i;
	double sum=0.;	
	for(i=begin; i<end; i++) sum += x[i];
	return sum;
}

int Base::sumArray(int x[], int end, int begin) { 
	int i, sum=0;	
	for(i=begin; i<end; i++) sum += x[i];
	return sum;
}

/********************************************
* Function: norm
* Input Parameter: array of double, int
* Output: Sqrt of the sum of the elements' square 
           sqrt(x0*x0 + x1*x1 + ...)
* Return Value: double
*********************************************/
double Base::norm(double x[], int n) {
	int i; 
	double t=0; 

	for(i=0; i<n; i++) t += square(x[i]);

	return sqrt(t);
}

/********************************************
* Function: scaleArray
* Input Parameter: double, array of double, int
* Output: Elements in array are mutipled by scale 
* Return Value: int
*********************************************/
int Base::scaleArray(double scale, double x[], int n) {	
	int i; 	
	for (i=0; i<n; i++) x[i] *= scale;

	return 1; 
}

/********************************************
* Function: copyArray
* Input Parameter: array, array, int
* Output: Copy array's values one by one: to[] = from[]
* Return Value: int
*********************************************/
int Base::copyArray(double from[], double to[], int n) {	
	int i; 
	for(i=0; i<n; i++) to[i] = from[i];
	
	return 1; 
}

/********************************************
* Function: innerp
* Input Parameter: array, array, int
* Output: Sum of 'n' products multiplied by 
			two elements in x[], y[].
* Return Value: int
*********************************************/
double Base::innerp(double x[], double y[], int n) {
	
	int i; 
	double t=0;

	for(i=0; i<n; t += x[i]*y[i], i++); 

	return t; 
}

/********************************************
* Function: initIdentityMatrix
* Input Parameter: array of double, int
* Output: Set x[i,j]=0 when x!=j and 
			  x[i,j]=1 when x=j 
* Return Value: int
*********************************************/
int Base::initIdentityMatrix(double x[], int n) {
	
	int i,j;

	for (i=0; i<n; i++)  {
		for(j=0; j<n; x[i*n+j]=0, j++);
		x[i*n+i] = 1; 
	} 

	return 0; 
}


/************************************************
* Function: writeFile
* Input Parameter: string, string
* Output: Write content into the given file.
* Return Value: True if succeed, otherwise false. 
*************************************************/
bool Base::writeFile(string output_filename, const char* result) {
	
	bool flag = true;
	try {
		//file name is ok
		if (output_filename!="" && output_filename.length()>0) {
			ofstream os(output_filename.c_str());
			if (!os.is_open()) throw 1;

			os<<result;
			os.close();					
		}
	}
	catch (...) {
		cout<<"Error in writing to file...";
		flag = false;
	}	

	return flag;
}

/* Get GCC of entire sequences GC[0] and of three codon positions GC[1,2,3] */
void Base::getGCContent(string str, int cds) {
	int i, j;

	initArray(GC, 4);
	if (cds == 1) {
		for (i = 0; i < str.length(); i += 3) {
			string codon = str.substr(i, 3);
			for (j = 0; j < 3; j++) {
				if (codon[j] == 'G' || codon[j] == 'C') GC[j + 1]++;
			}
		}
		GC[0] = sumArray(GC, 4, 1) / str.length()*1.0;
		for (i = 1; i < 4; i++) GC[i] /= (str.length() / 3.0);
	}
	else {
		for (i = 0; i < str.length(); i++) {
			if (str[i] == 'G' || str[i] == 'C') GC[0]++;
		}
		GC[0] /= 1.0* str.length();
	}

	return;
}

/****************************************************
 * Function: readAXTSeq
 * Input Parameter: string
 * Output: Read sequences and their names.
 * Return Value: True if succeed, otherwise false.

 * Note: Using axt file
 *****************************************************/
bool Base::readAXTSeq(string filename, vector<string> &vec_names, vector<string> &vec_seqs) {

	bool flag = true;
	vec_names.clear();
	vec_seqs.clear();

	try {

		ifstream is(filename.c_str());
		if (!is) {
			cout << "Error in opening file..." << endl;
			throw 1;
		}

		string temp = "";
		int count = 0;
		while (getline(is, temp, '\n')) {

			string name = "", seq = "";

			name = temp;

			getline(is, temp, '\n');
			while (temp != "") {
				seq += temp;
				getline(is, temp, '\n');
			}

			vec_names.push_back(name);
			vec_seqs.push_back(seq);
		}
		is.close();
		is.clear();

	}
	catch (...) {
		flag = false;
	}

	return flag;
}

/**************************************************
 * Function: checkPairwiseCoding
 * Input Parameter: string, string
 * Output: Check validity of pairwise sequences
 * Return Value: True if succeed, otherwise false.
 ***************************************************/
bool Base::checkPairwiseCoding(string &str, string &msg) {
	bool flag = true;
	try {

		//Check whether sequences are equal in length
		if (str.length() % 2 != 0) {
			msg = "[Error. The sequences are not equal in length.]";
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
		long i, j;
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
		}//end of for

		str = str1 + str2;
	}
	catch (...) {
		flag = false;
	}

	return flag;
}

bool Base::checkPairwiseNoncoding(string seq, string &msg) {
	bool flag = true;

	try {
		//Check whether pairwise sequences are equal in length
		if (seq.length() % 2 != 0) {
			msg = "[Error. The sequences are not equal in length.]";
			throw 1;
		}
	}
	catch (...) {
		flag = false;
	}
	return flag;
}

/*****************************************************
* Function: parseOutput
* Input Parameter: void
* Output: Parse estimated results for outputing
* Return Value: string

  Order: "Sequence", "Method", "Ka", "Ks", "Ka/Ks", 
		 "P-Value(Fisher)", "Length", "S-Sites", "N-Sites", "Fold-Sites(0:2:4)",
		 "Substitutions", "S-Substitutions", "N-Substitutions", "Fold-S-Substitutions(0:2:4)", "Fold-N-Substitutions(0:2:4)", 
		 "Divergence-Time", "Substitution-Rate-Ratio(rTC:rAG:rTA:rCG:rTG:rCA/rCA)", "GC(1:2:3)", "ML-Score", "AICc",
		 "Model"
******************************************************/		
string Base::parseOutput() {

	int i;
	string result = "", tmp;

	//Sequence name
	addString(result, seq_name);
	//Method name
	addString(result, name);

	//Ka
	if (Ka<SMALLVALUE) {
		tmp = "NA";
	}
	else {
		tmp = CONVERT<string>(Ka);
	}
	addString(result, tmp);
	
	//Ks
	if (Ks<SMALLVALUE) {
		tmp = "NA";
	}
	else {
		tmp = CONVERT<string>(Ks);
	}
	addString(result, tmp);
	
	//Ka/Ks
	if(Ks<SMALLVALUE || Ks==NA || Ka==NA) {
		tmp = "NA";
	}
	else {
		tmp = CONVERT<string>(Ka/Ks);
	}
	addString(result, tmp);

	//Fisher's test: p_value
	if (Sd<SMALLVALUE || Nd<SMALLVALUE || S<SMALLVALUE || N<SMALLVALUE)
		tmp = "NA";
	else {
		tmp = CONVERT<string>(fisher(Sd,Nd,S-Sd,N-Nd));
	}
	addString(result, tmp);

	//Length of compared pairwise sequences
	addString(result, CONVERT<string>(length));
	
	//Synonymous(S) sites
	if (S<SMALLVALUE) {
		tmp = "NA";
	}
	else {
		tmp = CONVERT<string>(S);
	}
	addString(result, tmp);

	//Nonsynonymous(N) sites
	if (N<SMALLVALUE) {
		tmp = "NA";
	}
	else {
		tmp = CONVERT<string>(N);
	}
	addString(result, tmp);

	//L[0], L[2], L[4] only for Prof.Li's series(LWL85, LPB93...)
	if (L[0]<SMALLVALUE && L[2]<SMALLVALUE && L[4]<SMALLVALUE) {
		tmp = "NA";		
	}
	else {		
		tmp = CONVERT<string>(L[0]);	tmp += ":";
		tmp += CONVERT<string>(L[2]);	tmp += ":";
		tmp += CONVERT<string>(L[4]);		
	}
	addString(result, tmp);

	
	//Substitutions
	addString(result, CONVERT<string>(snp));

	//Sysnonymous(Sd) Substitutions(Nd)	
	if (Sd>SMALLVALUE) {
		tmp = CONVERT<string>(Sd);		
	}	
	else {
		tmp = "NA";
	}
	addString(result, tmp);
	
	//Nonsysnonymous Substitutions(Nd)
	if (Nd>SMALLVALUE) {
		tmp = CONVERT<string>(Nd);		
	}	
	else {
		tmp = "NA";
	}
	addString(result, tmp);

	//Si for Li's series' methods(LWL85, LPB93...)
	if (Si[0]!=NULL || Si[2]!=NULL || Si[4]!=NULL) { //Si[0], Si[2], Si[4]
		tmp  = CONVERT<string>(Si[0]);	tmp += ":";
		tmp += CONVERT<string>(Si[2]);	tmp += ":";
		tmp += CONVERT<string>(Si[4]);		
	}
	else {
		tmp = "NA";
	}
	addString(result, tmp);

	//Vi for Li's series' methods(LWL85, LPB93...)
	if (Vi[0]!=NULL || Vi[2]!=NULL || Vi[4]!=NULL) { //Vi[0], Vi[2], Vi[4]
		tmp  = CONVERT<string>(Vi[0]);	tmp += ":";
		tmp += CONVERT<string>(Vi[2]);	tmp += ":";
		tmp += CONVERT<string>(Vi[4]);		
	}
	else {
		tmp = "NA";
	}
	addString(result, tmp);


	//Divergence time or distance t = (S*Ks+N*Ka)/(S+N)
	if(t<SMALLVALUE) {
		tmp = "NA";
	}
	else {
		tmp = CONVERT<string>(t);
	}
	addString(result, tmp);

	//Substitution-Rate-Ratio(rTC:rAG:rTA:rCG:rTG:rCA/rCA)
	for(i=0, tmp=""; i<NUMBER_OF_RATES-1; i++) {
		tmp += CONVERT<string>(KAPPA[i]); 
		tmp += ":";
	}
	tmp += CONVERT<string>(KAPPA[i]); 
	addString(result, tmp);

	//GC Content
	tmp  = CONVERT<string>(GC[0]);	tmp += "(";
	tmp += CONVERT<string>(GC[1]);	tmp += ":";
	tmp += CONVERT<string>(GC[2]);	tmp += ":";
	tmp += CONVERT<string>(GC[3]);	tmp += ")";
	addString(result, tmp);
	
	//Maximum Likelihood Value
	if (lnL==NA) tmp = "NA";
	else         tmp = CONVERT<string>(lnL);
	addString(result, tmp);
	
	//AICc
	if (AICc==NA) tmp = "NA";
	else		  tmp = CONVERT<string>(AICc);
	addString(result, tmp);

	//Akaike weight in model selection
	if (AkaikeWeight==NA) tmp = "NA";
	else		  tmp = CONVERT<string>(AkaikeWeight);
	addString(result, tmp);
	
	//Selected Model according to AICc
	if (model==""||model.length()==0) tmp = "NA";
	else tmp = model;
	addString(result, tmp, "\n");

/*
	//Standard Errors
	if (SEKa==NA) tmp = "NA";
	else tmp = CONVERT<string>(SEKa);
	addString(result, tmp, "\t");
	
	if (SEKs==NA) tmp = "NA";
	else tmp = CONVERT<string>(SEKs);
	addString(result, tmp, "\n");
*/
	
	return result;
}

/**************************************************
* Function: fisher
* Input Parameter: double, double, double, double
* Output: Compute p-value by Fisher exact test
* Return Value: double
***************************************************/
double Base::fisher(double sd, double nd, double s, double n) {
	double denominator, numerator, prob_total, prob_current, sum, fac_sum;
	double matrix[4],  R[2], C[2];
	int i, j;

	denominator = numerator = prob_total = prob_current = sum = fac_sum = 0.0;

	matrix[0]=sd;      matrix[1]=s;	    matrix[2]=nd;    matrix[3]=n;
	//Row & Column
	R[0]=matrix[0]+matrix[2];		R[1]=matrix[1]+matrix[3];
	C[0]=matrix[0]+matrix[1];		C[1]=matrix[2]+matrix[3];
	sum = R[0]+R[1];	
	
	//Calculate the numberator that is a constant
	numerator += factorial(R[0]);
	numerator += factorial(R[1]);
	numerator += factorial(C[0]);
	numerator += factorial(C[1]);
	
	//Log of Factorial of N
	fac_sum = factorial(sum);	
	for(i=0, denominator=fac_sum; i<4; i++) {
		denominator += factorial(matrix[i]);
	}
	//Probability of current situtation
	prob_current = exp(numerator-denominator);
	
	//Two-tail probabilities if less than prob_current
	for(i=0; i<=R[0]; i++) {
		matrix[0] = i;       
		matrix[1] = C[0]-i;
		matrix[2] = R[0]-i;
		matrix[3] = R[1]-C[0]+i;
		if (matrix[0]>=0 && matrix[1]>=0 && matrix[2]>=0 && matrix[3]>=0) {			
			for(j=0, denominator=fac_sum; j<4; j++) {				
				denominator += factorial(matrix[j]);
			}
			double temp = numerator-denominator;
			temp = exp(numerator-denominator);
			if (temp<=prob_current) {
				prob_total += temp;
			}
		}
	}
	
	return prob_total;
}

/**************************************************
* Function: factorial
* Input Parameter: n
* Output: Compute the factorial of 'n', then return 
          the log of it.
* Return Value: double
***************************************************/
double Base::factorial(double n) {
	double temp=1.0;
	if (n>0) {
		n = n + 1;
		double x = 0;
		x += 0.1659470187408462e-06/(n+7);
		x += 0.9934937113930748e-05/(n+6);
		x -= 0.1385710331296526    /(n+5);
		x += 12.50734324009056     /(n+4);
		x -= 176.6150291498386     /(n+3);
		x += 771.3234287757674     /(n+2);
		x -= 1259.139216722289     /(n+1);
		x += 676.5203681218835     /(n);
		x += 0.9999999999995183;
		temp = log(x) - 5.58106146679532777 - n +(n - 0.5)*log(n + 6.5);
	}
	
	return (temp);
}


/*

int matby (double a[], double b[], double c[], int n,int m,int k)
// a[n*m], b[m*k], c[n*k]  ......  c = a*b

{
   int i,j,i1;
   double t;
   FOR (i,n)  FOR(j,k) {
      for (i1=0,t=0; i1<m; i1++) t+=a[i*m+i1]*b[i1*k+j];
      c[i*k+j] = t;
   }
   return (0);
}





void PMatrixTaylor(double P[], double n, double t) {

// Get approximate PMat using polynomial approximation
//   P(t) = I + Qt + (Qt)^2/2 + (Qt)^3/3!

   int nterms=1000, i,j,k, c[2],ndiff,pos=0,from[3],to[3];
   double *Q=space, *Q1=Q+n*n, *Q2=Q1+n*n, mr, div;

   FOR (i,n*n) Q[i]=0;
   for (i=0; i<n; i++) FOR (j,i) {
      from[0]=i/16; from[1]=(i/4)%4; from[2]=i%4;
      to[0]=j/16;   to[1]=(j/4)%4;   to[2]=j%4;
      c[0]=GenetCode[com.icode][i];
      c[1]=GenetCode[com.icode][j];
      if (c[0]==-1 || c[1]==-1)  continue;
      for (k=0,ndiff=0; k<3; k++)  if (from[k]!=to[k]) { ndiff++; pos=k; }
      if (ndiff!=1)  continue;
      Q[i*n+j]=1;
      if ((from[pos]+to[pos]-1)*(from[pos]+to[pos]-5)==0)  Q[i*n+j]*=kappa;
      if(c[0]!=c[1])  Q[i*n+j]*=omega;
      Q[j*n+i]=Q[i*n+j];
   }
   FOR(i,n) FOR(j,n) Q[i*n+j]*=pi[j];
   for (i=0,mr=0;i<n;i++) { 
      Q[i*n+i]=-sum(Q+i*n,n); mr-=pi[i]*Q[i*n+i]; 
   }
   FOR(i,n*n) Q[i]*=t/mr;

   xtoy(Q,P,n*n);  FOR(i,n) P[i*n+i]++;   // I + Qt 
   xtoy(Q,Q1,n*n);
   for (i=2,div=2;i<nterms;i++,div*=i) {  // k is divisor 
      matby(Q, Q1, Q2, n, n, n);
      for(j=0,mr=0;j<n*n;j++) { P[j]+=Q2[j]/div; mr+=fabs(Q2[j]); }
      mr/=div;
      // if(debug) printf("Pmat term %d: norm=%.6e\n", i, mr); 
      if (mr<e) break;
      xtoy(Q2,Q1,n*n);
   }

   FOR(i,n*n) if(P[i]<0) P[i]=0;

}

*/
