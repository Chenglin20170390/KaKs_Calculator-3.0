/***************************************************************
* 2021, BIG of Chinese Academy of Sciences
* All rights reserved.
 
* Filename: KnKs_main.cpp
* Abstract: estimate seletive pressure on non-coding sequences.

* Version: 1.0
* Author: Zhang Zhang (zhanghzhang@genomics.org.cn)
* Date: Apr.2, 2006
****************************************************************/

#include "KnKs.h"

int main(int argc, const char* argv[]) {

	try {
		KNKS kk;
		if(!kk.Run(argc, argv)) throw 1;
	}
	catch (...) {
		
	}
	return 0;
}



