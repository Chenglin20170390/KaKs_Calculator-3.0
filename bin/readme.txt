KaKs_Calculator, Version 3.0

	1. The source codes of several methods for calculating Ka, Ks and Kn. Each method is defined as a C++ class. 
		(0) base(.h/.cpp): a base class, including several widely used functions
		(1) NG: NG86.h NG86.cpp
		(2) LWL, MLWL: LWL85.h LWL85.cpp
		(3) LPB, MLPB: LPB93.h LPB93.cpp
		(4) GY: GY94.h GY94.cpp
		(5) YN: YN00.h YN00.cpp
		(6) MYN: MYN.h MYN.cpp
		(7) MSMA: MSMA.h MSMA.cpp
		(8) ZZ: KnKs.h KnKs.cpp (Kn/Ks for estimating selection on non-coding sequences)
		
	2. AXTConverter.cpp
	AXTConverter is a program for converting Clustal/Msf/Nexus/Phylip/Pir format sequences to AXT ones.
	
	3. ConcatenatePairs.cpp
	Concatenate all pairs of sequences to a pair of sequences.
	
For compiling, just type 'make' command.
	
