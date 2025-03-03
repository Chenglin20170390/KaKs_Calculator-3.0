PRGS = KaKs KnKs AXTConvertor ConPairs
CC = g++ #gcc
CFLAGS = -O -w

all:	$(PRGS) 

KaKs: KaKs_main.cpp KaKs.cpp MSMA.cpp MYN.cpp base.cpp  NG86.cpp LWL85.cpp LPB93.cpp GY94.cpp YN00.cpp KaKs.h MSMA.h MYN.h base.h NG86.h LWL85.h LPB93.h GY94.h YN00.h
	$(CC) $(CFLAGS) -o  $@  KaKs_main.cpp KaKs.cpp MSMA.cpp MYN.cpp base.cpp NG86.cpp LWL85.cpp LPB93.cpp GY94.cpp YN00.cpp -lstdc++ -lm

KnKs: KnKs_main.cpp KnKs.cpp base.cpp GY94.cpp KnKs.h base.h GY94.h
	$(CC) $(CFLAGS) -o  $@  KnKs_main.cpp KnKs.cpp base.cpp GY94.cpp -lstdc++ -lm  

AXTConvertor: AXTConvertor.cpp
	$(CC) $(CFLAGS) -o  $@  AXTConvertor.cpp -lstdc++ -lm

ConPairs: ConcatenatePairs.cpp
	$(CC) $(CFLAGS) -o  $@  ConcatenatePairs.cpp -lstdc++ -lm

