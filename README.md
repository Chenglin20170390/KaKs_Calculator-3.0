# KaKs_Calculator-3.0)![github](https://img.shields.io/badge/3C-Certification-red)        
![github](https://img.shields.io/badge/Coding--sequence-hex)        ![github](https://img.shields.io/badge/Non--coding--sequence-green)         ![github](https://img.shields.io/badge/Dn/Ds-red)       


- Fork of the KaKs_Calculator-3.0 project https://ngdc.cncb.ac.cn/biocode/tools/BT000001.
- Calculating selective pressure on coding and non-coding sequences

## Building KaKs_Calculator-3.0
```
git clone https://github.com/Chenglin20170390/KaKs_Calculator-3.0.git
cd KaKs_Calculator-3.0/bin;make
```

## Run KaKs_Calculator-3.0
- for coding aligned sequences
```
./KaKs_Calculator-3.0/bin/KaKs -h
```

- for non-coding aligned sequences
```
./KaKs_Calculator-3.0/bin/KnKs -h
```

## Example
```
./bin/KaKs -i examples/coding.axt -o test.kaks
```

## Tutorial from the fasta file 
###  input cds and pep sequence 
- mafft alignment 
```
sample=ID
mafft --auto  test.pep.fa   >  $sample.pep.mafft.fa
```
- check alignment from pep and convert alignment to axt format (https://github.com/liaochenlanruo/PAL2NAL/blob/master/pal2nal.pl)
```
pal2nal.pl -nogap -nomismatch $sample.pep.mafft.fa  cds.list.$sample.fa  -output clustal > cds.$sample.aln.clustal.fa
AXTConvertor cds.$sample.aln.clustal.fa cds.$sample.aln.clustal.aux
```
- Ka/Ks analysis
```
KaKs -i cds.$sample.aln.clustal.aux -o cds.$sample.aln.clustal.aux.kaka.txt

##format to R
awk '{print $1,$5}'  test.cds.mafft.aux.kaks.txt | sed 's/&/-/g' |sed '1d'|sed '1iSample\tGene' > KaKs.txt
```
## R code (example of repeat genes)
```
library(ggplot2)
library(pheatmap)

dat <-read.table('KaKs.txt',header=T)  ## read input file
dat$Gene[dat$Gene > 5] <- 0   ## Filtering KaKs that large than 5

rownames(dat) <- dat$Sample  ##rename row
dat<-dat[,rep(2,8)]  ##repeat 8 times for the interest gene (example)
head(dat)  ## check data format

### heatmap figure
p<-pheatmap(dat, color = colorRampPalette(c( "white","blue"))(10),
            legend_labels = c("low", "median", "high"),
            border=T,display_numbers = F,
            clustering_distance_rows = "euclidean",
            cluster_cols=T,
            cluster_rows=T,
            clustering_method="ward.D")
p
ggsave(p, filename= "Kaks.heatmap.pdf" , width=30, height=30,  units =c( "cm" ),colormodel= "srgb" )  ##save pdf figure
```
![Kaks](https://github.com/user-attachments/assets/c8de14e6-a5d8-4182-bcde-91333598b66f)






# References

Zhang Z. KaKs_Calculator 3.0: calculating selective pressure on coding and non-coding sequences[J]. Genomics, Proteomics and Bioinformatics, 2022, 20(3): 536-540.
https://doi.org/10.1016/j.gpb.2021.12.002
