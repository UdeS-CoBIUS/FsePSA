# FsePSA

FsePSA compute pairwise alignments between coding sequences while accounting for the length of frameshift translations.


### Requirements

FsePSA is composed of a set of Python scripts. It requires the following to be available:

* python (2.7)
* python packages: numpy and biopython


### Usage
```
src/fse_main.py [-h] [-go GAPOPEN] [-ge GAPEXTEND] [-fso FSOPEN]
                [-fse FSEXTEND] [-aa AMINOACIDMATRIX]
                [-d DATASEQUENCE DATAFORMAT] [-o OUTFILE]
               
```


  *-h*, --help   show this help message and exit
  
  *-go , --gapopen \<GAPOPEN>*         gap openning cost   

  *-ge , --gapextend \<GAPEXTEND>*         gap extension cost   

  *-fso , --fsopen \<FSOPEN>*         frameshift openning cost   

  *-fso , --fsextend \<FSEXTEND>*         frameshift extension cost   

  *-aa , --aminoacidmatrix \<AMINOACIDMATRIX>*        amino acid substitution score matrix file  

  *-d , --datasequence \<DATASEQUENCE fasta> *      sequence file in fasta format    

  *-o , --outfile \<OUTFILE> *      output file   


#### Input files
##### Amino acid substitution score matrix.

Example:
```
#  Matrix made by matblas from blosum62.iij
#  * column uses minimum score
#  BLOSUM Clustered Scoring Matrix in 1/2 Bit Units
#  Blocks Database = /data/blocks_5.0/blocks.dat
#  Cluster Percentage: >= 62
#  Entropy =   0.6979, Expected =  -0.5209
A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  Z  X  *
A  4 -1 -2 -2  0 -1 -1  0 -2 -1 -1 -1 -1 -2 -1  1  0 -3 -2  0 -2 -1  0 -4 
R -1  5  0 -2 -3  1  0 -2  0 -3 -2  2 -1 -3 -2 -1 -1 -3 -2 -3 -1  0 -1 -4 
N -2  0  6  1 -3  0  0  0  1 -3 -3  0 -2 -3 -2  1  0 -4 -2 -3  3  0 -1 -4 
D -2 -2  1  6 -3  0  2 -1 -1 -3 -4 -1 -3 -3 -1  0 -1 -4 -3 -3  4  1 -1 -4 
C  0 -3 -3 -3  9 -3 -4 -3 -3 -1 -1 -3 -1 -2 -3 -1 -1 -2 -2 -1 -3 -3 -2 -4 
Q -1  1  0  0 -3  5  2 -2  0 -3 -2  1  0 -3 -1  0 -1 -2 -1 -2  0  3 -1 -4 
E -1  0  0  2 -4  2  5 -2  0 -3 -3  1 -2 -3 -1  0 -1 -3 -2 -2  1  4 -1 -4 
G  0 -2  0 -1 -3 -2 -2  6 -2 -4 -4 -2 -3 -3 -2  0 -2 -2 -3 -3 -1 -2 -1 -4 
H -2  0  1 -1 -3  0  0 -2  8 -3 -3 -1 -2 -1 -2 -1 -2 -2  2 -3  0  0 -1 -4 
I -1 -3 -3 -3 -1 -3 -3 -4 -3  4  2 -3  1  0 -3 -2 -1 -3 -1  3 -3 -3 -1 -4 
L -1 -2 -3 -4 -1 -2 -3 -4 -3  2  4 -2  2  0 -3 -2 -1 -2 -1  1 -4 -3 -1 -4 
K -1  2  0 -1 -3  1  1 -2 -1 -3 -2  5 -1 -3 -1  0 -1 -3 -2 -2  0  1 -1 -4 
M -1 -1 -2 -3 -1  0 -2 -3 -2  1  2 -1  5  0 -2 -1 -1 -1 -1  1 -3 -1 -1 -4 
F -2 -3 -3 -3 -2 -3 -3 -3 -1  0  0 -3  0  6 -4 -2 -2  1  3 -1 -3 -3 -1 -4 
P -1 -2 -2 -1 -3 -1 -1 -2 -2 -3 -3 -1 -2 -4  7 -1 -1 -4 -3 -2 -2 -1 -2 -4 
S  1 -1  1  0 -1  0  0  0 -1 -2 -2  0 -1 -2 -1  4  1 -3 -2 -2  0  0  0 -4 
T  0 -1  0 -1 -1 -1 -1 -2 -2 -1 -1 -1 -1 -2 -1  1  5 -2 -2  0 -1 -1  0 -4 
W -3 -3 -4 -4 -2 -2 -3 -2 -2 -3 -2 -3 -1  1 -4 -3 -2 11  2 -3 -4 -3 -2 -4 
Y -2 -2 -2 -3 -2 -1 -2 -3  2 -1 -1 -2 -1  3 -3 -2 -2  2  7 -1 -3 -2 -1 -4 
V  0 -3 -3 -3 -1 -2 -2 -3 -3  3  1 -2  1 -1 -2 -2  0 -3 -1  4 -3 -2 -1 -4 
B -2 -1  3  4 -3  0  1 -1  0 -3 -4  0 -3 -3 -2  0 -1 -4 -3 -3  4  1 -1 -4 
Z -1  0  0  1 -3  3  4 -2  0 -3 -3  1 -1 -3 -1  0 -1 -3 -2 -2  1  4 -1 -4 
X  0 -1 -1 -1 -2 -1 -1 -1 -1 -1 -1 -1 -1 -1 -2  0  0 -2 -1 -1 -1 -1 -1 -4 
* -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4  1 
```

##### Sequence file in fasta format.

Example:
```
>Seq1
ATGACCGAATCCAAGCAGCCCTGGCATAAGTGGGGGAACGATTGA
>Seq2
ATGACCGAATCCAAGCAGCCCTGGCATAATGGGGGAACGATTGAAGTAGGAACGATTTAA
>Seq3
ATGACCGAATCCAACAGCCCTGGCATAAGTGGGGGAACGATTGAAGTAGGAACGATTTAA
```

#### Output file

* *outfile* --- contains all pairwise alignments at the srspair format. The starting positions 
of frameshift translations in the alignments are indicated with the symbol *!* in the markup line.

Example:
```
########################################
# Program: fse
# Rundate: Mon Feb 22 15:00:09 2016
# Commandline: python fse_main.py
#    --datasequence examples/example_data.fasta fasta
#    --gapopen -2.0
#    --gapextend -1.0
#    --fsopen -2.0
#    --fsextend -1.0
#    --aminoacidmatrix resources/BLOSUM62.txt
#    --outfile examples/example_data.srspair
# Align_format: srspair
# Report_file: examples/example_data.srspair
########################################

#=======================================
#
# Aligned_sequences: 2
# 1: Seq1
# 2: Seq2
# Matrix: resources/BLOSUM62.txt
# GapOpen_penalty: -2.0
# GapExtend_penalty: -1.0
# FSopen_penalty: -2.0
# FSextend_penalty: -1.0
#
# Length: 61
# Identity:       43/61 70.5%
# Similarity:     43/61 70.5%
# Gaps:           17/61 27.9%
# FSopen:         1/61 1.6%
# FrameShift:     12/61 19.7%
# Score: 62.5
#
#
#=======================================

Seq1              1 ATGACCGAATCCAAGCAGCCCTGGCATAAGTGGGGGAACGAT--------     42
                    |||||||||||||||||||||||||||||!||||||||||||        
Seq2              1 ATGACCGAATCCAAGCAGCCCTGGCATAA-TGGGGGAACGATTGAAGTAG     49

Seq1             43 --------TGA     45
                            |.|
Seq2             50 GAACGATTTAA     60

#---------------------------------------
#---------------------------------------

#=======================================
#
# Aligned_sequences: 2
# 1: Seq1
# 2: Seq3
# Matrix: resources/BLOSUM62.txt
# GapOpen_penalty: -2.0
# GapExtend_penalty: -1.0
# FSopen_penalty: -2.0
# FSextend_penalty: -1.0
#
# Length: 61
# Identity:       43/61 70.5%
# Similarity:     43/61 70.5%
# Gaps:           17/61 27.9%
# FSopen:         1/61 1.6%
# FrameShift:     27/61 44.3%
# Score: 47.0
#
#
#=======================================

Seq1              1 ATGACCGAATCCAAGCAGCCCTGGCATAAGTGGGGGAACGAT--------     42
                    ||||||||||||||!|||||||||||||||||||||||||||        
Seq3              1 ATGACCGAATCCAA-CAGCCCTGGCATAAGTGGGGGAACGATTGAAGTAG     49

Seq1             43 --------TGA     45
                            |.|
Seq3             50 GAACGATTTAA     60

#---------------------------------------
#---------------------------------------

#=======================================
#
# Aligned_sequences: 2
# 1: Seq2
# 2: Seq3
# Matrix: resources/BLOSUM62.txt
# GapOpen_penalty: -2.0
# GapExtend_penalty: -1.0
# FSopen_penalty: -2.0
# FSextend_penalty: -1.0
#
# Length: 61
# Identity:       58/61 95.1%
# Similarity:     58/61 95.1%
# Gaps:           2/61 3.3%
# FSopen:         1/61 1.6%
# FrameShift:     12/61 19.7%
# Score: 80.5
#
#
#=======================================

Seq2              1 ATGACCGAATCCAAGCAGCCCTGGCAT-AATGGGGGAACGATTGAAGTAG     49
                    ||||||||||||||!|||||||||||| |.||||||||||||||||||||
Seq3              1 ATGACCGAATCCAA-CAGCCCTGGCATAAGTGGGGGAACGATTGAAGTAG     49

Seq2             50 GAACGATTTAA     60
                    |||||||||||
Seq3             50 GAACGATTTAA     60

#---------------------------------------
#---------------------------------------

#=======================================
#
# Score_matrix
#
#=======================================

Seq1	Seq2	Seq3	
62.5	47.0	
80.5	

#---------------------------------------
#---------------------------------------
```

#### Running FsePSA on an example

To run FsePSA on the example of sequence file *example_data.fasta* given in directory *examples*,
and using the *BLOSUM62* amino acid substitution matrix given in directory *resources*
you can use the following command:

```
python src/fse_main.py --datasequence examples/example_data.fasta fasta \
--gapopen -2.0 --gapextend -1.0 --fsopen -2.0 --fsextend -1.0 \
--aminoacidmatrix resources/BLOSUM62.txt \
--outfile examples/example_data.srspair
```




