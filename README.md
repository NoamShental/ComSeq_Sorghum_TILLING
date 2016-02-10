**ComSeq feasibility study results**
================================

This package contains the basic Matlab scripts used to analyze the experimental data described in the manuscript: "Highly efficient de novo mutant identification in a Sorghum bicolor TILLING population using the ComSeq approach" by Nida et al.

**Experimental Data**: Preprocessed experimental read data of 1024 Sorghum bicolor lines sequenced via 48 pools is provided in the file "experimentData.mat". The file contains a cell array named "data" of size 4x1. Each cell is a struct in itself and corresponds to one of the 4 amplicons tested.
The fields in this struct and their description appear below:

*name*: Amplicon name

*forward*: Sequence of forward primer

*reverse*: Sequence of reverse primer

*wholeAmpliconForward*: The amplicon sequence, forward

*wholeAmpliconReverse*: The amplicon sequence, reverse

*r1*: Is a cell of size 100x1 where 100 is the read length used. "r1" is the R1 part of the paired-end read. Each of the cell's entries corresponds to a position along the read, and holds a 48x4 matrix. Each column of this matrix corresponds to A/C/G/T nucleotides, respectively, and each row is a pool. The entry in this matrix corresponds to the number of reads that showed a specific nucleotide. For example, data{1}.r1{20}(10,3) is the number of reads that showed "G" in pool #10 at position 20 of the r1 amplicon of CST18.

*r2*: The same for the R2 part of the paired end read.

For details as to how these data were processed from the raw reads, see Analysis section in Methods.


**Code**: The code allows detecting rare SNPs and their carriers based on the ComSeq approach. The input is the abovementioned data file. 
The file "analysis_ComSeq_RS.m" contains the following parts:

a) Load the experimental data and the Reed-Solomon (RS) measurement matrix, M. 

b) Filter locations along the read that potentially hold a mutation: This is independently performed for R1 and R2 for each amplicon, using the function "filter_locations.m" (further description is provided in the function itself).  

c) Detect the carriers using ComSeq: For each of the potential locations in (b) the measurement vector y is calculated, and using the matrix M the carrier is detected via the function "findLine". The function applied the GPSR solver to find either homozygous or heterozygous carriers (further description is provided in the function itself). 

Running the file "analysis_ComSeq_RS.m" performs steps a-c and displays the carriers found for in each amplicon.

**Compiled Matlab script**: The directory includes a stand alone version of the script analysis_ComSeq_RS, called "analysis_ComSeq_RS.exe". The code was compiled for Win7 64bit. Running the code is similar to running the Matlab code and outputs the list of SNPs found. To run this on a Win7 machine do the following:

a) Download the file MCR_R2013a_win64_installer.exe from http://www.mathworks.com/products/compiler/mcr/
(look for release R2013a (8.1) for Windows 64-bit. It is free to download and install).

b) Following downling the file double click on the "MCR_R2013a_win64_installer.exe" file. This would install the Matlab runtime component needed. The software would ask for a directory name. For example, assume it is C:\runTime 

c) Type the following at the command line: 
       set %PATH%;C:\runTime

d) Change the directory to where you downloaded the whole "ComSeq_Sorthum_TILLING" directory and type "analysis_ComSeq_RS.exe".
       

**Data in text format**: The R1 and R2 data also appear in text format. Each amplicon has a separate subdirectory that holds a text file for each position along the amplicon. For example, the file CST18/CST18_pos_20_r1.txt holds a 48x4 data matrix as described above for R1 of the CST18 amplicon at position 20. 

Information about each amplicom appears in "GeneralData.xls"

**Acknowledgements**:
These scripts use a specific file taken from the GPSR package (relevant file included). 
We thank the Mario Figueiredo for allowing us to include this file in our package.

Please send any comments/bug reports to: Noam Shental, shental@openu.ac.il
