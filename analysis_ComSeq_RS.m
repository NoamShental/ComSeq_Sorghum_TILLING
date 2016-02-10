clear

%%%%%%%%%%%%%%%%
% a)	Load the experimental data and the Reed-Solomon (RS) measurement matrix, M. 

load experimentData 
% this file contains a single struct named "data" that contains experimetnal read
% counts for each amplicon at each position along the R1 and R2 read. 

% The file contains a cell array named "data" of size 4x1. 
% Each cell is a struct in itself and corresponds to one 
% of the 4 amplicons tested.
% The fields in this struct and their description appear below:

% name: Amplicon name
% forward: Sequence of forward primer
% reverse: Sequence of reverse primer
% wholeAmpliconForward: The amplicon sequence, forward
% wholeAmpliconReverse: The amplicon sequence, reverse

% r1: Is a cell of size 100x1 where 100 is the read length used. 
% "r1" is the R1 part of the paired-end read. 
% Each of the cell's entries corresponds to a position along the read, 
% and holds a 48x4 matrix. 
% Each column of this matrix corresponds to A/C/G/T nucleotides, respectively, 
% and each row is a pool. 
% The entry in this matrix corresponds to the number of reads that showed 
% a specific nucleotide. 
% For example, data{1}.r1{20}(10,3) is the number of reads that showed "G" 
% in pool #10 at position 20 of the r1 amplicon of CST18.

% r2: The same for the R2 part of the paired end read.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load RSdesign % the Reed-Solomom (RS) measurement matrix

nOrder = ['ACGT']; % nucleotide order in data file. The first column is A, 
% the second is C, etc.

readLength = 100; % The experimental read length was 100.
% end - part a)
%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%
% b) Filter locations along the read that potentially hold a mutation: 
% This is independently performed for R1 and R2 for each amplicon, 

% find potential loci along R1 
R1_flag = 1; % flag for R1 (1) or R2 (0). Needed since the first bp in R1 are barcode etc.
for i=1:length(data) % for each of the 4 amplicons 
  potentialPositions_R1{i}=filter_locations(data{i}.r1,RScode_design,readLength,R1_flag);
end

% find potential loci along R2
R1_flag = 0;
for i=1:length(data) % for each of the 4 amplicons 
  potentialPositions_R2{i}=filter_locations(data{i}.r2,RScode_design,readLength,R1_flag);
end

% end - part b)
%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%
% c) Detect the carriers using ComSeq
% For each of the potential locations in (b) the measurement 
% vector y is calculated and using the matrix M the carrier is detected.
% This is performed for R1 and R2

% normalize the measuremnt matrix by twice the number of samples per pool
% see "Detecting rare alleles and their carriers via ComSeq" in Methods.
M0 = RScode_design/256; % 256 is twice the number of alleles in a pool

% find lines for each potential position - R1
for i=1:length(data)
  for j=1:size(potentialPositions_R1{i},1)
    position = potentialPositions_R1{i}(j,1);
    major_allele_count = data{i}.r1{position}(:,potentialPositions_R1{i}(j,3));
    minor_allele_count = data{i}.r1{position}(:,potentialPositions_R1{i}(j,2));
    y = minor_allele_count./(major_allele_count+minor_allele_count);
    
    % perform compressed sensing
    [samples_heterozygous,samples_homozygous]=findLine(y,M0);
    
    % display results
    if ~isempty(samples_homozygous)
      disp(['Amplicon ',data{i}.name,'; Homozygous SNP at location: ',num2str(position),' in R1; major allele: ',nOrder(potentialPositions_R1{i}(j,3)),', minor allele: ',nOrder(potentialPositions_R1{i}(j,2))])
    end
    if ~isempty(samples_heterozygous)
      disp(['Amplicon ',data{i}.name,'; Heterozygous SNP at location: ',num2str(position),' in R1; major allele: ',nOrder(potentialPositions_R1{i}(j,3)),', minor allele: ',nOrder(potentialPositions_R1{i}(j,2))])
    end
  end
end

% find lines for each potential position - R2
for i=1:length(data)
  for j=1:size(potentialPositions_R2{i},1)
    position = potentialPositions_R2{i}(j,1);
    major_allele_count = data{i}.r2{position}(:,potentialPositions_R2{i}(j,3));
    minor_allele_count = data{i}.r2{position}(:,potentialPositions_R2{i}(j,2));
    y = minor_allele_count./(major_allele_count+minor_allele_count);

    % perform compressed sensing
    [samples_heterozygous,samples_homozygous]=findLine(y,M0);
    
    % display results
    if ~isempty(samples_homozygous)
      disp(['Amplicon ',data{i}.name,'; Homozygous SNP at location: ',num2str(position),' in R2; major allele: ',nOrder(potentialPositions_R2{i}(j,3)),', minor allele: ',nOrder(potentialPositions_R2{i}(j,2))])
    end
    if ~isempty(samples_heterozygous)
      disp(['Amplicon ',data{i}.name,'; Heterozygous SNP at location: ',num2str(position),' in R2; major allele: ',nOrder(potentialPositions_R2{i}(j,3)),', minor allele: ',nOrder(potentialPositions_R2{i}(j,2))])
    end
  end
end
