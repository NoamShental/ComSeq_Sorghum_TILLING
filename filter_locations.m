function [potentialSNPs]=filter_locations(curr_data,RScode_design,readLength,R1_flag)
% Find locations along the read that potentially hold a mutation 

% Input
% curr_data is the data of a single amplicon - either R1 or R2. 
% RScode_design is the measurement (pooling) matrix
% readLength is the read length
% R1_flag=1 for r1 data and 0 for R2.

% Output
% potentialSNPs = Each potential SNP is defined by its 
% [position,minor allele nucleotide index, major allele nucleotide index]
% The order of nucleotides is: ACGT, namely potentialSNPs=[10 1 3] means
% a potential SNP in location 10 where the major allele is G, 
% and minor allele is A.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% define the first and last positions to scan along the read. 
% These differ as r1 contains the barcode and the primer in this case
if R1_flag
  first_position = 8+18+1;
  last_position = readLength;
else
  first_position = 1;
  last_position = readLength;
end

% detect loci
% scan the relevant part of the read for positions in which the fraction
% of the minor allele is high enough
potentialSNPs = [];
for l=first_position:last_position
  
  % find the major and minor allele for the current location
  b = sum(curr_data{l});
  % ord(1) is the major allele, ord(2) is the minor allele
  [junk,ord] = sort(b,'descend'); 
  
  % find the fraction of the minor allele counts
  ratio = [curr_data{l}(:,ord(2))./(  curr_data{l}(:,ord(1))+curr_data{l}(:,ord(2))  )];
  
  % a location is set as a "potential SNP" if its fraction of minor allele
  % count is larger than 1/256/2 and there enough reads for this location
  % In addition - it has to be consistent with the RS code:
  % "a" is the set of pools whose minor allele fraction is high enough.
  % These "a" pools are checked using the RS code for consistency.
  % For a loci to be called, there should be at least one sample that
  % shares between 4 and 6 pools with the "a" pools detected.
  a = find(ratio>1/256/2 & curr_data{l}(:,ord(1))>10^3 & curr_data{l}(:,ord(2))>50);
  d = find(sum(RScode_design(a,:),1)>=4 & length(a)<=6);
 
  % Each potential SNP is defined by its [position,minor allele, major allele]
  if ~isempty(d) && median(ratio(a))>2*median(ratio)
      potentialSNPs = [potentialSNPs;l,ord(2),ord(1)];
  end

 
end








