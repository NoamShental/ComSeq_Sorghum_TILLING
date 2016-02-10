function [samples_heterozygous,samples_homozygous]=findLine(y,M0)
% Input:
% measurement vector y
% normalized measuremnt matrix M0

% Output:
% samples_heterozygous - the samples that hold a heterozygous SNP
% samples_homozygous - the samples that hold a homozygous SNP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%==========================================
% perform reconstruction
tau = 0.005*max(abs(M0'*y)); % regularization coefficient

% apply GPSR
fractionalOutput = applyGPSR(y,M0,tau);

% since the output of GPSR is fractional we "round" it to the 
% most likely integer. 
% We look for the solution with minimal squre error with the measurements
discreteOutput = modifySolution(fractionalOutput,M0,y);

% define the SNP type according to the output


%==========================================
% perform reconstruction
tau = 0.005*max(abs(M0'*y));

% fractional solution

fractionalOutput = applyGPSR(y,M0,tau);

% look for the solution with minimal squre error with the measurements
discreteOutput = modifySolution(fractionalOutput,M0,y);
samples_heterozygous = find(discreteOutput==1);
samples_homozygous = find(discreteOutput==2);