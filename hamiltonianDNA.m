%hamiltonianDNA creates Hamiltonian hh for a ds-DNA of a specified sequence using a tight binding ladder model.
% On-site energies and coupling parameters are defined here. 
% Because of nearest-neighbor interactions, three nucleotides are
% examined at a time. The nucleotides are coded as A=1, C=2, G=3, T=4
%Parameters are taken from: 
%[1] K. Senthilkumar,  F. C. Grozema, C. F. Guerra, F. M. Bickelhaupt, F. D. Lewis,  Y. A. Berlin,  M. A. Ratner, and L. D. A. Siebbeles, J. Am. Chem. Soc. 127, 14894 (2005).
% Energy parameters are averaged, as done in 
%[2] M. Zilly, O Ujsaghy, and D. E. Wolf, Phys. Rev. B 82, 125125 (2010).
%The Fermi energy is set at the energy of the G base

function hh=hamiltonianDNA (sequence)
seq_length=2*length(sequence);
seq = zeros(1,seq_length/2); %note: 'seq' is listed  3' to 5' 
%(unlike 'sequence', the input, which is listed 5' to 3')
for i=1:seq_length/2
    if sequence(i)=='a'||sequence(i)=='A'; seq(i)=4;
    elseif sequence(i)=='c'||sequence(i)=='C'; seq(i)=3;
    elseif sequence(i)=='g'||sequence(i)=='G'; seq(i)=2;
    elseif sequence(i)=='t'||sequence(i)=='T'; seq(i)=1;
    end
end

% Energy matrix and interactions 
energy = zeros(4,4,4); % energy of each site can be generalized to depend on neighboring sites as in Ref. [1]. 
energedge = [8.631 9.722 8.178 9.464]; %defining the energies on the edges (no neigbours)
intv = [-0.047 -0.055 -0.055 -0.047]; %vertical interactions-  between strands
int53 = [-0.004 0.042 -0.01 -0.063; -0.002 0.022 0.009 -0.055; -0.077 -0.114 0.053 0.141; -0.031 -0.028 0.018 0.072];
int35 =int53'; %matched by position
int55 = [0.031 -0.001 -0.013 0.007; -0.001 0.001 0.002 0.0003; -0.013 0.002 0.012 -0.009; 0.007 0.0003 -0.009 0.001];
int33 = [0.049 0.017 -0.011 -0.007; 0.017 0.01 0.022 0.004; -0.011 0.022 -0.032 -0.014; -0.007 0.004 -0.014 0.006];
energy(:,1,:)=8.631; energy(:,2,:)=9.722; energy(:,3,:)=8.178; energy(:,4,:)=9.464;
EFermi=energedge(3); %Setting Fermi energy to the energy of the G base
energynew = energy-EFermi;  % setting energies relative to the Fermi energy
energedgenew = energedge-EFermi;
%% Filling in hh
hh=zeros(seq_length);
for i=2:(seq_length/2-1)
    a=[seq(i-1) seq(i) seq(i+1)]; %the triple containing the nucleotide of interest and its two nearest neighbors
    for ii = 1:3 % adding another triple of complementary nucleotides, now a contains 6 numbers
        if a(ii)==1; a(ii+3)=4;
        elseif a(ii)==2; a(ii+3)=3;
        elseif a(ii)==3; a(ii+3)=2;
        elseif a(ii)==4; a(ii+3)=1;
        end
    end
    %filling in the hamiltonian based on the interactions with the
    %nearest neighbors (according to dobble ladder model)
    hh(2*i-1:2*i,2*i-1:2*i) = [energynew(a(1),a(2),a(3)) intv(a(2)); intv(a(2)) energynew(a(4),a(5),a(6))];
    hh(2*i-3:2*i-2,2*i-1:2*i)= [int35(a(1),a(2)) int33(a(1),a(5)); int55(a(2),a(4)) int53(a(4),a(5))];
    hh(2*i-1:2*i,2*i-3:2*i-2)= [int35(a(1),a(2)) int33(a(1),a(5)); int55(a(2),a(4)) int53(a(4),a(5))]';
end
begend = [seq(1) seq(seq_length/2)]; %first and last nucleotide in the sequence
for i = 1:2 % adding two complementary nucleotides
    begend(i+2)=1;
    if begend(i)==1; begend(i+2)=4;
    elseif begend(i)==2; begend(i+2)=3;
    elseif begend(i)==3; begend(i+2)=2;
    end
end
%filling in the first and last pairs of nucleotides into hamiltonian
hh(1:2,1:2)=[energedgenew(begend(1)) intv(begend(1)); intv(begend(1)) energedgenew(begend(3))];
hh(seq_length-1:seq_length,seq_length-1:seq_length)=[energedgenew(begend(2)) intv(begend(2)); intv(begend(2)) energedgenew(begend(4))];
hh(seq_length-3:seq_length-2,seq_length-1:seq_length)= [int35(a(2),a(3)) int33(a(2),a(6)); int55(a(3),a(5)) int53(a(5),a(6))];
hh(seq_length-1:seq_length,seq_length-3:seq_length-2)= [int35(a(2),a(3)) int33(a(2),a(6)); int55(a(3),a(5)) int53(a(5),a(6))]';
end

