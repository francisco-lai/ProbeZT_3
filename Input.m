
calctype=1; % Model (1=linear-chain   2=DNA) 
probe=2; %Probe toggle (1=dephasing  2=voltage  3=Voltage-Temperature) 

%% PARAMETERS
%simulation variables
if calctype==1 %linear model
    N=[3:3:12]; % number of sites in the molecule; dimension of the system Hamiltonian [loop variable]
    tn1=0.2; % inter-site tunneling element (eV), off-diagonal entries in the Hamiltonian 
    tn0=0.05;
    epn=0.5; %site energies of the localized molecular orbitals (eV) set relative to the Fermi energy (0). 
    % for accuracy, one should set epn off-resonance (at least kBT off the Fermi energy)
    sequence_cell=NaN;
    block_pattern = [1 1 1 0 0 0]; % 1 corresponds to Gamma1, 0 corresponds to Gamma0
    tn_pattern = [1 1 1 0 0 0]; % 1 corresponts to tn1, 0 corresponds to tn0
end

if calctype==2 % double-stranded DNA model
    sequence_cell = {'acgcgcgt','acgcagcgt', 'acgcatgcgt',  'acgcatagcgt', 'acgcatatgcgt', 'acgcatatagcgt', 'acgcatatatgcgt', 'acgcatatatagcgt', 'acgcatatatatgcgt' }; %enter the 5' to 3' single-strand sequence
    %the matching antiparallel strand will be built in automatically
    epn=NaN; tn=NaN;
end

Temperature = [ 300 ]; %temperature of the metal leads (K) [loop variable]  

Gamma1 =[ 0.05 ]; %hybridization energy between system and probes (eV) [loop variable]
Gamma0 = 0.09; %hybridization energy between system and probes (eV) NOT a loop variable

% we assume that incoherent-probe effects are local, uncorrelated, and uniform

GammaL=[ 0.05 ]; % hybridization of the first site in the chain to the left electrode (eV)  [loop variable]
GammaR=[ 0.05 ]; % hybridization of the last site in the chain to the right electrode (eV) [loop variable]
                     

%integration parameters
D = 5; % bandwidth; the metals extend from -D to +D (eV)
de = 1e-3; % energy step in numerical integration; smallest energy parameter (eV). 
% These parameters critically determine run time 
