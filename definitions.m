%% Physical constants
ee=1.60217662e-19; % charge of electron in C, conversion factor from Joule (J) to eV
kb=1.3806485e-23; %Boltzmann constant in J/K  
KB=kb/ee; %Boltzmann constant in eV/K

%% DNA-specific inputs
if calctype==2 %DNA model
    N=zeros(1,length(sequence_cell)); % vector, stores the lengths of the different sequences
    for i=1:length(sequence_cell) %counting the number of sites for DNA calculation
        N(i) = 2*length(cell2mat(sequence_cell(i)));
    end
end

%% Calculation parameters
muM = 0; % Fermi energy, set here at 0.
voltage=5e-3; % Voltage bias across the junction (eV); In linear response calculations this value does not influence results.

TT=Temperature*KB; %Temperature*KB; % temperature in eV (k_B*T)

EE=-D:de:D; % discretized metallic bands

lenN=length(N); lenE=length(EE); lenG=length(Gamma1); lenGLR=length(GammaL); lenT=length(TT);

% creating four dimensional matrices for output
% (sequence_length), (coupling_to_probes) , (coupling_to_terminals), (temperature)
outMat=[lenN lenG lenGLR lenT];
GR=zeros(outMat); GL=zeros(outMat);
IL=zeros(outMat); IR=zeros(outMat);

if probe==3 % more output matrices for VT-probe
    SR=zeros(outMat); SL=zeros(outMat);
    KR=zeros(outMat); KL=zeros(outMat);
    ZTR=zeros(outMat); ZTL=zeros(outMat);
end
