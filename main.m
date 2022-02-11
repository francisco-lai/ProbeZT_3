clear
tic % check runtime

%Input
file_in = input('What file would you like to open? ', 's');
run(file_in);

definitions

%Output name
filename = strcat('seq10/test2');

% make the block sequence for gamma
block_sequence = make_seq(block_pattern, N);

% make the block sequence for tn
tn_sequence = make_seq(tn_pattern, N);


for nn=1:1:lenN %loop over chains of different length (linear, calctype=1) or (DNA sequences, calctype=2)
    % Build the hamiltonian matrix hh
    if calctype==1
        hh=hamiltonianLinear(N(nn),epn,tn0,tn1,tn_sequence);
    elseif calctype==2
        hh=hamiltonianDNA(sequence_cell{nn});
    end
    
    for igp=1:1:lenG %loop over probe coupling-strength
        
        for iglr=1:lenGLR %loop over hybridization to electrodes
            [TRpp, TRlr, TRpl, TRpr, Gr1, Ga1] = transmission_zero(hh, EE, GammaL(iglr), GammaR(iglr), Gamma1(igp), Gamma0(igp), block_sequence,nn);
            for tt=1:1:lenT %loop over temperature
                if probe==1 %dephasing probe
                    [GL(nn,igp,iglr,tt), GR(nn,igp,iglr,tt)] = ...
                        dprobe(N(nn), EE, TT(tt), muM, voltage, TRpp, TRpl, TRpr, TRlr);
                    
                elseif probe==2 %voltage probe
                    [GL(nn,igp,iglr,tt), GR(nn,igp,iglr,tt), MU] = ...
                        vprobe(N(nn), EE, TT(tt), muM, voltage, TRpp, TRpl, TRpr, TRlr);
                elseif probe == 3 %voltage-temperature probe
                    [GL(nn,igp,iglr,tt), GR(nn,igp,iglr,tt), SL(nn,igp,iglr,tt), SR(nn,igp,iglr,tt), ...
                        KL(nn,igp,iglr,tt), KR(nn,igp,iglr,tt), ZTL(nn,igp,iglr,tt), ZTR(nn,igp,iglr,tt)] ...
                        = vtprobe(N(nn), EE, TT(tt), muM, voltage, TRpp, TRpl, TRpr, TRlr);
                end
            end
            clear TRpl TRpr TRlr TRpp
        end
    end
    density_mat = calculate_density(Gr1, Ga1, MU, D, de, block_sequence, Gamma1, Gamma0, GammaL, GammaR, nn, lenN, EE, TT);
    if probe == 1
        save (filename, 'GL', 'GR', 'IL', 'IR', 'probe', 'D', 'de', 'TT', 'N', 'Gamma1', 'Gamma0', 'GammaL', 'GammaR',...
            'epn', 'tn1','tn0', 'sequence_cell', 'block_sequence','tn_sequence','density_mat')
    elseif probe ==2
        save (filename, 'GL', 'GR', 'IL', 'IR', 'probe', 'D', 'de', 'TT', 'N', 'Gamma1', 'Gamma0', 'GammaL', 'GammaR',...
            'epn', 'tn1','tn0', 'sequence_cell', 'block_sequence','tn_sequence', 'voltage','density_mat')
    elseif probe == 3
        save (filename,'dataVT', 'GL', 'GR', 'IL', 'IR', 'SL', 'SR', 'KL', 'KR', 'ZTL', 'ZTR', 'probe', 'D', 'de', 'TT', 'N', 'Gamma1', 'Gamma0',...
            'GammaL', 'GammaR', 'epn', 'tn1', 'tn0', 'sequence_cell', 'block_sequence','tn_sequence','density_mat')
    end
end
toc
