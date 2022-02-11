function [condl, condr] = dprobe(length_Hamiltonian, energy_grid, temp_eV, fermi_energy, voltage, TRpp, TRpl, TRpr, TRlr)
% dprobe, dephasing probe. It calculates conductance for the dephasing probe.
% If the probe coupling is zero, coherent current is calculated. 
% If it is non-zero, the matrix MM is computed. A linear equation provides the 
% probes' distribution functions. The calculation of the current follows. 
% The output is the conductance of the junction.

% Definitions
b=1/temp_eV; Ln=length_Hamiltonian;
EE=energy_grid; lenE=length(energy_grid); de=energy_grid(2)-energy_grid(1);
muL = fermi_energy+voltage/2; muR = fermi_energy-voltage/2;
fl=1./(exp(b.*(EE-muL))+1); fr=1./(exp(b.*(EE-muR))+1);
% Fermi Dirac functions for the L and R electrodes

if isnan(TRpp); % if probes are off
    Il=sum(TRlr.*(fl-fr).*de); Ir=-Il; %current in the case of no probes
else % if probes are on
    
    MM=zeros(lenE,Ln,Ln); % create the matrix MM based on the transmission functions
    for aa=1:1:Ln
        for aap=1:1:Ln
            MM(:,aa,aap)=-TRpp(:,aa,aap);
        end
        MM(1:lenE,aa,aa)=0;
        for aap=1:1:Ln
            if aap~=aa
                MM(:,aa,aa)=TRpp(:,aa,aap) + MM(:,aa,aa);
            end
        end
        MM(:,aa,aa)=MM(:,aa,aa)+TRpl(:,aa)+TRpr(:,aa);
    end
    clear TRpp
    
    % the probes' functions are computed based on the probe condition 
    vec=zeros(lenE,Ln); fp=zeros(lenE,Ln);
    Ilp=zeros(1,Ln);    Irp=zeros(1,Ln);% create empty arrays
    %build inhomogeneous vector
    for aa=1:1:Ln
        vec(:,aa)=TRpl(:,aa).*fl' + TRpr(:,aa).*fr';
    end
    %calculate unknown probes' functions
    for ii=1:1:lenE
        Mn(1:Ln,1:Ln)=MM(ii,1:Ln,1:Ln);
        if (rcond(Mn) < 1e-12)
            fp(ii,1:Ln)=pinv(Mn)*(vec(ii,1:Ln)');
        else
            fp(ii,1:Ln)=Mn\(vec(ii,1:Ln)');
        end
    end
    Ilr=sum(TRlr.*(fl-fr))*de;
    for aa=1:1:Ln
        Ilp(aa)= sum(TRpl(:,aa).*(fl'-fp(:,aa)))*de;
        Irp(aa)= sum(TRpr(:,aa).*(fr'-fp(:,aa)))*de;
    end
    Il=Ilr+sum(Ilp);        Ir=-Ilr+sum(Irp);
end
condl = Il/voltage; condr = Ir/voltage;
disp(Il)
end
