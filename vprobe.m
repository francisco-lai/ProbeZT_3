function [condl, condr, MU] = vprobe(length_Hamiltonian, energy_grid, temp_eV, fermi_energy, voltage, TRpp, TRpl, TRpr, TRlr)
%VPROBE voltage probe. It calculates linear response conductance.
% If the probe coupling is zero, coherent current is calculated. 
% If it is non-zero, the matrix MM is computed.  
% A linear equation provides the  probes' local potentials.
%The calculation of the current follows. 
% Output is the conductance of the junction.

Ln=length_Hamiltonian;
EE=energy_grid; de=energy_grid(2)-energy_grid(1); 
b=1/temp_eV; muL=fermi_energy+voltage/2; muR=fermi_energy-voltage/2;

if isnan(TRpp) % if probes' coupling is zero 
    fl=1./(exp(b.*(EE-muL))+1); %Fermi Dirac functions at the L and R leads
    fr=1./(exp(b.*(EE-muR))+1);
    Il=sum(TRlr.*(fl-fr).*de); Ir=-Il; %current when coupling to probes is zero
    MU = 0
     
else % if probes are on
    Fpl=zeros(1,Ln); Fpr=zeros(1,Ln); Fpp=zeros(Ln); MM=zeros(Ln); 
    dfdE=-b./(4*cosh(b*EE/2).^2); %derivative of the fermi function with respect to energy
    Flr=-sum(TRlr.*dfdE)*de;  % elements of matrix MM
    
    for aa=1:Ln
        Fpl(aa)=-sum(TRpl(:,aa).*dfdE')*de;
        Fpr(aa)=-sum(TRpr(:,aa).*dfdE')*de;

        for aap=1:1:Ln
            Fpp(aa,aap)=-sum(TRpp(:,aa,aap).*dfdE')*de;
            MM(aa,aap)=-Fpp(aa,aap);
        end
    end
    
    for aa=1:1:Ln
        MM(aa,aa)=0;
        for aap=1:1:Ln
            if aap~=aa
                MM(aa,aa)=MM(aa,aa)+Fpp(aap,aa);
            end
        end
        MM(aa,aa)=MM(aa,aa)+Fpl(aa)+Fpr(aa);
    end
    clear Fpp dfdE
    
    vec=zeros(1,Ln);
    for aa=1:Ln
        vec(aa)=muL*Fpl(aa)+muR*Fpr(aa); % inhomogeneous term for the linear equation.
    end
    
    if (rcond(MM) < 1e-12); MU=pinv(MM)*(vec');
    else
        MU=MM\(vec'); % stores probes' chemical potential
    end
    
    Ilr=Flr*(muL-muR);
    Ilp=sum(Fpl(1:Ln)*(muL-MU(1:Ln)));
    Il=Ilr+Ilp;

    % electric current due to applied voltage is computed from left to right terminal and vice versa
    Irp=sum(Fpr(1:Ln)*(muR-MU(1:Ln)));
    Ir=-Ilr+Irp;

end

%disp(MU)

condl = Il/(voltage); condr = Ir/(voltage);

end
