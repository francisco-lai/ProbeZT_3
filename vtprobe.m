function [condl, condr, sl, sr, kl, kr, ztl, ztr] = vtprobe ...
    (length_Hamiltonian, energy_grid, temp_eV, fermi_energy, voltage, TRpp, TRpl, TRpr, TRlr)
%VTPROBE voltage temperature probe. It calculates conductance, thermopower, electronic thermal conductance, and ZT figure of merit 
% If the probe coupling is zero, coherent current is calculated. 
% If it is non-zero, the matrix MM is computed. 
% A linear equation provides the probes' local potentials and temperatures.
% Calculation of the charge and heat currents follows.


% Definitions
Ln=length_Hamiltonian; EE=energy_grid; de=energy_grid(2)-energy_grid(1);
b=1/temp_eV; t=temp_eV; tl=999e-3*t; tr=1001e-3*t; bl=1/tl; br=1/tr;  % (tr-tl)/t=0.002 
muM=fermi_energy; muL=fermi_energy+voltage/2; muR=fermi_energy-voltage/2;

if isnan(TRpp); % if coupling to probes is zero  
    fl=1./(exp(b.*(EE-muL))+1); fr=1./(exp(b.*(EE-muR))+1);
    Il=sum(TRlr.*(fl-fr).*de); Ir=-Il; %charge current  due to voltage bias 
    El=sum(TRlr.*(fl-fr)*de.*EE); Er=-El; %energy current due to voltage bias 

    fl=1./(exp(bl.*(EE-muM))+1); fr=1./(exp(br.*(EE-muM))+1);
    Il_t=sum(TRlr.*(fl-fr).*de); Ir_t=-Il_t; %charge current due to temperature difference
    El_t=sum(TRlr.*(fl-fr)*de.*EE); Er_t=-El_t; %energy current due to temperature difference

else % if probes are on
    Fpl=zeros(1,Ln); FEpl=zeros(1,Ln); FTpl=zeros(1,Ln); FETpl=zeros(1,Ln);
    Fpr=zeros(1,Ln); FEpr=zeros(1,Ln); FTpr=zeros(1,Ln); FETpr=zeros(1,Ln);
    Fpp=zeros(Ln); FEpp=zeros(Ln); FTpp=zeros(Ln); FETpp=zeros(Ln);
    %derivatives of Fermi function with respect to energy and temperaure. 
    %The negative sign of dfdT is missing here, but it is added in the calculation of the 'Fp' and 'FE' functions below.
    dfdE=-b./(4.*cosh(b.*EE./2).^2); dfdT=dfdE.*(EE).*b;
    
    %% Computing the integrals
    Flr=-sum(TRlr.*dfdE)*de;
    FTlr = -sum(TRlr.*dfdT)*de; FETlr = -sum(TRlr.*dfdT.*EE)*de;
    for aa=1:1:Ln
        Fpl(aa)=-sum(TRpl(:,aa).*dfdE')*de;     FEpl(aa)=-sum(TRpl(:,aa).*dfdE'.*EE')*de;
        FTpl(aa)=-sum(TRpl(:,aa).*dfdT')*de;    FETpl(aa)=-sum(TRpl(:,aa).*dfdT'.*EE')*de;
        Fpr(aa)=-sum(TRpr(:,aa).*dfdE')*de;     FEpr(aa)=-sum(TRpr(:,aa).*dfdE'.*EE')*de;
        FTpr(aa)=-sum(TRpr(:,aa).*dfdT')*de;    FETpr(aa)=-sum(TRpr(:,aa).*dfdT'.*EE')*de;
        for aap=1:1:Ln;
            Fpp(aa,aap)=-sum(TRpp(:,aa,aap).*dfdE')*de;        FEpp(aa,aap)=-sum(TRpp(:,aa,aap).*dfdE'.*EE')*de;
            FTpp(aa,aap)=-sum(TRpp(:,aa,aap).*dfdT')*de;       FETpp(aa,aap)=-sum(TRpp(:,aa,aap).*dfdT'.*EE')*de;
        end
    end
    
    %% Costruction of the MM-matrix by 4 blocks
    MM = zeros (2*Ln); %dimension is  (2 Ln x 2 Ln),
    % first Ln rows: contributions to electric current, last Ln rows - contributions to thermal energy current.
    % first Ln columns: contributions from voltage bias, last Ln columns - contributions from temperature difference.
    for aa=1:1:Ln
        for aap=1:1:Ln
            MM(aa,aap)=-Fpp(aa,aap);            MM(aa,aap+Ln) = -FTpp(aa,aap);
            MM(aa+Ln,aap) = -FEpp(aa,aap);      MM(aa+Ln,aap+Ln) = -FETpp(aa,aap);
        end
        MM(aa,aa)=0;        MM(aa,aa+Ln)=0;
        MM(aa+Ln,aa)=0;     MM(aa+Ln, aa+Ln) = 0;
        
        for aap=1:1:Ln
            if aap~=aa
                MM(aa,aa)=MM(aa,aa)+Fpp(aap,aa);            MM(aa,aa+Ln)=MM(aa,aa+Ln)+FTpp(aap,aa);
                MM(aa+Ln,aa)=MM(aa+Ln,aa)+FEpp(aap,aa);     MM(aa+Ln,aa+Ln)=MM(aa+Ln,aa+Ln)+FETpp(aap,aa);
            end
        end
        
        MM(aa,aa)=MM(aa,aa)+Fpl(aa)+Fpr(aa);                MM(aa,aa+Ln)=MM(aa,aa+Ln)+FTpl(aa)+FTpr(aa);
        MM(aa+Ln,aa)=MM(aa+Ln,aa)+FEpl(aa)+FEpr(aa);        MM(aa+Ln,aa+Ln)=MM(aa+Ln,aa+Ln)+FETpl(aa)+FETpr(aa);
    end
    
    %% (1) Voltage-driven current, constant temperature t
    vec=zeros(1,2*Ln);
    for aa=1:1:Ln
        vec(aa)=muL*Fpl(aa)+muR*Fpr(aa)+t*FTpl(aa)+t*FTpr(aa); % inhomogeneous term for the linear equation.
        vec(aa+Ln)=muL*FEpl(aa)+muR*FEpr(aa)+t*FETpl(aa)+t*FETpr(aa); % inhomogeneous term for the linear equation.
    end
    if (rcond(MM) < 1e-12); MU=pinv(MM)*(vec');
    else
        MU=MM\(vec'); % stores both mu's and T's for probes (dT=0, only voltage bias is applied)
    end

    Ilr=Flr*(muL-muR);
    Ilp=sum(Fpl(1:Ln)*(muL-MU(1:Ln))+FTpl(1:Ln)*(t-MU(1+Ln:2*Ln)));
    Irp=sum(Fpr(1:Ln)*(muR-MU(1:Ln))+FTpr(1:Ln)*(t-MU(1+Ln:2*Ln)));
    % electric current due to applied voltage is computed from left to right terminal and vice versa
    Il=Ilr+Ilp;     Ir=-Ilr+Irp;
     
    
    Elr=FETlr*(muL-muR);
    Elp=sum(FEpl(1:Ln)*(muL-MU(1:Ln))+FETpl(1:Ln)*(t-MU(1+Ln:2*Ln)));
    Erp=sum(FEpr(1:Ln)*(muR-MU(1:Ln))+FETpr(1:Ln)*(t-MU(1+Ln:2*Ln)));
    El=Elr+Elp;     Er=-Elr+Erp;

 
    %% (2) Currents due to temperature difference, constant chemical potential muM
    vec=zeros(1,2*Ln);
    
    for aa=1:1:Ln
        vec(aa) = muM*Fpl(aa) + muM*Fpr(aa) + tl*FTpl(aa) + tr*FTpr(aa); 
        vec(aa+Ln) = muM*FEpl(aa) + muM*FEpr(aa)+tl*FETpl(aa) + tr*FETpr(aa); 
    end
    if (rcond(MM) < 1e-12); MU=pinv(MM)*(vec');
    else
    MU=MM\(vec'); %stores both mu's and T's  (dV=0, no potential difference applied)
    end
    
    Ilr_t=FTlr*(tl-tr);
    Ilp_t=sum(Fpl(1:Ln)*(muM-MU(1:Ln))+FTpl(1:Ln)*(tl-MU(1+Ln:2*Ln)));
    Irp_t=sum(Fpr(1:Ln)*(muM-MU(1:Ln))+FTpr(1:Ln)*(tr-MU(1+Ln:2*Ln)));
    % current to both L and R terminals is calculated. They are supposed to be of the same magniture, but with opposite sign ---if the calculation is consistent.
    Il_t=Ilr_t+Ilp_t;    Ir_t=-Ilr_t+Irp_t;
    
    Elr_t=FETlr*(tl-tr);
    Elp_t=sum(FEpl(1:Ln)*(muM-MU(1:Ln))+FETpl(1:Ln)*(tl-MU(1+Ln:2*Ln)));
    Erp_t=sum(FEpr(1:Ln)*(muM-MU(1:Ln))+FETpr(1:Ln)*(tr-MU(1+Ln:2*Ln)));
    % similarly, energy current is calculated at both the left and right terminals
    El_t=Elr_t+Elp_t;     Er_t=-Elr_t+Erp_t;
end


%% output % 
condl = Il/(voltage);                   condr = Ir/(voltage);        % conductances in units of G0
sl=(Il_t./condl./(tl-tr));              sr=(Ir_t./condr./(tl-tr));    % thermopower in multiples of kb/e 
Pil=El/Il;                              Pir=Er/Ir;
kl=(El_t-Pil*Il_t)./(tl-tr) ;            kr=(Er_t-Pir*Ir_t)/(tl-tr); % electronic thermal conductance, units of energy (eV).  
ztl=abs(condl*sl^2*t/kl);               ztr=abs(condr*sr^2*t/kr); %dimensionless figure of merit ZT

end
