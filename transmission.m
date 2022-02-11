function [TRpp, TRlr, TRpl, TRpr, Ga1, Gr1]= transmission (hamiltonian,energy_grid,gamma_l,gamma_r,gamma_p1,gamma_p2, block_sequence, nn)
%TRANSMISSION builds Green's function and computes the transmission

Ln=length(hamiltonian); II=eye(Ln);
lenE=length(energy_grid);
gamma=zeros(Ln); gamma(1,1)=1i/2*gamma_l; 
gamma(Ln,Ln)=gamma(Ln,Ln)+1i/2*gamma_r;

for qq=1:1:Ln
    if block_sequence(qq) == 1
        gamma_p(qq,qq)= gamma_p1;
    elseif block_sequence(qq) == 0
        gamma_p(qq,qq)= gamma_p2;
    end
end

gamma=gamma + 1i/2*gamma_p;
% disp(gamma_p)
%building the gamma matrix: coupling to the terminals and probes
%disp(gamma_p)
HH=hamiltonian-gamma; clear gamma
Gamma_P= zeros(Ln,Ln,Ln); Gamma_L= zeros(Ln); Gamma_R=zeros(Ln);
Gamma_L(1,1)=gamma_l; Gamma_R(Ln,Ln)=gamma_r;

% hybridization matrices describing coupling to the terminals and each of the probes
if gamma_p(nn,nn)~=0
    for aa=1:1:Ln
        if block_sequence(aa) == 1
            Gamma_P(aa,aa,aa)=gamma_p1;
        elseif block_sequence(aa) == 0
            Gamma_P(aa,aa,aa)=gamma_p2;
        end
    end
    disp(Gamma_P(aa,aa,aa))
    for ii=1:1:lenE
        Gr1(:,:,ii)=inv(energy_grid(ii)*II-HH); Ga1(:,:,ii)=permute(conj(Gr1(:,:,ii)),[2 1 3]); %Green's function  matrices
        Ga = Ga1(:,:,ii);
        Gr = Gr1(:,:,ii);
%         Gr(ii)=inv(energy_grid(ii)*II-HH); Ga=Gr'; %Green's function  matrices
        for aa=1:1:Ln
            for aap=1:1:Ln
%                 if aap>=aa
                    TRpp(ii,aa,aap)=trace(Gamma_P(:,:,aa)*Ga*Gamma_P(:,:,aap)*Gr);  % transmission from probe p to probe p'
%                 else
%                     TRpp(ii,aa,aap)=TRpp(ii,aap,aa);
%                 end
            end
            
            TRpl(ii,aa)=trace(Gamma_L*Ga*Gamma_P(:,:,aa)*Gr);   % transmission from probe p to the left terminal l
            TRpr(ii,aa)=trace(Gamma_R*Ga*Gamma_P(:,:,aa)*Gr);   % transmission from probe p to the right terminal r
        end
        TRlr(ii)=trace(Gamma_L*Ga*Gamma_R*Gr);  % transmission from the left terminal to the right directly (not through probes)
    end
    
else %if no probes, transmission is calculated at each energy from left to right directly
    for ii=1:1:lenE
        Gr=inv(energy_grid(ii)*II-HH); Ga=Gr';
        TRlr(ii)=trace(Gamma_L*Ga*Gamma_R*Gr);
        TRpp=NaN;TRpl=NaN;TRpr=NaN;
    end
end

