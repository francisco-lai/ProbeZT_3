function [TRpp, TRlr, TRpl, TRpr, Ga1, Gr1]= transmission_zero(hamiltonian,energy_grid,gamma_l,gamma_r,gamma_p1,gamma_p2, block_sequence, nn)
%TRANSMISSION builds Green's function and computes the transmission

Ln=length(hamiltonian); II=eye(Ln);
lenE=length(energy_grid);
gamma=zeros(Ln); gamma(1,1)=1i/2*gamma_l; 
gamma(Ln,Ln)=gamma(Ln,Ln)+1i/2*gamma_r;

TRpp = zeros(lenE, Ln, Ln);
TRpl = sparse(lenE, Ln);
TRpr = sparse(lenE, Ln);
TRlr = sparse(lenE);

for qq=1:1:Ln
    if block_sequence(qq) == 1
        gamma_p(qq,qq)= gamma_p1;
    elseif block_sequence(qq) == 0
        gamma_p(qq,qq)= gamma_p2;
    end
end

gamma=gamma + 1i/2*gamma_p;
%disp(gamma_p)
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
    %disp(Gamma_P(aa,aa,aa))
    parfor ii=1:1:lenE
        Gr=(energy_grid(ii)*II-HH)\eye(Ln); Ga=Gr'; %Green's function  matrices
        GLGa = Gamma_L*Ga;
        GRGa = Gamma_R*Ga;
        for aa=1:1:Ln
            GmGa=Gamma_P(:,:,aa)*Ga;
            GmGr2=Gamma_P(:,:,aa)*Gr;
            for aap=1:1:Ln
                GmGr=Gamma_P(:,:,aap)*Gr;
                TRpp(ii,aa,aap)=GmGa(:).'*reshape(GmGr.',[],1);
            end
            TRpl(ii,aa)= GLGa(:).'*reshape(GmGr2.',[],1);   % transmission from probe p to the left terminal l
            TRpr(ii,aa)= GRGa(:).'*reshape(GmGr2.',[],1);   % transmission from probe p to the right terminal r
        end
        TRlr(ii)=trace(Gamma_L*Ga*Gamma_R*Gr);  % transmission from the left terminal to the right directly (not through probes)
    end

    for ii=1:1:lenE
        Gr1(:,:,ii)=inv(energy_grid(ii)*II-HH); Ga1(:,:,ii)=Gr1(:,:,ii)';
    end


else %if no probes, transmission is calculated at each energy from left to right directly
    for ii=1:1:lenE
        Gr=inv(energy_grid(ii)*II-HH); Ga=Gr';
        TRlr(ii)=trace(Gamma_L*Ga*Gamma_R*Gr);
        TRpp=NaN;TRpl=NaN;TRpr=NaN;
    end
end

