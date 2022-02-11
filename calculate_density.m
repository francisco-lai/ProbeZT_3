function density_mat = calculate_density(Gr1, Ga1, MU, D, de, block_sequence, Gamma1, Gamma2, GammaL, GammaR, nn, lenN, EE, TT)
%CALCULATE_DENSITY Summary of this function goes here
%   Detailed explanation goes here
Ln = lenN;
Ga = Ga1;
Gr = Gr1;
ee = EE;
bb = 1/TT;
if nn ~= Ln
    density_mat = 0;
elseif nn == Ln
    
    lenE = 2*D/de;
    
    for zz = 1:size(Ga,2)
        for aa = 1:1:size(Ga,2)
            if block_sequence(1,zz) == 1
                gamma(aa,aa,aa) = Gamma1;
            elseif block_sequence(1,zz) == 0
                gamma(aa,aa,aa) = Gamma2;
            end
            if aa==1
                gamma(aa,aa,aa) = gamma(aa,aa,aa)+GammaL;
            end
            if aa==size(Ga,2)
                gamma(aa,aa,aa) = gamma(aa,aa,aa)+GammaR;
            end
        end
    end
    
    
    for aa=1:size(Ga,2)
        for ii = 1:1:lenE
            fermi_func = 1/(exp(bb*(ee(ii)-MU(aa)))+1);
            mult = Gr(:,:,ii)*gamma(:,:,aa)*Ga(:,:,ii)*fermi_func;
%             mult(mult<1e-15) = 0;
            tfda(:,:,ii) = mult.*de;
        end
        dens(:,:,aa) = sum(tfda,3);
    end
    density_mat = abs(real((1/2*pi)*sum(dens,3)));
    
end

