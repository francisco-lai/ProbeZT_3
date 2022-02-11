load seq10\seq10_lowE.mat

Ga = Ga1;
Gr = Gr1;

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
        mult = Gr(:,:,ii)*gamma(:,:,aa)*Ga(:,:,ii);
        mult(mult<1e-15) = 0;
        tfda(:,:,ii) = mult*de;
    end
    dens(:,:,aa) = sum(tfda,3);
end

full_density = real((1/2*pi)*sum(dens,3));

