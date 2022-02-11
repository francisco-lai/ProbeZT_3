%create Hamiltonian hh for linear model
% uniform chain with energies given by site_energy 
% and site-to-site nearest neighbor couplings are assumed.

function hh = hamiltonianLinear(number_of_sites, site_energy, tn0, tn1, tn_seq)
Ln=number_of_sites; hh=zeros(Ln); epn=site_energy;
    for jj=1:1:Ln-1
        if tn_seq(jj) == 0
            tn = tn0;
        elseif tn_seq(jj) == 1
            tn = tn1;
        end
        hh(jj,jj)=epn; hh(jj,jj+1)=tn; hh(jj+1,jj)=tn;
    end
    hh(end,end)=epn;
end
