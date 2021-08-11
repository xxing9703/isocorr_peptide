%pep2mass returns monotopic m/z and atoms
function [mz,atoms]=pep2mass(pep)
H=1.00728;
aa='ARNDCEQGHOILKMFPUSTWYVm'; %m stands for oxidized M
aa_formula={'C3H7NO2';'C6H14N4O2';'C4H8N2O3';'C4H7NO4';'C3H7NO2S';'C5H9NO4';'C5H10N2O3';'C2H5NO2';'C6H9N3O2';'C5H9NO3';'C6H13NO2';'C6H13NO2';'C6H14N2O2';'C5H11NO2S';'C9H11NO2';'C5H9NO2';'C5H7NO3';'C3H7NO3';'C4H9NO3';'C11H12N2O2';'C9H11NO3';'C5H11NO2';'C5H11NO3S'};
aa_formula{5}='C9H14N2O4S';  % revision: cystine protection on (this is requested by user)
for i=1:length(pep)
    [a(i),~,c(:,i)]=formula2mass(aa_formula{aa==pep(i)});
end
mz=sum(a)-formula2mass('H2O')*(length(pep)-1)+H;
atoms=sum(c,2)';
atoms(3)=atoms(3)-2*(length(pep)-1)+1;
atoms(4)=atoms(4)-(length(pep)-1);



