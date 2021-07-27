% this script is for single peptide test run, with plot 
fname='m27c_Liver.csv';
T=readtable(fname);
%%
i=27;  % pick an item
pep=T{i,13}; %peptide sequence location column 13
pep=pep{1};

pep=strrep(pep ,'*' ,'L' );  % replace * with L
out=peptide_mid(pep);  
MID_sim=[out.pct];   %simulated MID 

MID_measure=T{i,16:26};
sz=min(length(MID_sim),length(MID_measure));
MID_sim=MID_sim(1:sz);
MID_measure=MID_measure(1:sz);

MID_measure=MID_measure/sum(MID_measure)*100;
[MID_corr,fval]=pepcorr(pep,MID_measure)

corr=conv(MID_corr,MID_sim);
corr=corr(1:length(MID_sim));

figure
bar([MID_measure;corr;-MID_sim]');
%%
figure,
stem(MID_sim,'.')
hold on
stem(-MID_measure,'.')
hold on
stem(conv(MID_corr,MID_sim),'.g')

MID_sim(7)/MID_sim(1)
MID_measure(7)/MID_measure(1)