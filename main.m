% this script is for single peptide test run, with plot 
fname='m27c_Liver.csv';
T=readtable(fname);

%%
for i=1:size(T,1);  % loop over peptides
    i
pep=T{i,13}; %peptide sequence location column 13
pep=pep{1};

MID_measure=T{i,16:26};
MID_measure=MID_measure/sum(MID_measure); %normalize to 0-1
[MID_sim,MID_corr,fval]=pepcorr(pep,MID_measure);

MID_simcorr=conv(MID_corr,MID_sim);
MID_simcorr=MID_simcorr(1:length(MID_measure));
out(i).corr=MID_corr;
out(i).err=fval;
out2(i).sim=MID_sim;
out2(i).simcorr=MID_simcorr;
end

writetable(struct2table(out),'output.csv')

%%
% figure
% bar([MID_measure(:),MID_simcorr(:)]);
