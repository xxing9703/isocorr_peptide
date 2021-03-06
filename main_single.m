%This function takes a single filename as input, and output the corrected
function main_single(fname)
%fname='m27c_Liver_test.csv';
[~,fn]=fileparts(fname);
T=readtable(fname);
for i=1:size(T,1)  % loop over peptides
   fprintf(['... <',fn,'>  i=',num2str(i),'\n']);
pep=T{i,13}; %peptide sequence location column 13
pep=pep{1};
mask = startsWith(T.Properties.VariableNames, 'Area');
MID_measure=T{i,mask};
poolsize=sum(MID_measure);
MID_measure=MID_measure/sum(MID_measure); %normalize to 0-1
try
 [MID_sim,MID_corr,fval]=pepcorr(pep,MID_measure);
 MID_simcorr=conv(MID_corr,MID_sim);
 MID_simcorr=MID_simcorr(1:length(MID_measure));
 out(i).sum=poolsize;
 out(i).corrM_0=MID_corr(1);
 out(i).corrM=MID_corr(2:end);
 out(i).err=fval;
 out2(i).sim=MID_sim;
 out2(i).simcorr=MID_simcorr;
 catch
   out(i).sum=nan;  
 end
end
Tout=struct2table(out);
Tout=[T(:,~mask),Tout];
writetable(Tout,[fn,'_corrected.csv']);

%
% figure
% bar([MID_measure(:),MID_simcorr(:)]);

