% plot simulated results
offset=32; %starting peptide = 1+offset
load output
fname='m27c_Liver.csv';
T=readtable(fname);
for i=1:16
    k=i+offset;
    subplot(4,4,i)
    MID_measure=T{k,16:26};
MID_measure=MID_measure/sum(MID_measure);
MID_simcorr=out2(k).simcorr;
    bar([MID_measure(:),MID_simcorr(:)])
end