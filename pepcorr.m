%input peptide sequence, measured MID, returns simulated MID (natural abundance)
%and natural ab corrected MID_corr, fval returns the err 
%constrain is optional, set the upper bounds (0-1) for corrected MID
function [MID_sim,MID_corr,fval]=pepcorr(pep,MID_measure,constrain)
pep=strrep(pep ,'*' ,'L' ); % replace * with L in sequence

out=peptide_mid(pep);
MID_sim=[out.pct];  % get simulated MID for natural abundance(0-100)
MID_sim=MID_sim/sum(MID_sim); %(0-1)
% sz=min(length(MID_sim),length(MID_measure));  % trancate to the shorter one sz
% MID_sim=MID_sim(1:sz);
% MID_measure=MID_measure(1:sz);

% optimization setup to find MID_corr
sz=length(MID_measure);
x0=zeros(sz,1)/sz;x0(1)=1;  %initial x
lb=zeros(sz,1);
ub=ones(sz,1);
if nargin==3
    ub=constrain;
end
if sum(MID_measure)>0
MID_measure=MID_measure/sum(MID_measure);  % normalize to (0-1)
[xopt,fval]=fmincon(@(x)myfit(x,MID_sim,MID_measure),x0,[],[],[],[],lb,ub);
MID_corr=xopt'/sum(xopt);
else
MID_corr=zeros(1,length(MID_measure));
fval=nan;
   
end

end

function fit=myfit(x,MID_sim,MID_measure)
y=conv(x,MID_sim);  %convolution
if length(y)>length(MID_measure)
  y=y(1:length(MID_measure)); %truncate
else
  y(length(MID_measure))=0;   %padding 0
end
y=y/sum(y);
fit=sqrt(sum((y(:)-MID_measure(:)).^2));
end
