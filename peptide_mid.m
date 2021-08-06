function [out,tb]=peptide_mid(peptide,maxM,ext,type)
%##input 
% maxM is the maximum M to calculate up to.  default is 8(if not specified)
% ext and type are optional inputs, the combination of these two inputs are for labeled MID simulation. 
% type =1,2 or 3 representing (C, N, or H) 
% ext is the corresponding MID (measured from labeled)specified by type 

%##output
% out is a structure containing the m/z and abundance of M+0, M+1,....
% tb is a structure containing the detailed info for each individual monoisotopic
% peaks (absolute abundance, name in str such as '13C2 15N 18O', m/z, dm/z, type:a flag used internally; pct: normalized to each M+i)  

%##isotopes considered in this simulation: 13C; 15N; 2H; 17O,18O; 33S,34S,36S

if nargin==1
    maxM=8;
end
%----------------------------------------------------
[mass,atoms]=pep2mass(peptide); %ordering: C,N,H,O,S
str={'13C','15N','2H','17O','33S','','','','18O','34S','','','','','36S'};
n=5;
A0=[0.98893,0.996337,0.99985,0.9975904,0.9502]; %abundance for M0
A1=[0.0110694,0.003663,0.00015,3.7387E-04,0.0074962]; %abundance type A for +1
B1=[0,0,0,0.0020358,0.042099]; %abundance type B for +2
C1=[0,0,0,0,0.00020529]; %abundance type C for +4

mzA=[1.00335,0.99703,1.00630,1.00422,0.99940]; %mz diff for +1
mzB=[0,0,0,2.00424,1.99580]; %mz diff for +2
mzC=[0,0,0,0,3.9950]; %mz diff for +4

%------------------------------------------

ab_0=prod(A0.^atoms(1:n)); %abundance for M0
% calculate abundance up to m
m=maxM;
for i=1:n
   for j=1:m
    ab_A(i,j)=ab_0*(A1(i)/A0(i))^j*nck(atoms(i),j);
    ab_B(i,j)=ab_0*(B1(i)/A0(i))^j*nck(atoms(i),j);
    ab_C(i,j)=ab_0*(C1(i)/A0(i))^j*nck(atoms(i),j);
    dmz_A(i,j)=mzA(i)*j;
    dmz_B(i,j)=mzB(i)*j;  
    dmz_C(i,j)=mzC(i)*j; 
   end
 end
M=[ab_A;ab_B;ab_C]; %MID of singly labeled(e.g, "33S" is single, but "13C 33S" is not single)
ct=0; % 
%-------------------- revise the MID of singly labeled 
if nargin==4
    %M1=M; %for debug
    ab_0_new=ab_0/prod(A0(type)^atoms(type))*(ext(1)/sum(ext)); %revise ab_0
    for k=1:size(M,1)
     for i=1:size(M,2)    
          M(k,i)=M(k,i)/ab_0*ab_0_new;  %revise M table by rescaling
     end
    end
    for i=1:length(ext)-1 
        M(type,i)=ext(i+1)/ext(1)*ab_0_new; %modify type
    end
    ab_0=ab_0_new;  
    %M2=M; %for debug
end
% ----------- write into tb
N=[dmz_A;dmz_B;dmz_C];
for i=1:size(M,1)
  for j=1:size(M,2)
      if M(i,j)>0
       ct=ct+1;
       tb(ct).ab=M(i,j);
       tb(ct).str=[str{i},num2str(j)];
       tb(ct).dmz=N(i,j); 
       tb(ct).mz=N(i,j)+mass; 
       tb(ct).type=i;
      end
  end
end

[~,ind]=sort([tb.ab],'descend'); %sorting
tb=tb(ind);
cutoff=0.00001;%apply cutoff
ind=find([tb.ab]>cutoff);  
tb=tb(ind);
%tb_save1=tb; % for debug
%% ----tb append to consider all combinations
ct=length(tb);
for i=1:length(tb)
    for j=i:length(tb)
        if isempty(intersect(tb(i).type,tb(j).type))
            tp=tb(i).ab*tb(j).ab/ab_0;
          if tp>cutoff
              ct=ct+1;
              tb(ct).ab=tp;
               str_in=[tb(i).str,' ',tb(j).str];
              tb(ct).str=strjoin( sort(strsplit(str_in,' ')));
              tb(ct).dmz=tb(i).dmz+tb(j).dmz; 
              tb(ct).mz=tb(i).dmz+tb(j).dmz+mass; 
              tb(ct).type=union(tb(i).type,tb(j).type);
          end              
        end
    end
end
 [~,ind]=unique({tb.str});  %remove duplicates
 tb=tb(ind);

 tb=decouple(atoms,tb); %decouple, fix coeffients of exclusive items (e.g, "17O 18O")
 
 [~,ind]=find([tb.ab]>cutoff); %apply cutoff, remove ~zero items
 tb=tb(ind);
%% ---------- sorting option1 (not used)
 [~,ind]=sort([tb.ab],'descend'); %sorting by abudance
tb=tb(ind);
%tb_save2=tb; %for debug
%% -------------sorting option2 (used)
 [~,ind]=sort([tb.dmz],'ascend'); %sorting by dmz
tb=tb(ind);
%tb_save3=tb; %for debug
%% -------------- sum up for M+i and normalize
out(1).IsotopeNumber=0;
out(1).mz=mass;
out(1).pct=ab_0*100;
out(1).pctMax=100;
for i=1:round(max([tb.dmz]))
  out(i+1).IsotopeNumber=i;
  [~,ind]=find(round([tb.dmz])==i);
  out(i+1).mz=sum([tb(ind).mz].*[tb(ind).ab])/sum([tb(ind).ab]);
  out(i+1).pct=sum([tb(ind).ab])*100;
  for j=1:length(ind)
    tb(ind(j)).pct=tb(ind(j)).ab/max([tb(ind).ab])*100;
  end  
end
for i=1:round(max([tb.dmz]))+1
   out(i).pctMax=out(i).pct/max([out.pct])*100;
end