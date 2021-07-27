%if a combination contains items from the same element, e.g., "17O 18O2"
%the factor needs to be revised e.g.,num_O=6, nck(6,1)*nck(6,2) --> nck(6,1)*nck(5,2)

function tb=decouple(atoms,tb)
for i=1:length(tb)
  item=tb(i);
  tp=strfind(item.str,'O');
  if length(tp)==2
    n1=str2num(item.str(tp(1)+1));
    n2=str2num(item.str(tp(2)+1));
    nS=atoms(4);
    factor=nck(nS,n1)*nck(nS-n1,n2)/(nck(nS,n1)*nck(nS,n2));
    item.ab=item.ab*factor;    
  end

  tp=strfind(item.str,'S');
  if length(tp)==2
    n1=str2num(item.str(tp(1)+1));
    n2=str2num(item.str(tp(2)+1));
    nS=atoms(5);
    factor=nck(nS,n1)*nck(nS-n1,n2)/(nck(nS,n1)*nck(nS,n2));
    item.ab=item.ab*factor;  
  elseif length(tp)==3
    n1=str2num(item.str(tp(1)+1));
    n2=str2num(item.str(tp(2)+1));
    n3=str2num(item.str(tp(3)+1));
    nS=atoms(5);
    factor=nck(nS,n1)*nck(nS-n1,n2)*nck(nS-n1-n2,n3)/(nck(nS,n1)*nck(nS,n2)*nck(nS,n3));
    item.ab=item.ab*factor;      
  end
  tb(i)=item;
end




 