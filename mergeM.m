%input: M is a matrix, each row is measured isotopemer (13C or 15N or 2H) of one amino acid.
%conv() is used to combine different amino acids into one final isotopemer distribution. 
function output=mergeM(M)
M=M./sum(M,2);  %normalize
tp=M(1,:);
if size(M,1)>1
 for i=2:size(M,1)
   tp=conv(tp,M(i,:));
 end
end
output=tp(1:size(M,2));