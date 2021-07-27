% modified function of nchoosek, accepting n<k which returns 0
function out=nck(n,k)
if n<k 
    out=0;
else
    out=nchoosek(n,k);
end