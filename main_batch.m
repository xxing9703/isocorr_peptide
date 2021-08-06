[flist,path]=uigetfile('.csv','selected files','MultiSelect','on');
if ~iscell(flist)
    cc{1}=flist;
    flist=cc;
end
for i=1:length(flist)
    main_single(fullfile(path,flist{i}))
end
