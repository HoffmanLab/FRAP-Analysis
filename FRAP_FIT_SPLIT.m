function FRAP_FIT_SPLIT(dname,npblch,keywords)

if ~isfield(keywords,'dt')
    keywords.dt=1;
end
if ~isfield(keywords,'ctrain')
    keywords.ctrain=[0,0,0];
end

split_files = file_search(['split_FRAP_blch_' dname '\w+.dat'],keywords.folder);

for k = 1:length(split_files)
    data = load(split_files{k});
    [r,~] = size(data);
    for i = 1:r
        res = data(i,:);
        res=res(npblch:end);
        nele=length(res);
        t=(1:nele).*keywords.dt;
        if ~isfield(keywords,'sfo')
            fo=min(res);
        else
            fo=keywords.sfo;
        end
        if ~isfield(keywords,'sff')
            ff=max(res);
        else
            ff=keywords.sff;
        end
        w=find(res > fo & res < ff);
        if ~isfield(keywords,'sk')
            k=-log((ff-res(w))./(ff-fo))./t(w);
            k=mean(k);
        else
            k=keywords.sk;
        end
        a=[fo,ff,k];
        
        if isfield(keywords,'tmax')
            w=find(t < keywords.tmax);
            t=t(w);
            res=res(w);
        end
        
        ft = fittype('FRAP_FIT_FUNC(x,a,b,c)');
        f = fit(t',res',ft,'StartPoint',a);
        disp(f)
        
        t12=log(2)/f.c;
        
        if keywords.showfit
            figure; hold on;
            plot(f,t,res);
            legend(['TC = ' num2str(t12)])
            title(strrep(dname,'_',' '))
            axis([0 max(t) 0 1])
        end
    end
end

end