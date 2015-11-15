function FRAP_FIT(dname,npblch,keywords)

% f = frap_fit('avg.control','FRAP_FIT_FUNC',2,keywords)

% A function to read in a FRAP recovery curve and fit to an arbitrary
% function defined by the variable fname.
% Currently, the data is fit to the following function:
%   f = a(2)-(a(2)-a(1))*exp(-(a(3))*x)
% Guesses for the parameters are taken from the data
% a(1) is the intensity directly after bleach and is initialized to the
%       minimum of the data
% a(2) is the recovery intensity and is initiallized to the maximum of the
%       data
% a(3) is the recovery rate and is estimated from the data by taking the
%       mean of the following formula:
%       k = -a*log((ff(w)-res(w))/ff(w)-fo(w)))/t(w)

% INPUTS:
% dname - regular expression representing name of file that contains data
% npblch - the number of images take before the bleach
% keywords - a structure containing the following fields:
%   dt - contains the time between images
%   showfit - if set to 1, it will plot the data and the fit
%   ctrain - constrains the fit parameters to the guess, set a vector the
%       size of the number of fit parameters. 0 is optimized, 1 is constrained
%   sfo - contains the initial intensity parameter guess
%   sff - contains the final intensity parameter guess
%   tmax - contains the maximum time to be used

% OUTPUTS:
% a - contains the fit parameters (3 elements)

% 9/28/12 - File created by Katheryn Rotheberg

if ~isfield(keywords,'dt')
    keywords.dt=1;
end
if ~isfield(keywords,'ctrain')
    keywords.ctrain=[0,0,0];
end
% get files
% fn = file_search(dname,keywords.folder);
bkgf = file_search(['FRAP_bkg_' dname '.dat'],keywords.folder);
bkg = load(fullfile(keywords.folder,'FRAP Curve Files',bkgf{1}));
conf = file_search(['FRAP_con_' dname '.dat'],keywords.folder);
con = load(fullfile(keywords.folder,'FRAP Curve Files',conf{1}));
blchf = file_search(['FRAP_blch_' dname '.dat'],keywords.folder);
blch = load(fullfile(keywords.folder,'FRAP Curve Files',blchf{1}));
fn = file_search(['norm_FRAP_blch_' dname '.dat'],keywords.folder);

h1 = figure; hold on; set(h1,'Visible','off');
plot(bkg./max(max(bkg)),'g','LineWidth',3)
plot(con./max(max(con)),'r','LineWidth',3')
plot(blch./max(max(blch)),'b','LineWidth',3')
legend('Bkg','Con','Blch')
saveas(h1,fullfile(keywords.folder,'FRAP Curve Figures',['FRAPCurves_' dname '.png']))

data = load(fullfile(keywords.folder,'FRAP Curve Files',fn{1}));
[r,c] = size(data);
for i = 1:r
    res = data(i,:);
    res=res(npblch:end);
    nele=length(res);
    t=(1:nele).*keywords.dt;
%     t=t';
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
%     disp(f)
    coeffs = [f.a f.b f.c];
    save(fullfile(keywords.folder,'FRAP Curve Figures',['Coeffs_' dname '.txt']),'coeffs','-ascii')
    t12=log(2)/f.c;
    
    if keywords.showfit
        h2 = figure; hold on; set(h2,'Visible','off')
        plot(f,t,res);
        legend(['TC = ' num2str(t12)])
        title(strrep(dname,'_',' '))
        axis([0 max(t) 0 1])
        saveas(h2,fullfile(keywords.folder,'FRAP Curve Figures',['NormFRAP_' dname '.png']))
    end
end

end
