% clear
close all
clc

keywords.folder = input('Enter the name of the folder containing your images: ','s');
numexp = input('How many constructs did you use? ');
expcell = cell(1,numexp);
for k = 1:numexp
    expcell{k} = input('Enter an experimental group name: ','s');
end
keywords.blchtime = input('What was the last frame before bleaching began? ');
keywords.ntime = input('What is the last time point? ');

% for FRAP_FIT
keywords.dt = input('How many seconds between each time point? ');
keywords.showfit = 1;

% for each data set, check if need to generate polygons & parameters -> do it

keywords.read = input('Use previously generated polygons? Enter (1) if yes, (0) if no: ');
if ~keywords.red
    for k = 1:numexp
        frap_files = file_search([expcell{k} '\w+_t66.TIF'],keywords.folder);
        analyze = zeros(1,1);
        ncon = zeros(1,1);
        nblch = zeros(1,1);
        pinit = cell(1,1);
        for i = 1:length(frap_files)
            disp(frap_files{i})
            analyze(i) = input('Would you like to analyze this set? Enter (1) if yes, (0) if no: ');
            if analyze(i) ~= 0
                keywords.pix = input('How many pixels to zoom out each time (5 is default)? ');
                ncon(i) = input('How many control adhesions to select? ');
                nblch(i) = input('How many adhesions did you bleach? ');
                pstart = input('What FA parameters to start with? ');
                pinit{i,1} = poly_param_gen(strrep(frap_files{i},'t66','t\d+'),ncon(i),nblch(i),pstart,keywords.pix,keywords.folder);
            end
        end
        pinit_mat = cell2mat(pinit);
        save(fullfile(keywords.folder,['pinit_' expcell{k} '.dat']),'pinit_mat','-ascii')
        save(fullfile(keywords.folder,['ncon_' expcell{k} '.dat']),'ncon','-ascii')
        save(fullfile(keywords.folder,['nblch_' expcell{k} '.dat']),'nblch','-ascii')
        save(fullfile(keywords.folder,['analyze_' expcell{k} '.dat']),'analyze','-ascii')
    end
end


% run through FA_FRAP & FRAP_FIT for each -> take out parts that generate
% initial polygons
for k = 1:numexp
    frap_files = file_search([expcell{k} '\w+_t66.TIF'],keywords.folder);
    pinit_mat = load(fullfile(keywords.folder,['pinit_' expcell{k} '.dat']));
    ncon = load(fullfile(keywords.folder,['ncon_' expcell{k} '.dat']));
    nblch = load(fullfile(keywords.folder,['nblch_' expcell{k} '.dat']));
    pinit = mat2cell(pinit_mat,ncon+nblch,3);
    analyze = load(fullfile(keywords.folder,['analyze_' expcell{k} '.dat']));
    for i = 1:length(frap_files)
        if analyze(i) ~=0
            keywords.recf = 0.4;
            FA_FRAP_RECENT(strrep(frap_files{i},'t66','t\d+'),[ncon(i) nblch(i) reshape(pinit{i,1},1,numel(pinit{i,1}))],keywords)
            s = strfind(frap_files{i},'t');
            s = s(end);
            nname = frap_files{i}(1:s-2);
            disp('FRAP Fitting')
            FRAP_FIT_NEW(nname,keywords.blchtime+2,keywords)
        end
    end
    %     FRAP_split(expcell{k},keywords.folder)
    %     disp('FRAP Split Fitting')
    %     FRAP_FIT_SPLIT(expcell{k},keywords.blchtime+2,keywords)
    FA_state(expcell{k},keywords.folder)
end