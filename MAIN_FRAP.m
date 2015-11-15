clear
close all
clc

keywords.folder = input('Enter the name of the folder containing your images: ','s');
addpath(keywords.folder)
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
if ~keywords.read
    for k = 1:numexp
        frap_files = file_search([expcell{k} '\w+_t66.TIF'],keywords.folder);
        analyze = zeros(length(frap_files));
        ncon = zeros(length(frap_files));
        nblch = zeros(length(frap_files));
        conFA = zeros(length(frap_files));
        blchFA = zeros(length(frap_files));
        pinit = cell(1,1);
        for i = 1:length(frap_files)
            disp(frap_files{i})
            analyze(i) = input('Would you like to analyze this set? Enter (1) if yes, (0) if no: ');
            if analyze(i) ~= 0
                keywords.pix = input('How many pixels to zoom out each time (5 is default)? ');
                ncon(i) = input('How many control adhesions to select? ');
                nblch(i) = input('How many adhesions did you bleach? ');
                pstart = input('What FA parameters to start with? ');
                conFA(i) = input('Enter the FA# for the control FA: ');
                blchFA(i) = input('Enter the FA# for the bleached FA: ');
                pinit{i,1} = poly_param_gen(strrep(frap_files{i},'t66','t\d+'),ncon(i),nblch(i),pstart,keywords.pix,keywords.folder);
            end
        end
        pinit_mat = cell2mat(pinit);
        mkdir(keywords.folder,'Accessory Files')
        save(fullfile(keywords.folder,'Accessory Files',['pinit_' expcell{k} '.txt']),'pinit_mat','-ascii')
        save(fullfile(keywords.folder,'Accessory Files',['ncon_' expcell{k} '.txt']),'ncon','-ascii')
        save(fullfile(keywords.folder,'Accessory Files',['nblch_' expcell{k} '.txt']),'nblch','-ascii')
        save(fullfile(keywords.folder,'Accessory Files',['analyze_' expcell{k} '.txt']),'analyze','-ascii')
        save(fullfile(keywords.folder,'Accessory Files',['conFAs_' expcell{k} '.txt']),'conFA','-ascii')
        save(fullfile(keywords.folder,'Accessory Files',['blchFAs_' expcell{k} '.txt']),'blchFA','-ascii')
    end
end


% run through FA_FRAP & FRAP_FIT for each -> take out parts that generate
% initial polygons
for k = 1:numexp
    % for k = 2
    frap_files = file_search([expcell{k} '\w+_t66.TIF'],keywords.folder);
    pinit_mat = load(fullfile(keywords.folder,['pinit_' expcell{k} '.dat']));
    ncon = load(fullfile(keywords.folder,['ncon_' expcell{k} '.dat']));
    nblch = load(fullfile(keywords.folder,['nblch_' expcell{k} '.dat']));
    pinit = mat2cell(pinit_mat,ncon+nblch,3);
    analyze = load(fullfile(keywords.folder,['analyze_' expcell{k} '.dat']));
    for i = 1:length(frap_files)
        if analyze(i) ~=0
            keywords.recf = 0.25;
            FA_FRAP(strrep(frap_files{i},'t66','t\d+'),[ncon(i) nblch(i) reshape(pinit{i,1},1,numel(pinit{i,1}))],keywords)
            s = strfind(frap_files{i},'t');
            s = s(end);
            nname = frap_files{i}(1:s-2);
            disp(['FRAP Fitting: ' nname])
            FRAP_FIT(nname,keywords.blchtime+2,keywords)
        end
    end
    FA_state(expcell{k},keywords.folder)
end