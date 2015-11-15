clear
close all
clc

% get general info about the data sets
Get_FRAP_Info

% for each data set, check if need to generate polygons & parameters -> do it
Init_Poly_Gen

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