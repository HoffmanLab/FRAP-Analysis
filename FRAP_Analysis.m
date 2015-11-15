for k = 1:numexp
    % for k = 2
    frap_files = file_search([expcell{k} '\w+_t66.TIF'],keywords.folder);
    pinit_mat = load(fullfile(keywords.folder,'Accessory Files',['pinit_' expcell{k} '.txt']));
    ncon = load(fullfile(keywords.folder,'Accessory Files',['ncon_' expcell{k} '.txt']));
    nblch = load(fullfile(keywords.folder,'Accessory Files',['nblch_' expcell{k} '.txt']));
    pinit = mat2cell(pinit_mat,ncon+nblch,3);
    analyze = load(fullfile(keywords.folder,'Accessory Files',['analyze_' expcell{k} '.txt']));
    recf = load(fullfile(keywords.folder,'Accessory Files',['recf_' expcell{k} '.txt']));
    pix = load(fullfile(keywords.folder,'Accessory Files',['pix_' expcell{k} '.txt']));
    for i = 1:length(frap_files)
        if analyze(i) ~=0
            keywords.recf = recf(i);
            keywords.pix = pix(i);
%             FA_FRAP(strrep(frap_files{i},'t66','t\d+'),[ncon(i) nblch(i) reshape(pinit{i,1},1,numel(pinit{i,1}))],keywords)
            s = strfind(frap_files{i},'t');
            s = s(end);
            nname = frap_files{i}(1:s-2);
%             disp(['FRAP Fitting: ' nname])
            mkdir(keywords.folder,'FRAP Curve Figures')
            FRAP_FIT(nname,keywords.blchtime+2,keywords)
        end
    end
    FA_state(expcell{k},keywords.folder)
end