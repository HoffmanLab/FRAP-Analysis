for k = 1:numexp
    analyze = load(fullfile(keywords.folder,['analyze_' expcell{k} '.dat']));
    frap_files = file_search([expcell{k} '\w+_t66.TIF'],keywords.folder);
    for i = 1:length(frap_files)
        if analyze(i) ~=0
            s = strfind(frap_files{i},'t');
            s = s(end);
            nname = frap_files{i}(1:s-2);
            FRAP_PostProcess(nname,keywords.folder)
            FRAP_PostProcess_Con(nname,keywords.folder)
            plot_cents(nname,keywords.folder);
            blbn = input('Enter the control FA #: ');
            plot_motion(nname,blbn,keywords.folder);
        end
    end
    FRAP_split(expcell{k},keywords.folder)
    disp('FRAP Split Fitting')
    FRAP_FIT_SPLIT(expcell{k},keywords.blchtime+2,keywords)
    FRET_split(expcell{k},['blb_anl_pre_' strrep(expcell{k},'FRAP','FRET') '.txt'],keywords.folder)
end