for k = 1:numexp
    analyze = load(fullfile(keywords.folder,'Accessory Files',['analyze_' expcell{k} '.txt']));
    conFA = load(fullfile(keywords.folder,'Accessory Files',['conFAs_' expcell{k} '.txt']));
    frap_files = file_search([expcell{k} '\w+_t66.TIF'],keywords.folder);
    for i = 1:length(frap_files)
        if analyze(i) ~=0
            s = strfind(frap_files{i},'t');
            s = s(end);
            nname = frap_files{i}(1:s-2);
            mkdir(keywords.folder,'FRAP PostProcess Files')
            FRAP_PostProcess(nname,keywords.folder)
            FRAP_PostProcess_Con(nname,keywords.folder)
            blbn = conFA(i);      
            mkdir(keywords.folder,'FRAP Venus Max Displacement')
            plot_cents(nname,keywords.folder);
            plot_motion(nname,blbn,keywords.folder);
        end
    end
    mkdir(keywords.folder,'FRAP FA State')
    FA_state(expcell{k},keywords.folder)
    mkdir(keywords.folder,'FRAP Split Curve Files')
    mkdir(keywords.folder,'FRAP Split Curve Figures')
    FRAP_split(expcell{k},keywords.folder)
    FRAP_FIT_SPLIT(expcell{k},keywords.blchtime+2,keywords)
    FRET_split(expcell{k},['blb_anl_rp_' strrep(expcell{k},'FRAP','FRET') '.txt'],keywords.folder)
end