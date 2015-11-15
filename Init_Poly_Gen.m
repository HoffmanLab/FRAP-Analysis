keywords.read = input('Use previously generated polygons? Enter (1) if yes, (0) if no: ');
if ~keywords.read
    for k = 1:numexp
        frap_files = file_search([expcell{k} '\w+_t66.TIF'],keywords.folder);
        analyze = zeros(length(frap_files));
        ncon = zeros(length(frap_files));
        nblch = zeros(length(frap_files));
        conFA = zeros(length(frap_files));
        blchFA = zeros(length(frap_files));
        recf = zeros(length(frap_files));
        pinit = cell(1,1);
        for i = 1:length(frap_files)
            disp(frap_files{i})
            analyze(i) = input('Would you like to analyze this set? Enter (1) if yes, (0) if no: ');
            if analyze(i) ~= 0
                keywords.pix = input('How many pixels to zoom out each time (5 is default)? ');
                ncon(i) = input('How many control adhesions to select? ');
                nblch(i) = input('How many adhesions did you bleach? ');
                pstart = input('What FA parameters to start with? ');
                recf(i) = input('Enter the recovery fraction (0.25 is default): ');
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
        save(fullfile(keywords.folder,'Accessory Files',['recf_' expcell{k} '.txt']),'recf','-ascii')
    end
end
