clear
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

for k = 1:numexp
    frap_files = file_search([expcell{k} '\w+_t66.TIF'],keywords.folder);
    for i = 1:length(frap_files)
        fprintf('\nCurrent image set: %s\n',frap_files{i})
        skip = input('Would you like to skip this set? Enter (1) if yes, (0) if no: ');
        if skip == 0
            keywords.read = input('Use previously generated polygons? Enter (1) if yes, (0) if no: ');
            keywords.pix = input('How many pixels to zoom out each time (5 is default)? ');
            ncon = input('How many control adhesions would you like to select? ');
            nblch = input('How many adhesions did you bleach? ');
            pinit = input('What FA parameters would you like to start with? ');
            FA_FRAP_NEW(strrep(frap_files{i},'t66','t\d+'),ncon,nblch,pinit,keywords);
            s = strfind(frap_files{i},'t');
            s = s(end);
            nname = frap_files{i}(1:s-2);
            FRAP_FIT(['norm_FRAP_blch_' nname '.dat'],keywords.blchtime+2,keywords)
        end
    end
end