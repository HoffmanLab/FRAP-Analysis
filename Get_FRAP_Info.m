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