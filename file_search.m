function filecell = file_search(exp,folder)

% a function to find files of a given expression in a particular folder

listing = dir(folder);
listing(1:2)=[];
files = {listing.name};

match = regexp(files,exp);
indi = [];
for i = 1:length(match)
    if match{i}==1
        indi = [indi i];
    end
end
filecell = files(indi);