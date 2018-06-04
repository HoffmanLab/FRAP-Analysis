function blob_decay(fname,bimg,outname,varargin)
%% testing
folder = varargin{end};
addpath(folder)
addpath(fullfile(folder,'FA Images'))
f = file_search(fname,folder);
nimg = length(f);
mask = single(imread(bimg));
if nargin == 5
    blbs = load(varargin{1});
    nblbs = length(blbs);
else
    blbs = 1:max(max(mask));
    nblbs = max(max(mask));
end

if nimg == 0
    disp('No files found')
else
    blbarr = zeros(nimg+1,nblbs+1);
    for i = 1:nimg
        img = single(imread(f{i}));
%         sz = size(img);
%         crop = [round(0.0246*sz(1)) round(0.0246*sz(1)) round(0.9509*sz(1)) round(0.9509*sz(1))];
%         img = imcrop(img,crop);
        for j = 1:nblbs
            blbarr(1,j+1) = blbs(j);
            blbarr(i+1,1) = i;
            blbarr(i+1,j+1) = mean(img(mask == blbs(j)));
        end
    end
end

% blb_titles = mat2cell(blbs
save(fullfile(folder,outname),'blbarr','-ascii')
rmpath(folder)
rmpath(fullfile(folder,'FA Images'))

end