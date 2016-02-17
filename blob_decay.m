function blob_decay(fname,bimg,outname,varargin)

folder = varargin{end};
f = file_search(fname,folder);
nimg = length(f);
mask = single(imread(bimg));
if nargin == 4
    blbs = load(varargin{1});
    nblbs = length(blbs);
else
    blbs = 1:max(max(mask));
    nblbs = max(max(mask));
end

if nimg == 0
    disp('No files found')
else
    blbarr = zeros(nblbs,nimg+1);
    for i = 1:nimg
        img = single(imread(f{i}));
        sz = size(img);
        crop = [round(0.0246*sz(1)) round(0.0246*sz(1)) round(0.9509*sz(1)) round(0.9509*sz(1))];
        img = imcrop(img,crop);
        for j = 1:nblbs
            blbarr(j,1) = blbs(j);
            blbarr(j,i+1) = mean(img(mask == blbs(j)));
        end
    end
end
save(fullfile(folder,outname),'blbarr','-ascii')

end