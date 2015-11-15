function FA_FRAP_RECENT(fname,parray,keywords)

% FA_FRAP(imgstruct,[1,1,0,65,30],keywords)

% A function to analyze FRAP data, specifically for moving blobs (focal
% adhesions). It uses the water program to automatically find the FAs of
% interest. User must supply approximate region for bleached FAs,
% unbleached FAs, and a background region outside the cell.
% Can handle up to 9 bleached and unbleached regions.
% Uses the external function FRAP_naming to handle the naming of all the
% files.

% INPUTS:
% imfile - a structure containing the .tif files extracted from the orinial
%       lsm file
% parray - an array of numerical parameters used to control the program
%       parray(1) - number of FAs bleached in the picture
%       parray(2) - number of control FAs to find
%       parray(3) - width of smoothing filter for water (0 for no
%       smoothing, 25 for high background)
%       parray(4) - minimal pixel intensity for water
%       parray(5) - minimum area for merger for water (set to something quite large)
% keywords - a structure containing the following fields:
%   movie - if set to 1, show the images with polygons for each FA and
%       background
%   dt - contains the timestep between the frames of the movie
%   read - if set to 1, read in previously drawn polygons
%   recf - contains the percentage that a FA must recover before a new FA
%       is located
%   pix - contains the number of pixels to enlarge the FA by when shifting,
%       should be larger than the distance traveled
%   shift - if set to 1, initiates auto finding of FAs using the water
%       routine. if set to zero, uses user-defined
%   silent - if set to 1, turns off warnings

if ~isfield(keywords,'dt')
    keywords.dt = 1;
end
if ~isfield(keywords,'pix')
    keywords.pix = 5;
end
if ~isfield(keywords,'recf')
    keywords.recf = 0.5;
end
if ~isfield(keywords,'bit')
    keywords.bit = 16;
end

f = file_search(fname,keywords.folder);

nblch = parray(2);
ncon = parray(1);
parray = reshape(parray(3:end),nblch+ncon,3);
nele = length(f);
backv = zeros(1,nele);
blchv = zeros(nblch,nele);
conv = zeros(ncon,nele);
nfas = nblch+ncon;

% 0 is bkg, 1 is con, 2 is bleached
% Create vector containing the type of polygon we will need
tyvec = zeros(1,nfas+1)+2;
tyvec(2:ncon+1) = 1;
tyvec(1) = 0;

for i = 1:nele
    img = double(imread(f{i}));
    [xdim,ydim] = size(img);
    imgix = repmat(1:xdim,ydim,1);
    imgiy = repmat([1:ydim]',1,xdim);
    for j = 1:length(tyvec);
        type = tyvec(j);
        if type == 0 % Background first
            bname = FRAP_naming(f{1},1,tyvec,ncon);
            p = load(fullfile(keywords.folder,'FRAP Poly Files',bname));
            in = inpolygon(imgix,imgiy,p(:,1),p(:,2));
            masked = img.*in;
            backv(i) = mean(masked(masked~=0));
        else % Control or Bleached
            bname = FRAP_naming(f{i},j,tyvec,ncon);
            p = load(fullfile(keywords.folder,'FRAP Poly Files',bname));
            in = inpolygon(imgix,imgiy,p(:,1),p(:,2));
            tb = mean(mean(img((img.*in)~=0)));
            if type == 2
                move = (tb>(keywords.recf*blchv(1,j-ncon-1)));
            else
                move = 0;
            end
            if type == 1 || move || i == 1
                str = ones(keywords.pix,keywords.pix);
                tmask = zeros(xdim,ydim);
                tmask(in) = 1;
                tmask = imdilate(tmask,str);
                p = mask_2_con(tmask);
                in = FIND_FA(img,p,parray(j-1,:),xdim,ydim);
                intemp = zeros(xdim,ydim);
                intemp(in) = 1;
                p = mask_2_con(intemp); 
                in = inpolygon(imgix,imgiy,p(:,1),p(:,2));
            end
        end
        % Store data into array
        if type == 2
            pos = j-ncon-1;
            masked = img.*in;
            blchv(pos,i) = mean(masked(masked~=0));
        elseif type == 1
            pos = j-1;
            masked = img.*in;
            conv(pos,i) = mean(masked(masked~=0));
        end
        if i ~= nele
            wname = FRAP_naming(f{i+1},j,tyvec,ncon);
            save(fullfile(keywords.folder,'FRAP Poly Files',wname),'p','-ascii')
        end
    end
end

% Write out data
nname = f{i};
s = strfind(nname,'t');
s = s(end);
nname = nname(1:s-2);

mkdir(keywords.folder,'FRAP Curve Files')
save(fullfile(keywords.folder,'FRAP Curve Files',['FRAP_bkg_' nname '.dat']),'backv','-ascii')
save(fullfile(keywords.folder,'FRAP Curve Files',['FRAP_blch_' nname '.dat']),'blchv','-ascii')
save(fullfile(keywords.folder,'FRAP Curve Files',['FRAP_con_' nname '.dat']),'conv','-ascii')

% Normalize & Write Out Data

bsblch = zeros(nblch,nele);
for i = 1:nblch
bsblch(i,:) = blchv(i,:) - mean(backv);
end

bscon = zeros(ncon,nele);
for i = 1:ncon
bscon(i,:) = conv(i,:) - mean(backv);
bscon(i,:) = bscon(i,:)/mean(bscon(i,1:keywords.blchtime));
end

if ncon > 1
   bscon = mean(bscon); 
end

normblch = zeros(nblch,nele);
for i = 1:nblch
normblch(i,:) = bsblch(i,:)./bscon;
normblch(i,:) = normblch(i,:)/mean(normblch(i,1:keywords.blchtime));
end

save(fullfile(keywords.folder,'FRAP Curve Files',['norm_FRAP_blch_' nname '.dat']),'normblch','-ascii')

end

% SUBFUNCTIONS

function infa = FIND_FA(img,pfa,parray,xdim,ydim)

%Create a box
xmin=floor(min(pfa(:,1)));
xmax=ceil(max(pfa(:,1)));
ymin=floor(min(pfa(:,2)));
ymax=ceil(max(pfa(:,2)));
temp=img(ymin:ymax,xmin:xmax);

%Find and aanalyze the adhesion
fao=water(temp,parray);

%Find the largest one
u = unique(fao);
big = 0;
blength = 0;
for i = 2:length(u)
    l = length(find(fao==u(i)));
    if l > blength
        big = u(i);
        blength = l;
    end
end
temp(:,:) = 0;
temp(fao == big) = 1;
temp2 = zeros(xdim,ydim);
temp2(ymin:ymax,xmin:xmax) = temp;
infa = find(temp2);

end

function pfa = mask_2_con(mask)
figure;
[C,h] = contour(mask,1,'Visible','off');
close;
a = [];
aloc = [];
a(1) = C(2,1);
aloc(1) = 1;
count = 1;
while(sum(a)+count < length(C))
    count = count+1;
    aloc(count) = sum(a)+count;
    a(count) = C(2,sum(a)+count);
end
bigai = find(a == max(a));
biga = a(bigai);
locbiga = aloc(bigai);
pfa = [C(:,locbiga+1:biga+1)';C(:,locbiga+1)'];

end