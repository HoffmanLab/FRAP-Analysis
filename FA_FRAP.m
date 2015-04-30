function FA_FRAP(fname,parray,keywords)

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
    keywords.bit = 8;
end

f = file_search(fname,keywords.sourcefile);

nblch = parray(1);
ncon = parray(2);
parray = parray(3:end);
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

% Read in initial image and get the polygons for the image
img = double(imread(f{1}));
img(img > 2^keywords.bit-1) = 0;
[xdim,ydim] = size(img);
img = (255+0.999).*(img-min(min(img)))/(max(max(img))-min(min(img)));
if ~keywords.read
    for i = 1:nfas+1
        type = tyvec(i);
        p = get_FA_poly(img,type,keywords.shift);
        ti = 1;
        name = FRAP_naming(f{ti},ti,i,tyvec,ncon);
        save([pwd '/' keywords.destfile '/' name] ,'p','-ascii')
    end
end

% For each type of polygon, find the average value
imgix = repmat(1:xdim,ydim,1);
imgiy = repmat([1:ydim]',1,xdim);

count = 0;

for i = 1:nele
    img = double(imread(f{i}));
    for j = 1:length(tyvec);
        type = tyvec(j);
        if type == 0 % Background first
            bname = FRAP_naming(f{1},1,1,tyvec,ncon);
            p = load(bname);
            in = inpolygon(imgix,imgiy,p(:,1),p(:,2));
        else % Control or Bleached
            bname = FRAP_naming(f{i},i,j,tyvec,ncon);
            p = load(bname);
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
                if type == 1 || (type == 2 && i <= (keywords.blchtime + 1))
                    in = FIND_FA(img,p,parray,xdim,ydim);
                elseif type == 2 && i > keywords.blchtime
                    in = FIND_FA(img,p,keywords.newparray,xdim,ydim); 
                    count = count+1;
                end
                intemp = zeros(xdim,ydim);
                intemp(in) = 1;
                p = mask_2_con(intemp); 
                in = inpolygon(imgix,imgiy,p(:,1),p(:,2));
            end
        end
        % Store data into array
        if type == 0
            masked = img.*in;
            backv(i) = mean(masked(masked~=0));
        elseif type == 2
            pos = j-ncon-1;
            masked = img.*in;
            blchv(pos,i) = mean(masked(masked~=0));
        elseif type == 1
            pos = j-1;
            masked = img.*in;
            conv(pos,i) = mean(masked(masked~=0));
        end
        if i ~= nele
            wname = FRAP_naming(f{i+1},i+1,j,tyvec,ncon);
            save([pwd '/' keywords.destfile '/' wname],'p','-ascii')
        end
    end
end

% Write out data
nname = f{i};
s = strfind(nname,'.');
s = s(end);
nname = nname(1:s-1);

save(fullfile(pwd,keywords.destfile,['FRAP_bkg_' nname '.dat']),'backv','-ascii')
save(fullfile(pwd,keywords.destfile,['FRAP_blch_' nname '.dat']),'blchv','-ascii')
save(fullfile(pwd,keywords.destfile,['FRAP_con_' nname '.dat']),'conv','-ascii')

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

save([pwd '/' keywords.destfile '/' 'norm.FRAP_blch_' nname '.dat'],'normblch','-ascii')

end

% SUBFUNCTIONS

function pfa = get_FA_poly(img,type,shift)
switch type
    case 2
        text='Draw a approximate polygon around the bleached FA';
    case 1
        text='Draw a approximate polygon around the control FA';
    case 0
        text='Draw a polygon in the background';
end
if shift
    cont_image(img)
end
if type == 0
    disp(text)
    [pmask, px, py] = roipoly;
    pfa = [px py];
else
    if type == 1
        disp(text)
    end
    if type == 2
        disp(text)
    end
    [appmask, ax, ay]=roipoly;
    app = [ax ay];
    
    if ~shift
        mins=[min(app(1,:)),min(app(2,:))];
        maxs=[max(app(1,:)),max(app(2,:))];
        img2=img(mins(1):maxs(1),mins(2):maxs(2));
        cont_image(img2)
        disp('draw a refined polygon around the FA')
        pfa=roipoly;
        pfa(1,:)=pfa(1,:)+mins(1);
        pfa(2,:)=pfa(2,:)+mins(2);
    else
        pfa=app;
    end
end
close
end

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