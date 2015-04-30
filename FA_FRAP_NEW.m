function FA_FRAP_NEW(fname,ncon,nblch,pinit,keywords)

%% set default values of keywords fields
if ~isfield(keywords,'movie')
    keywords.movie = 0;
end
if ~isfield(keywords,'read')
    keywords.read = 0;
end
if ~isfield(keywords,'pix')
    keywords.pix = 5;
end
if ~isfield(keywords,'bit')
    keywords.bit = 16;
end
if ~isfield(keywords,'ntime')
    keywords.ntime = 66;
end
if ~isfield(keywords,'blchtime')
    keywords.blchtime = 4;
end

%% initialize vectors to hold values
backv = zeros(keywords.ntime,1);
conv = zeros(keywords.ntime,ncon);
blchv = zeros(keywords.ntime,nblch);
nfas = ncon + nblch;
tyvec = zeros(nfas+1,1)+2;
tyvec(2:ncon+1) = 1;
tyvec(1) = 0;

%% get initial polygons for each & save
f = file_search(fname,keywords.folder);
img = double(imread(f{1}));
img(img > 2^keywords.bit-1) = 0;
[xdim,ydim] = size(img);
img = (255+0.999).*(img-min(min(img)))/(max(max(img))-min(min(img)));
if ~keywords.read
    for i = 1:nfas+1
        type = tyvec(i);
        p = get_FA_poly(img,type); %#ok<NASGU>
        name = FRAP_naming(f{1},i,tyvec,ncon);
        save(fullfile(pwd,keywords.folder,name),'p','-ascii')
    end
end

%% initialize indices for inpolygon

imgix = repmat(1:xdim,ydim,1);
imgiy = repmat([1:ydim]',1,xdim); %#ok<NBRAK>

%% loop through time points
for i = 1:keywords.ntime
    img = double(imread(f{i}));
    for j = 1:length(tyvec)
        type = tyvec(j);
        if type == 0
           bname = FRAP_naming(f{1},1,tyvec,ncon);
           p = load(bname);
           in = inpolygon(imgix,imgiy,p(:,1),p(:,2));
           masked = img.*in;
           backv(i) = mean(masked(masked~=0));
        else
           bname = FRAP_naming(f{i},j,tyvec,ncon);
           p = load(bname);
           if rem(i-1,5) == 0 % if this is every 5th frame, starting with #1
               [p,pinit] = FIND_FA(img,p,pinit,xdim,ydim,keywords.pix);
           end
           in = inpolygon(imgix,imgiy,p(:,1),p(:,2));  
           if i ~= keywords.ntime
              bname = FRAP_naming(f{i},j,tyvec,ncon);
              save(fullfile(pwd,keywords.folder,bname),'p','-ascii')
              wname = FRAP_naming(f{i+1},j,tyvec,ncon);
              save(fullfile(pwd,keywords.folder,wname),'p','-ascii')
           end
           masked = img.*in;
           if type == 1
               pos = j-1;
               conv(i,pos) = mean(masked(masked~=0));
           elseif type == 2
               pos = j-ncon-1;
               blchv(i,pos) = mean(masked(masked~=0));
           end
        end
    end
end

%% write out raw data
nname = f{i};
s = strfind(nname,'t');
s = s(end);
nname = nname(1:s-2);
save(fullfile(pwd,keywords.folder,['FRAP_bkg_' nname '.dat']),'backv','-ascii')
save(fullfile(pwd,keywords.folder,['FRAP_blch_' nname '.dat']),'blchv','-ascii')
save(fullfile(pwd,keywords.folder,['FRAP_con_' nname '.dat']),'conv','-ascii')

%% normalize data & output final vector
bsblch = blchv - backv;
bscon = zeros(keywords.ntime,ncon);
for i = 1:ncon
   bscon(:,i) = conv(:,i) - backv;
   bscon(:,i) = bscon(:,i)/mean(bscon(1:keywords.blchtime,i));
end

if ncon > 1
   bscon = mean(bscon,2); 
end

normblch = bsblch./bscon;
normblch = normblch/mean(normblch(1:keywords.blchtime,1)); %#ok<NASGU>
save(fullfile(pwd,keywords.folder,['norm_FRAP_blch_' nname '.dat']),'normblch','-ascii')
end

function p = get_FA_poly(img,type)
    switch type
        case 2
            text='Draw a approximate polygon around the bleached FA';
        case 1
            text='Draw a approximate polygon around the control FA';
        case 0
            text='Draw a polygon in the background';
    end
    cont_image(img)
    disp(text)
    [~, px, py] = roipoly;
    p = [px py];
    close 
end

function p = mask_2_con(mask)
    figure;
    [C,h] = contour(mask,1,'Visible','off'); %#ok<NASGU>
    close;
    a = [];
    aloc = [];
    a(1) = C(2,1);
    aloc(1) = 1;
    count = 1;
    while(sum(a)+count < length(C))
        count = count+1;
        aloc(count) = sum(a)+count; %#ok<AGROW>
        a(count) = C(2,sum(a)+count); %#ok<AGROW>
    end
    bigai = find(a == max(a));
    biga = a(bigai);
    locbiga = aloc(bigai);
    p = [C(:,locbiga+1:biga+1)';C(:,locbiga+1)'];
end

function [p,pinit] = FIND_FA(img,p,pinit,xdim,ydim,pix)
   xmin=floor(min(p(:,1)));
   xmax=ceil(max(p(:,1)));
   ymin=floor(min(p(:,2)));
   ymax=ceil(max(p(:,2)));
   temp=img(ymin-2*pix:ymax+2*pix,xmin-2*pix:xmax+2*pix);
   figure; imagesc(img,[0 max(max(temp))]); hold on;
   patch(p(:,1),p(:,2),'w','EdgeColor','k','FaceColor','none')
   axis([xmin-pix,xmax+pix,ymin-pix,ymax+pix])
   choice = menu('Would you like to continue with this polygon or generate a new one?','Continue','Re-generate');
   close;
   if choice == 2
       pinit = ParameterSelectorFunction(temp,[0,100],[0,10000],[0,200],pinit);
       fao = water(temp,pinit);
       u = unique(fao);
       blength = 0;
       big = 0;
       for i = 2:length(u)
           l = length(find(fao == u(i)));
           if l > blength
               big = u(i);
               blength = l;
           end
       end
       temp(:,:) = 0;
       temp(fao == big) = 1;
       infa = zeros(xdim,ydim);
       infa(ymin-2*pix:ymax+2*pix,xmin-2*pix:xmax+2*pix) = temp;
       p = mask_2_con(infa);
   end
end