function pinit = poly_param_gen(fname,ncon,nblch,pstart,pix,folder)

bit = 16;

nfas = ncon+nblch;
tyvec = zeros(nfas+1,1)+2;
tyvec(2:ncon+1) = 1;
tyvec(1) = 0;

% generate initial polygons & initial parameters
f = file_search(fname,folder);
img = double(imread(f{1}));
img(img > 2^bit-1) = 0;
[xdim,ydim] = size(img);
img = (255+0.999).*(img-min(min(img)))/(max(max(img))-min(min(img)));
pinit = zeros(nfas,3);
for i = 1:nfas+1
    type = tyvec(i);
    p = get_FA_poly(img,type); %#ok<NASGU>
    name = FRAP_naming(f{1},i,tyvec,ncon);
    if i ~=1
        [p,pinit(i-1,:)] = FIND_FA(img,p,pstart,xdim,ydim,pix);
    end
    save(fullfile(pwd,folder,name),'p','-ascii')
end
% pinit = mean(pinit);
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