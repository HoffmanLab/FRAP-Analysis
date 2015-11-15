function FRAP_PostProcess_Con(exp,folder)

% Outputs file w/ information about Control FA over time
% Column 1: time point
% Column 2: geometric centroid x
% Column 3: geometric centroid y
% Column 4: area
% Column 5: axis ratio
% Column 6: intensity centroid x
% Column 7: intensity centroid y
% Column 8: x position of max intensity
% Column 9: y position of max intensity
% Column 10: instantaneous velocity

% find relevant files
init_img = [exp '_t01.TIF'];
con_polys = file_search(['poly_FRAP_con_' exp '_\w+_0.dat'],folder);
[r,c] = size(imread(init_img));
len = length(con_polys);
output = zeros(len,10);

% calculate centriod, area, axis ratio, etc at each point
imgix = repmat(1:r,c,1);
imgiy = repmat([1:c]',1,r);

for i = 1:len
    output(i,1) = 5*(i-1);
    p = load(fullfile(folder,'FRAP Poly Files',con_polys{i}));
    if i == 1
        pxmin = floor(min(p(:,1)));
        pxmax = ceil(max(p(:,1)));
        pymin = floor(min(p(:,2)));
        pymax = ceil(max(p(:,2)));
    end
    eval(sprintf('img = double(imread(fullfile(folder,[exp ''_t%02d.TIF''])));',i))
    in = inpolygon(imgix,imgiy,p(:,1),p(:,2));
    [y1,x1] = ind2sub([r,c],find(in==1));
    output(i,3) = mean(x1); % centroid x
    output(i,2) = mean(y1); % centroid y
    output(i,4) = length(y1); % area
    uxx = sum((x1-output(i,3)).^2)./length(x1);
    uyy = sum((y1-output(i,2)).^2)./length(y1);
    uxy = sum((x1-output(i,3)).*(y1-output(i,2)))./length(x1);
    qrot = sqrt((uxx-uyy).^2+4*(uxy.^2));
    mjra = sqrt(2)*sqrt(uxx+uyy+qrot);
    mnra = sqrt(2)*sqrt(uxx+uyy-qrot);
    output(i,5) = mjra/mnra; % axis ratio
    output(i,6) = sum(img(in).*y1)./sum(img(in)); % intensity centroid x
    output(i,7) = sum(img(in).*x1)./sum(img(in)); % intensity centroid y
    [suby,subx] = find(img == max(img(in)));
    subin = inpolygon(subx,suby,p(:,1),p(:,2));
    testy = suby(subin);
    output(i,8) = testy(1); % location of max intensity x
    testx = subx(subin);
    output(i,9) = testx(1); % location of max intensity y
    if i > 1
        output(i,10) = sqrt((output(i,2)-output(i-1,2))^2+(output(i,3)-output(i-1,3))^2)/5; % velocity pixels/s
    end
%     if rem(i-1,5) == 0
%         subimg = img(floor(min(p(:,1))):ceil(max(p(:,1))),floor(min(p(:,2))):ceil(max(p(:,2))));
% %         figure; 
%         subplot(ceil(len/5),1,(i-1)/5+1); 
%         imagesc(img,[0 max(max(subimg))]);
%         hold on; patch(p(:,1),p(:,2),'k','FaceColor','none','EdgeColor','k')
%         axis([pxmin-50 pxmax+50 pymin-30 pymax+30])
%     end
end

save(fullfile(folder,'FRAP PostProcess Files',['con_' exp '_poly_data.txt']),'output','-ascii')
end