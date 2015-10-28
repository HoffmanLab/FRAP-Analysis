function FRAP_split(imgexp,folder)

% load images & get poly file structure
% imgs = file_search([imgexp '_t01.TIF'],folder);
imgs = file_search([imgexp '_\d+_t01.TIF'],folder);

for i=1:length(imgs)
    % check that FRAP analysis was completed for this image
    im = double(imread(imgs{i}));
    imset = file_search(strrep(imgs{i},'t01','t\d+'),folder);
    if ~isempty(file_search(['norm_FRAP_blch_' strrep(imgs{i},'_t01.TIF','.dat')],folder))
        % plot image & ask user to click approximately in center
        [r,c] = size(im);
        imgix = repmat(1:r,c,1);
        imgiy = repmat([1:c]',1,r);
        f1 = figure; imagesc(im); title('Select the approximate center of the cell and press Enter')
        [x,y] = getpts(f1);
        close(f1);
        polys = file_search(['poly_FRAP_blch_' strrep(imgs{i},'t01.TIF','t\d+_0.dat')],folder);
        backpolys = file_search(['poly_FRAP_bkg_' strrep(imgs{i},'t01.TIF','t\d+_0.dat')],folder);
        conpolys = file_search(['poly_FRAP_con_' strrep(imgs{i},'t01.TIF','t\d+_0.dat')],folder);
        mean_front_int = zeros(1,length(polys));
        mean_back_int = zeros(1,length(polys));
        backv = zeros(1,length(polys));
        conv = zeros(1,length(polys));
        for j = 1:length(polys)
            % load polygon, find center & minor axis end points
            p = load(polys{j});
            in = inpolygon(imgix,imgiy,p(:,1),p(:,2));
            [y1,x1] = ind2sub([r,c],find(in==1));
            cx = mean(x1);
            cy = mean(y1);
            uxx = sum((x1-cx).^2)./length(x1);
            uyy = sum((y1-cy).^2)./length(y1);
            uxy = sum((x1-cx).*(y1-cy))./length(x1);
            MM = [uxx uxy; uxy uyy];
            [U,S,V] = svd(MM);
            W = V(:,1)/sign(V(1,1));
            H = V(:,2);
            aminx1 = cx+W(1);
            aminy1 = cy+H(1);
            aminx2 = cx-W(1);
            aminy2 = cy-H(1);
            % decide which side is front/back, decide which points are front/back
            side = sign((aminx1-aminx2)*(y1-aminy2)-(aminy1-aminy2)*(x1-aminx2));
            cdist = sqrt((x1-x).^2+(y1-y).^2);
            in_side = side(cdist == min(cdist));
            % avg front & back intensities
            y_in_side = y1(side == in_side);
            x_in_side = x1(side == in_side);
            ind_in_side = sub2ind([r,c],y_in_side,x_in_side);
            img = double(imread(imset{j}));
            mean_front_int(j) = mean(img(ind_in_side));
            y_out_side = y1(side ~= in_side);
            x_out_side = x1(side ~= in_side);
            ind_out_side = sub2ind([r,c],y_out_side,x_out_side);
            mean_back_int(j) = mean(img(ind_out_side));
            % get background & control intensities
            backp = load(backpolys{1});
            in = inpolygon(imgix,imgiy,backp(:,1),backp(:,2));
            masked = img.*in;
            backv(j) = mean(masked(masked~=0));
            conp = load(conpolys{j});
            in = inpolygon(imgix,imgiy,conp(:,1),conp(:,2));
            masked = img.*in;
            conv(j) = mean(masked(masked~=0));
        end
        % normalize to control & background adhesions
        bs_front_blch = mean_front_int - backv;
        bs_back_blch = mean_back_int - backv;
        bscon = conv - backv;
        bscon = bscon/mean(bscon(1:4));
        norm_front_blch = bs_front_blch./bscon;
        norm_front_blch = norm_front_blch/mean(norm_front_blch(1:4));
        norm_back_blch = bs_back_blch./bscon;
        norm_back_blch = norm_back_blch/mean(norm_back_blch(1:4));
        mean_split = [norm_front_blch; norm_back_blch];
        % output file with 2 columns: front & back intensities
        save(fullfile(pwd,folder,['split_FRAP_blch_' strrep(imgs{i},'_t01.TIF','.dat')]),'mean_split','-ascii')
    end
end
end