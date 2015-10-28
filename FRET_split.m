function FRET_split(expname,bfile,folder)

% This function adds columns to an existing blb file that has had closed regions
% selected to represent the "front" (closest to cell center) and "rear"
% (closest to cell edge) as determined by the center of mass of the acceptor 
% intensity. It also calculates the geometric centroid and the center of
% mass of the FRET signal.

% Sample Call:
%   FRET_split('VinTS','blb_anl_pre_VinTS.txt','120914 VinTS MEF FRET')

% Inputs:
%   expname - expression representing the experimental group to analyze
%   bfile - name of the blob file matching the experimental data
%   folder - folder containing the images

% Outputs:
%   creates new blob file starting with "split_" with the following
%   additional columns:
%       mean FRET signal "front"
%       mean FRET signal "rear"
%       mean donor signal "front"
%       mean donor signal "rear"
%       mean acceptor signal "front"
%       mean acceptor signal "rear"
%   the FAs are divided into halves based on their geometric centroid

fa_files = file_search(['fa_bsa_pre_' expname '_\d+_w1Venus.TIF'],folder);
fret_files = file_search(['masked_on_Venus_eff_pre_' expname '_\d+_w2TVFRET.TIF'],folder);
bsa_files = file_search(['masked_on_Venus_bsa_pre_' expname '_\d+_w1Venus.TIF'],folder);
% bsd_files = file_search(['bsd_pre_' expname '_\d+_w3Teal.TIF'],folder);
b = load(bfile);
newfile = zeros(size(b,1),7);
blbcol = 30; % SHOULDN'T HARD CODE THESE...
imgcol = 31;
cellcol = 32;
szcol = 27;

% [r,c] = size(b);
% b(:,end+1:end+6) = zeros(size(b,1),6);

imgs = unique(b(:,imgcol));
nimgs = length(imgs);

if length(imgs) ~= length(fa_files)
    disp('Number of image files doesn''t match image IDs')
else
    for i = 1:nimgs
        fa_img = double(imread(fa_files{i}));
        fret_img = double(imread(fret_files{i}));
        bsa_img = double(imread(bsa_files{i}));
%         bsd_img = double(imread(bsd_files{i}));
        btemp = b(b(:,imgcol) == imgs(i),:);
        cells = unique(btemp(:,cellcol));
        ncells = length(cells);
        for j = 1:ncells
            if cells(j) ~=0
                btemp2 = btemp(btemp(:,cellcol) == cells(j),:);
                blobs = unique(btemp2(:,blbcol));
                nblobs = length(blobs);
                for k = 1:nblobs
                    index = find(btemp2(:,blbcol) == blobs(k));
                    tot_index = find((b(:,imgcol) == imgs(i)) & (b(:,cellcol) == cells(j)) & (b(:,blbcol) == blobs(k)));
                    sz = btemp2(index,szcol);
                    if blobs(k)~=0 && sz > 3
%                         x = btemp2(index,18);
%                         y = btemp2(index,19);
                        cx = btemp2(index,1);
                        cy = btemp2(index,2);
                        [y1,x1] = find(fa_img == blobs(k));
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
                        side = sign((aminx1-aminx2)*(y1-aminy2)-(aminy1-aminy2)*(x1-aminx2));
                        cdist = sqrt((x1-cx).^2+(y1-cy).^2);
                        in_side = side(cdist == min(cdist));
                        in_side = in_side(1);
                        y_in_side = y1(side == in_side);
                        x_in_side = x1(side == in_side);
                        mean_fret_front = mean(mean(fret_img(y_in_side,x_in_side)));
                        mean_bsa_front = mean(mean(bsa_img(y_in_side,x_in_side)));
%                         mean_bsd_front = mean(mean(bsd_img(y_in_side,x_in_side)));
                        y_out_side = y1(side ~= in_side);
                        x_out_side = x1(side ~= in_side);
                        mean_fret_back = mean(mean(fret_img(y_out_side,x_out_side)));
                        mean_bsa_back = mean(mean(bsa_img(y_out_side,x_out_side)));
%                         mean_bsd_back = mean(mean(bsd_img(y_out_side,x_out_side)));
                        
                        add_cols = [mean_fret_front mean_fret_back mean_bsa_front mean_bsa_back];
                        newfile(tot_index,1) = imgs(i);
                        newfile(tot_index,2) = cells(j);
                        newfile(tot_index,3) = blobs(k);
                        newfile(tot_index,4:end) = add_cols;
                    end
                end
            end
        end
    end
end

save(fullfile(pwd,folder,['FRETsplit_' bfile '.txt']),'newfile','-ascii')

col_names = {'Image ID','Cell ID','Blob ID','FRET Proximal','FRET Distal','Venus Proximal','Venus Distal'};
cell_final_data = num2cell(newfile);
cell_final_file = [col_names; cell_final_data];
xlswrite(fullfile(pwd,folder,['FRETsplit_' bfile '.xlsx']),cell_final_file)