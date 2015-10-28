function plot_motion(exp,blbn,folder)

% find geometric centroid @ initial point
b = file_search(['con_' exp '_poly_data.txt'],folder);
data = load(b{1});
geox = data(:,2);
geoy = data(:,3);

% find max venus over time
maxvx = data(:,8);
maxvy = data(:,9);

% find orientation of blob
u = strfind(exp,'_');
short = exp(1:u(1)-1);
a = file_search(['blb_anl_pre_' short '_FRET.txt'],folder);
blb = load(a{1});
blbcol = 30; % SHOULDN'T HARD CODE THESE...
imgcol = 31;
cellcol = 32;
orcol = 29;
imgn = str2num(exp(u(2)+1:end));
celln = 1;
blbrow = find(blb(:,imgcol)==imgn & blb(:,cellcol) == celln & blb(:,blbcol) == blbn);
or = pi - blb(blbrow,orcol);
cellx = blb(blbrow,end-1);
celly = blb(blbrow,end);

% rotate points according to orientation
rot_mat = [cos(or) -1*sin(or) ; sin(or) cos(or)];
rmaxcoords = rot_mat * [maxvy';maxvx'];
rcent = rot_mat * [geoy';geox'];
rcell = rot_mat * [cellx;celly];

% plot motion along major vs. minor axis, make positive toward cell edge
tocentx = sign(rcell(1) - rcent(1));
tocenty = sign(rcell(2) - rcent(2));

rmaxx = rmaxcoords(1,:)' - rcent(1,:)';
rmaxy = rmaxcoords(2,:)' - rcent(2,:)';

matchx = double(tocentx == sign(rmaxx));
matchx(matchx == 0) = -1;
matchx = -1*matchx;
matchy = double(tocenty == sign(rmaxy));
matchy(matchy == 0) = -1;
matchy = -1*matchy;

rmaxx = abs(rmaxx).*matchx;
rmaxy = abs(rmaxy).*matchy;

figure; plot([0:5:325],rmaxx,'k','LineWidth',3); axis([0 325 -50 50])
figure; plot([0:5:325],rmaxy,'k','LineWidth',3); axis([0 325 -50 50])

rmaxcoords = [rmaxx rmaxy];

save(fullfile(pwd,folder,['Venus Max Coords_' exp '.txt']),'rmaxcoords','-ascii')

end