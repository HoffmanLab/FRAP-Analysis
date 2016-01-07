function plot_motion(exp,blbn,folder)

% find geometric centroid @ initial point
b = file_search(['con_' exp '_poly_data.txt'],folder);
data = load(fullfile(folder,'FRAP PostProcess Files',b{1}));
geox = data(:,2);
geoy = data(:,3);

% find max venus over time
maxvx = data(:,8);
maxvy = data(:,9);

% find orientation of blob
u = strfind(exp,'_');
short = exp(1:u(end-1)-1);
a = file_search(['blb_anl_rp_' short '_FRET.txt'],folder);
blb = load(fullfile(folder,a{1}));
blbcol = 27; % SHOULDN'T HARD CODE THESE...
imgcol = 28;
cellcol = 29;
orcol = 26;
imgn = str2num(exp(u(end)+1:end));
% imgn = str2num(exp(u(end)+1:end))-1;
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

h1=figure; set(h1,'Visible','off')
plot([0:5:325],rmaxx,'k','LineWidth',3); axis([0 325 -50 50])
saveas(h1,fullfile(folder,'FRAP Venus Max Displacement',['Maj Axis Venus Disp_' exp '.png']))
h2=figure; set(h2,'Visible','off')
plot([0:5:325],rmaxy,'k','LineWidth',3); axis([0 325 -50 50])
saveas(h2,fullfile(folder,'FRAP Venus Max Displacement',['Min Axis Venus Disp_' exp '.png']))

rmaxcoords = [rmaxx rmaxy];

save(fullfile(folder,'FRAP Venus Max Displacement',['Venus Max Coords_' exp '.txt']),'rmaxcoords','-ascii')

end