function plot_cents(exp,folder)

% find polygons
a = file_search(['poly_FRAP_con_' exp '_t02_0.dat'],folder);
p = load(fullfile(folder,'FRAP Poly Files',a{1}));
a2 = file_search(['poly_FRAP_con_' exp '_t66_0.dat'],folder);
p2 = load(fullfile(folder,'FRAP Poly Files',a2{1}));

% find geometric centroid @ initial point
b = file_search(['con_' exp '_poly_data.txt'],folder);
data = load(fullfile(folder,'FRAP PostProcess Files',b{1}));
geox = data(2,2);
geoy = data(2,3);
geox2 = data(66,2);
geoy2 = data(66,3);

% find max venus over time
maxvx = data(:,8);
maxvy = data(:,9);

% plot everything
h1 = figure; hold on; set(h1,'Visible','off')
patch(p(:,1),p(:,2),'w','EdgeColor','r','FaceColor','none','LineWidth',3)
patch(p2(:,1),p2(:,2),'w','EdgeColor','g','FaceColor','none','LineWidth',3)

plot(maxvy,maxvx,'b','LineWidth',2)
plot(geoy,geox,'kx','LineWidth',3)
plot(geoy2,geox2,'kx','LineWidth',3)
axis([round(min([p(:,1); p2(:,1)])-5) round(max([p(:,1); p2(:,1)])+5) round(min([p(:,2); p2(:,2)])-5) round(max([p(:,2); p2(:,2)])+5)])
saveas(h1,fullfile(folder,'FRAP Venus Max Displacement',['Venus Max_' exp '.png']))

end