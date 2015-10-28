function plot_cents(exp,folder)

% find polygon
a = file_search(['poly_FRAP_con_' exp '_t02_0.dat'],folder);
p = load(a{1});

% find geometric centroid @ initial point
b = file_search(['con_' exp '_poly_data.txt'],folder);
data = load(b{1});
geox = data(2,2);
geoy = data(2,3);

% find max venus over time
maxvx = data(:,8);
maxvy = data(:,9);

% plot everything
figure; hold on;
patch(p(:,1),p(:,2),'w','EdgeColor','g','FaceColor','none','LineWidth',3)
plot(maxvy,maxvx,'b','LineWidth',2)
plot(geoy,geox,'kx','LineWidth',3)
axis([round(min(p(:,1))-5) round(max(p(:,1))+5) round(min(p(:,2))-5) round(max(p(:,2))+5)])

end