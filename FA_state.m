function FA_state(exp_name,folder)

% Find all poly_FRAP_blch that are completely analyzed
a = file_search(['poly_FRAP_blch_' exp_name '_\d+_t66_0.dat'],folder);
% a = file_search(['poly_FRAP_blch_' exp_name '_\d+_t66_66_2.dat'],folder);
r = 2048;
c = 2048;
imgix = repmat(1:r,c,1);
imgiy = repmat([1:c]',1,r); %#ok<NBRAK>

for i = 1:length(a)
   % find centroids
   p1 = load(strrep(a{i},'t66','t02'));
%    temp = strrep(a{i},'t66','t02');
%    temp = strrep(temp,'66','2');
%    p1 = load(temp);
   p2 = load(a{i});
   in1 = inpolygon(imgix,imgiy,p1(:,1),p1(:,2));
   in2 = inpolygon(imgix,imgiy,p2(:,1),p2(:,2));
   [y1,x1] = ind2sub([r c],find(in1==1));
   [y2,x2] = ind2sub([r c],find(in2==1));
   cx1 = mean(x1);
   cy1 = mean(y1);
   cx2 = mean(x2);
   cy2 = mean(y2);
   
   figure; hold on;
   patch(p1(:,1),p1(:,2),'w','EdgeColor','r','FaceColor','none')
   patch(p2(:,1),p2(:,2),'w','EdgeColor','g','FaceColor','none')
   plot(cx1,cy1,'r*')
   plot(cx2,cy2,'g*')
   plot([cx1 cx2],[cy1,cy2],'k')
   axis([min([p1(:,1);p2(:,1)])-5,max([p1(:,1);p2(:,1)])+5,min([p1(:,2);p2(:,2)])-5,max([p1(:,2);p2(:,2)])+5]) 
   a1 = length(x1);
   a2 = length(x2);
   perarea = abs(a2-a1)/a1 * 100;
   title(strrep(a{i},'_',' '))
   
   % find major/minor axis
   xsh1 = x1-cx1;
   ysh1 = y1-cy1;
   uxx1 = sum(xsh1.^2)/length(x1);
   uyy1 = sum(ysh1.^2)/length(y1);
   uxy1 = sum(xsh1.*ysh1)/length(x1);
   qrot1 = sqrt((uxx1-uyy1).^2+4*(uxy1.^2));
   mjra1 = sqrt(2)*sqrt(uxx1+uyy1+qrot1);
   mnra1 = sqrt(2)*sqrt(uxx1+uyy1-qrot1); %#ok<*NASGU>
   xsh2 = x2-cx2;
   ysh2 = y2-cy2;
   uxx2 = sum(xsh2.^2)/length(x2);
   uyy2 = sum(ysh2.^2)/length(y2);
   uxy2 = sum(xsh2.*ysh2)/length(x2);
   qrot2 = sqrt((uxx2-uyy2).^2+4*(uxy2.^2));
   mjra2 = sqrt(2)*sqrt(uxx2+uyy2+qrot2);
   mnra2 = sqrt(2)*sqrt(uxx2+uyy2-qrot2);
   
   % find displacement of centroid
   cntdist = sqrt((cx1-cx2)^2+(cy1-cy2)^2);
   
   % find overlapping area
   pnts1 = find(in1==1);
   pnts2 = find(in2==1);
   overlap = intersect(pnts1,pnts2);
   perover = length(overlap)/a1*100;
   
%    classify - stable, growing, sliding
   state = [];
   if perarea < 30 && perover > 60
       state = 'stable';
   elseif perarea > 30 && perover > 60
       state = 'growing';
   elseif perover < 60
       state = 'sliding'; 
   end
   
   if a1 > a2
       state = [state '-SHRINKING']; %#ok<AGROW>
   end
   
   legend(['Area 1 = ' num2str(a1) '; Major Axis = ' num2str(mjra1)],['Area 2 = ' num2str(a2) '; Major Axis = ' num2str(mjra2)],...
       ['Centroid Displaced = ' num2str(cntdist)],['Percent Overlap = ' num2str(perover) '%'],...
       ['State = ' state])   
   
end

end