folder = input('Enter folder containing images: ','s');
blob_images = file_search('fa_\w+.TIF',folder);
exp_name = input('Enter the experimental group: ','s');
totdata = [];
for i = 1:length(blob_images)
    outname = ['blbdecay_' blob_images{i}(8:end-4) '.txt'];
    blb_file = (file_search(['blbd_' exp_name '_' num2str(i) '.txt'],folder));
    if isempty(blb_file)
%         blob_decay([exp_name '_FRAP_' num2str(i,'%02d') '_t\d+.TIF'],blob_images{i},outname,folder)
        blob_decay([exp_name '_s' num2str(i) '_t\d+.TIF'],blob_images{i},outname,folder);
    else
        blob_decay([exp_name '_FRAP_' num2str(i,'%02d') '_t\d+.TIF'],blob_images{i},outname,blb_file{1},folder)
    end
    
    addpath(folder)
%     figure; hold on; box on;
    data = load(outname);
    data(1,:) = [];
    data(:,1) = [];
%     slopes = zeros(1,size(data,2)-1);
    for j = 1:size(data,2)
        data(:,j) = data(:,j)./max(data(:,j));
%         plot(data(5:end,1),data(5:end,j)./max(data(5:end,j)))
%         pfit = polyfit(data(5:end,1),data(5:end,j)./max(data(5:end,j)),1);
%         slopes(j-1) = pfit(1);
    end
    totdata = [totdata data];
%     title([exp_name ' ' num2str(i)])
%     xlabel('Time')
%     ylabel('Intensity')
    rmpath(folder)
%     pos = find(slopes > (mean(slopes)+abs(.1*mean(slopes))));
%     neg = find(slopes < (mean(slopes)-abs(.1*mean(slopes))));
%     disp([exp_name ' ' num2str(i)])
%     fprintf('The percent of FAs with increasing intensity is %.2f\n',length(pos)/length(slopes))
%     fprintf('The percent of FAs with decreasing intensity is %.2f\n',length(neg)/length(slopes))
end

meandata = mean(totdata,2);
stddata = std(totdata,0,2);
stedata = stddata./sqrt(length(blob_images));