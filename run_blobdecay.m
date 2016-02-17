folder = input('Enter folder containing images: ','s');
blob_images = file_search('fa_\w+.TIF',folder);
exp_name = input('Enter the experimental group: ','s');
for i = 1:length(blob_images)
    outname = ['blbdecay_' blob_images{i}(8:end-4) '.txt'];
    blb_file = (file_search(['blbd_' exp_name '_' num2str(i) '.txt'],folder));
    if isempty(blb_file)
        blob_decay([exp_name '_FRAP_' num2str(i,'%02d') '_t\d+.TIF'],blob_images{i},outname,folder)
    else
        blob_decay([exp_name '_FRAP_' num2str(i,'%02d') '_t\d+.TIF'],blob_images{i},outname,blb_file{1},folder)
    end
    
%     addpath(folder)
%     figure; hold on; box on;
%     data = load(outname);
%     data(1,:) = [];
%     for j = 2:size(data,2)
%         plot(data(5:end,1),data(5:end,j))
%     end
%     title([exp_name ' ' num2str(i)])
%     xlabel('Time')
%     ylabel('Intensity')
%     rmpath(folder)
end