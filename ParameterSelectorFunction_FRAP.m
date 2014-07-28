function [ParameterValues] = ParameterSelectorFunction(Image,WidthRange,ThreshRange,MergeRange,ParameterValues)

% A program that allows the user to manually optimize and select fa_gen
% parameters.

    % Determine colors
    numColors = 16;
    Colors = distinguishable_colors(numColors,'k');
    
    % Change contrast
    topval = prctile(Image(:),99.9);
    ImagePlot = mat2gray(Image, [0 topval]);
    
    % Determine Focal Adhesions
    FA = water(Image,ParameterValues);
    FAone = FA ~= 0;
    FAoutline = bwperim(FAone);
    FAcolor = FAoutline .* FA;

    % Change contrast
    topval = prctile(Image(:),99.9);
    ImagePlot = mat2gray(Image, [0 topval]);
    
    % Focal adhesion figure
    hFig = figure;
    hAx = axes('Parent',hFig);

    imagesc(ImagePlot,'Parent',hAx)
    colormap(gray)
    
    % Show FA outlines
    hold on
    for i = 1:numColors
        colorMask{i} = cat(3, Colors(i,1)*ones(size(FA)), Colors(i,2)*ones(size(FA)), Colors(i,3)*ones(size(FA)));
        h{i} = imshow(colorMask{i},'Parent',hAx);
        x = rem(FAcolor,numColors);
        set(h{i},'AlphaData',rem(FAcolor,numColors)+1==i & FAcolor~=0);
    end
    hold off

    % Slider 1: Width
    uicontrol('Parent',hFig,'Style','slider','Min',WidthRange(1),'Max',WidthRange(2),'Value',ParameterValues(1),'Position',[0 350 120 20],'SliderStep',[.05 .05],'Callback',@ChangeFA1); 
    hText1 = uicontrol('Style','text','Position',[0 370 120 20],'String',sprintf('Width: %d',ParameterValues(1)));

    % Slider 2: Threshold
    uicontrol('Parent',hFig,'Style','slider','Min',ThreshRange(1),'Max',ThreshRange(2),'Value',ParameterValues(2),'Position',[0 290 120 20],'SliderStep',[.002 .002],'Callback',@ChangeFA2); 
    hText2 = uicontrol('Style','text','Position',[0 310 120 20],'String',sprintf('Thresh: %d',ParameterValues(2)));

    % Slider 3: Merge
    uicontrol('Parent',hFig,'Style','slider','Min',MergeRange(1),'Max',MergeRange(2),'Value',ParameterValues(3),'Position',[0 230 120 20],'SliderStep',[.125 .125],'Callback',@ChangeFA3); 
    hText3 = uicontrol('Style','text','Position',[0 250 120 20],'String',sprintf('Merge: %d',ParameterValues(3)));
    
    w = waitforbuttonpress;
    close all
    
    function ChangeFA1(hObj,event)
        % Change width value
        val = get(hObj,'value');
        ParameterValues(1) = val;
        
        % Re-determine focal adhesions
        FA = water(Image,ParameterValues);
        FAone = FA ~= 0;
        FAoutline = bwperim(FAone);
        FAcolor = FAoutline .* FA;
        
        % Re-plot focal adhesions
        imagesc(ImagePlot,'Parent',hAx)
        colormap(gray)
        
        % Show FA outlines
        hold on
        for i = 1:numColors
            colorMask{i} = cat(3, Colors(i,1)*ones(size(FA)), Colors(i,2)*ones(size(FA)), Colors(i,3)*ones(size(FA)));
            h{i} = imshow(colorMask{i},'Parent',hAx);
            x = rem(FAcolor,numColors);
            set(h{i},'AlphaData',rem(FAcolor,numColors)+1==i & FAcolor~=0);
        end
        hold off
        
        % Change text value
        set(hText1,'String',sprintf('Width: %d',ParameterValues(1)))
    end

    function ChangeFA2(hObj,event)
        % Change width value
        val = get(hObj,'value');
        ParameterValues(2) = val;
        
        % Re-determine focal adhesions
        FA = water(Image,ParameterValues);
        FAone = FA ~= 0;
        FAoutline = bwperim(FAone);
        FAcolor = FAoutline .* FA;
        
        % Re-plot focal adhesions
        imagesc(ImagePlot,'Parent',hAx)
        colormap(gray)
        
        % Show FA outlines
        hold on
        for i = 1:numColors
            colorMask{i} = cat(3, Colors(i,1)*ones(size(FA)), Colors(i,2)*ones(size(FA)), Colors(i,3)*ones(size(FA)));
            h{i} = imshow(colorMask{i},'Parent',hAx);
            x = rem(FAcolor,numColors);
            set(h{i},'AlphaData',rem(FAcolor,numColors)+1==i & FAcolor~=0);
        end
        hold off
        
        % Change text value
        set(hText2,'String',sprintf('Thresh: %d',ParameterValues(2)))
    end

    function ChangeFA3(hObj,event)
        % Change width value
        val = get(hObj,'value');
        ParameterValues(3) = val;
        
        % Re-determine focal adhesions
        FA = water(Image,ParameterValues);
        FAone = FA ~= 0;
        FAoutline = bwperim(FAone);
        FAcolor = FAoutline .* FA;
        
        % Re-plot focal adhesions
        imagesc(ImagePlot,'Parent',hAx)
        colormap(gray)
        
        % Show FA outlines
        hold on
        for i = 1:numColors
            colorMask{i} = cat(3, Colors(i,1)*ones(size(FA)), Colors(i,2)*ones(size(FA)), Colors(i,3)*ones(size(FA)));
            h{i} = imshow(colorMask{i},'Parent',hAx);
            x = rem(FAcolor,numColors);
            set(h{i},'AlphaData',rem(FAcolor,numColors)+1==i & FAcolor~=0);
        end
        hold off
        
        % Change text value
        set(hText3,'String',sprintf('Merge: %d',ParameterValues(3)))
    end

end