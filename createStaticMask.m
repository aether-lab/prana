function MASK = createStaticMask(IMAGE, WRITEPATH)

 MASK = ones(size(IMAGE,1),size(IMAGE,2));
    roiwindow = CROIEditor(IMAGE);
    while isempty(roiwindow.labels)
        addlistener(roiwindow,'MaskDefined',@your_roi_defined_callback);
        drawnow
    end
    
    function your_roi_defined_callback(h,e)
        [MASK, labels, n] = roiwindow.getROIData;


    end

% Set masked portion to zero
MASK(roiwindow.labels == 0) = 0;
delete(roiwindow);

% Write the mask if a file path is specified
if nargin > 1
    imwrite(MASK, WRITEPATH);
end

end



   