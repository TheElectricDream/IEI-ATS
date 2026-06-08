function exportTightRasterPDF(fig, filename)
    % EXPORTTIGHTRASTERPDF Exports a tightly cropped rasterized plot to PDF.
    %
    % Inputs:
    %   fig      - Handle to the figure (e.g., gcf)
    %   filename - String or char array of the output filename (must end in .pdf)
    
    % 1. Force the backgrounds to white
    set(fig, 'Color', 'w');
    axList = findall(fig, 'Type', 'axes');
    set(axList, 'Color', 'none'); 
    
    % 2. Render directly to a high-res RGB matrix in memory
    img = print(fig, '-RGBImage', '-r300');
    
    % 3. Create the background mask (white = 255 across all channels)
    bgMask = all(img == 255, 3);
    [rows, cols] = find(~bgMask);
    
    % Safety check
    if isempty(rows) || isempty(cols)
        warning('Figure appears completely white. Aborting export.');
        return;
    end
    
    % 4. Crop the matrix to the exact mathematical bounds of the data
    croppedImg = img(min(rows):max(rows), min(cols):max(cols), :);
    
    % 5. Get the dimensions of the cropped image matrix
    [imgH, imgW, ~] = size(croppedImg);
    
    % 6. Create the temporary, invisible figure
    tempFig = figure('Visible', 'off');
    
    % --- THE FIX ---
    % Resize the figure window to exactly match the image's aspect ratio.
    % This stops MATLAB from adding "letterboxing" margins to fit a default window.
    baseWidth = 800; 
    scaledHeight = round(baseWidth * (imgH / imgW));
    tempFig.Position(3:4) = [baseWidth, scaledHeight];
    
    % Create an axes that takes up the perfectly proportioned figure
    tempAx = axes(tempFig, 'Position', [0 0 1 1]);
    
    % Display the perfectly cropped image
    imshow(croppedImg, 'Parent', tempAx);
    
    % 7. Export as a rasterized PDF
    % Added the 300 Resolution flag so exportgraphics doesn't downsample the matrix
    exportgraphics(tempAx, filename, 'ContentType', 'image', 'Resolution', 300);
    
    % 8. Clean up the temporary figure from memory
    close(tempFig);
end