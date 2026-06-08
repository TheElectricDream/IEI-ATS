function [] = showImageAsVideo(image)
    % Plot the global hot mask so changes can be watched each iteration
    hFig = figure(20);
    try
        imagesc(image);
        colormap(gca, hot);
        colorbar;
        axis image off;
        drawnow;
    catch
        % If plotting fails, ensure the figure handle is closed to avoid clutter
        if ishghandle(hFig)
            close(hFig);
        end
    end
end
