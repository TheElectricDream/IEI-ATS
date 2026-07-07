function [h_figs, h_axs, h_imgs, video_writers] = ...
    initializeEventVideos(coh_out, ats_out, image_size, video_out_path)

    %   This function initializes off-screen figures and AVI
    %   video writers for recording the LATS surface and coherence map
    %   outputs during processing.
    %
    %   Inputs:
    %     coh_out - Logical. If true, initialize coherence video output.
    %     ats_out - Logical. If true, initialize ATS video output.
    %     image_size  - [1 x 2] Image dimensions [nRows, nCols].
    %     video_out_path - String. Where the output video should be saved.
    %
    %   Outputs:
    %     h_figs        - {1 x 2} Cell array of figure handles.
    %     h_axs         - {1 x 2} Cell array of axes handles.
    %     h_imgs        - {1 x 2} Cell array of image handles.
    %     video_writers - {1 x N} Cell array of VideoWriter objects.

    % Initialize output cells
    h_figs  = {[], []};
    h_axs   = {[], []};
    h_imgs  = {[], []};
    video_writers = {};

    % Create the LATS output video
    if ats_out
        hFigATS = figure('Visible', 'off', ...
            'Position', [100 100 image_size]);
        hAxATS  = axes('Parent', hFigATS);
        colormap(hAxATS, 'gray');
        colorbar(hAxATS);
        set(hAxATS, 'FontSize', 16, 'Color', 'white', ...
            'XTick', [], 'YTick', []);
        initial_data = nan(image_size(2), image_size(1));
        hImgATS = imagesc(hAxATS, initial_data);
        set(hImgATS, 'AlphaData', ~isnan(initial_data));

        videoFileName = fullfile(video_out_path);
        video_writers{1} = VideoWriter(videoFileName);
        video_writers{1}.FrameRate = 60;
        open(video_writers{1});

        h_figs{1} = hFigATS;
        h_axs{1}  = hAxATS;
        h_imgs{1} = hImgATS;
    end

    % Create the coherence output video
    if coh_out
        hFigCOH = figure('Visible', 'off', ...
            'Position', [100 100 image_size]);
        hAxCOH  = axes('Parent', hFigCOH);
        colormap(hAxCOH, 'gray');
        colorbar(hAxCOH);
        set(hAxCOH, 'FontSize', 16, 'Color', 'white');
        initial_data = nan(image_size(2), image_size(1));
        hImgCOH = imagesc(hAxCOH, initial_data);
        set(hImgCOH, 'AlphaData', ~isnan(initial_data));

        videoFileName = fullfile(pwd, ...
            'coherence_map_output.avi');
        video_writers{2} = VideoWriter(videoFileName);
        video_writers{2}.FrameRate = 60;
        open(video_writers{2});

        h_figs{2} = hFigCOH;
        h_axs{2}  = hAxCOH;
        h_imgs{2} = hImgCOH;
    end

end