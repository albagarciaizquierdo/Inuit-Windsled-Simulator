function fun_download_fig(foldername,figures,width,height)
%% Description:
% This function downloads figures in .jpg format to an specific folder with a 
% resolution of 600 ppi and specific width and height
%% Inputs:
% foldername --> string of the desired folder name 
% figures --> figures as a variable coming from: figures = findall(0, 'Type', 'figure');
% width --> width of the desired image
% height --> height of the desired image
%% Ouputs:
% Downloads the image
%%
% Folder to save figure
if ~exist(foldername, 'dir')
    mkdir(foldername);
end

% Save figure
for i = 1:length(figures)
    figure(figures(i));
    filename = fullfile(foldername, sprintf('figure_%d.jpg', i));
    set(figures(i), 'PaperPositionMode', 'auto');
    set(figures(i), 'Position', [0, 0, width, height]);
    print(figures(i), filename, '-djpeg', '-r600'); % Resolution at 600 ppi
end
end