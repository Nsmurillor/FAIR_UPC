figFiles = dir(fullfile(info.plot.vid_dir, '*.png'));
fig_names = {figFiles.name};

frame_files = {};

for k = 1:length(fig_names)
    if startsWith(fig_names{k}, 'Fig1_')
        frame_files{end+1} = fig_names{k};
    end
end

nums = zeros(size(frame_files));
for k = 1:length(frame_files)
    tok = regexp(frame_files{k}, '_\d+', 'match');  
    nums(k) = str2double(tok{1}(2:end));           
end

[~, idx] = sort(nums);
files_sorted = frame_files((idx));

vid = VideoWriter([folder_results,'Animation.mp4'],'MPEG-4');
vid.FrameRate = 30;   
open(vid);

img0 = imread([info.plot.vid_dir,files_sorted{1}]);
targetSize = [size(img0,1), size(img0,2)];
for ii_fig=1:size(files_sorted,2)
    
        frame = imread([info.plot.vid_dir,files_sorted{ii_fig}]);
        bool_resize=size(frame,1) ~= targetSize(1) || size(frame,2) ~= targetSize(2);
        if bool_resize
            frame = imresize(frame, targetSize);
        end
        for kk=1:3
            writeVideo(vid, frame);   
        end
      
    
end

rmdir(info.plot.vid_dir,'s');
close(vid)

