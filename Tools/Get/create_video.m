function create_video(mainDirName,imageDir,videoFileName)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

videoDir ='Videos';  
mkdir(mainDirName,videoDir);
outputVideo = VideoWriter(fullfile([mainDirName,'\',videoDir,...
    '\',videoFileName]));
outputVideo.FrameRate = 8;
open(outputVideo);
imageFileStr = get_files_in_folder(strcat(mainDirName,'\',imageDir));
hWait = waitbar(0);
nn = length(imageFileStr);
for ii = 1:nn
   custom_waitbar(hWait,ii,nn,'Generating Video');
   img = imread(fullfile(mainDirName,imageDir,imageFileStr{ii}));
   writeVideo(outputVideo,img)
end
close(outputVideo);
delete(hWait);

end

