function create_video(mainDirName,imageDir,videoFileName)
% create_video Creates .avi file with images in imageDir
%   Detailed explanation goes here
multiWaitbar('Generating Video',0);
videoDir ='Videos';
mkdir(mainDirName,videoDir);
outputVideo = VideoWriter(fullfile([mainDirName,filesep,videoDir,...
    filesep,videoFileName]));
outputVideo.FrameRate = 8;
open(outputVideo);
imageFileStr = get_files_in_folder(strcat(mainDirName,filesep,imageDir));

nn = length(imageFileStr);
di = 1./nn;
for ii = 1:nn
   multiWaitbar('Generating Video','Increment',di);
   img = imread(fullfile(mainDirName,imageDir,imageFileStr{ii}));
   writeVideo(outputVideo,img);
end
close(outputVideo);
multiWaitbar('Generating Video','Close');

end
