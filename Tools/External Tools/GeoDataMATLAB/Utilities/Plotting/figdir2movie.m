function figdir2movie(fig_dir,outname,varargin)
% figdri2movie
% figdir2movie(fig_dir,outname,s)
% This will take a directory with a set of figures and make them into a
% movie.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs
% fig_dir - A string that is the location of the directory or a string that
% the dir command can find fig files.
% outname - The name of the video
% s - A struct that is set up with the same variable names as the
% videoWriter class
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if strcmp(fig_dir(end-3:end),'.fig')
    fig_info = dir(fig_dir);
    fig_str = fig_dir;
    [fig_dir,~,~] = fileparts(fig_dir);
else
    fig_info = dir(fullfile(fig_dir,'*.fig'));
end
writerObj = VideoWriter(outname);
% populate the writerObj object
if nargin >2
    s = varargin{1};
    names = fieldnames(s);
    for k = 1:length(names)
        set(writerObj,names{k},getfield(s,names{k}));
    end
end
% start saving things
open(writerObj);
% Make all the positions the same
h = open(fullfile(fig_dir,fig_info(1).name));
pause(1)
%posvec = get(h,'Position');
posvec = h.Position;
close(h)
for k = 1:length(fig_info)
    
    h = open(fullfile(fig_dir,fig_info(k).name)); 
    set(h,'Position',posvec);
    pause(1)
    frame = getframe(h);
    writeVideo(writerObj,frame);
    close(h)
end
close(writerObj);
