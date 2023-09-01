% source path of the video
vid_path = 'C:\source\free-viewing\videos\leftfront1.mp4';
% destination folder in which to save the clip
output_path = 'C:\source\free-viewing\videos\convert';

% clip timestamps
start_minute = 8;
start_second = 34;
clip_dur_s = 6.5;

start_s = start_minute * 60 + start_second;
start_end_times = [ start_s, start_s + clip_dur_s ];

make_video_clips( vid_path, output_path, start_end_times, [] );