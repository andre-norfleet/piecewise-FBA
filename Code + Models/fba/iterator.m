folder_path = '/storage/scratch1/5/dnorfleet7/my_simdre2/piecewise-FBA/Code + Models/fba/input2'; % Replace with the actual folder path
folder_contents = dir(folder_path);
dir_flags = [folder_contents.isdir];
subdirectories = folder_contents(dir_flags);
subdirectory_names = {subdirectories(3:end).name}; % Exclude '.' and '..'
ghh = string(subdirectory_names);

% Display the list of subdirectories
disp(subdirectory_names)
gg = {};
V = 22:42;
ahh = compose("/storage/scratch1/5/dnorfleet7/my_simdre2/piecewise-FBA/Code + Models/fba/runfba%d.m",V(:));
filename = '/storage/scratch1/5/dnorfleet7/my_simdre2/piecewise-FBA/Code + Models/fba/runfba.m';
% Read the file contents into a cell array, one line per cell
fileLines = readlines(filename);

% Replace line 7 with '{new input folder name}'
if numel(fileLines) >= 3
    fileLines(7) = 'input_folder = "test1_25"';
else
    % If file has fewer than 3 lines, pad with empty lines
    fileLines(numel(fileLines)+1:3) = "";
    fileLines(7) = "ln = 33";
end
%writelines(fileLines, ahh(1,1));




for i = 1:length(subdirectory_names)
    fileLines(7,1) = strrep(fileLines(7,1), 'test1_25', ghh(1,i));
    writelines(fileLines, ahh(i,1));
%     filename = 'runfba.m';
%     fileContents = fileread(filename);
%     lines = splitlines(fileContents);
%     lines = strrep(lines{2}, 'test1_25', 'gg');
%     newContents = strjoin(lines, newline);
%     fid = fopen(filename, 'w');
%     fwrite(fid, lines);
%     fclose(fid);
    %writelines(lines, ahh(1,i));
end