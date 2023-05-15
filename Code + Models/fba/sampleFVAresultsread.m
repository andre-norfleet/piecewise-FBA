filename = 'Z:/Andre/From Andre/FBA-pipeline-master-Riya/Code + Models/fba/results/WTC11_D0_Zhangetal_ATPallmTeSRplusAlbuMAX_relaxconstralltest_4/_TEMPLATE_1/fva/WTC_-dox_avg.tsv';
delimiter = '\t';
startRow = 1;

formatSpec = '%s%f%[^\n\r]';
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'TextType', 'string', 'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
fclose(fileID);

columnB = dataArray{:, 1};
columnC = dataArray{:, 2};

disp(columnB)
disp(columnC)

