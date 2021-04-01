function [t,Ca,Ct] = import_voistat(path,folder)

file_IF = 'IF.voistat';
file_TAC = 'tumor_tac.voistat';

delimiter = '\t';
startRow = 4;
formatSpec = '%s%s%f%f%f%f%f%f%f%f%f%f%f%f%f%s%*s%*s%*s%[^\n\r]';

%% import IF
filename_IF = strcat(path,folder,file_IF);

fileID = fopen(filename_IF,'r');
dataArray_IF = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'HeaderLines' ,startRow-1, 'ReturnOnError', false);
fclose(fileID);

t = dataArray_IF{:,4}./60;
Ca = dataArray_IF{:,7};

%% import tumor TAC

filename_TAC = strcat(path,folder,file_TAC);

fileID = fopen(filename_TAC,'r');
dataArray_TAC = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'HeaderLines' ,startRow-1, 'ReturnOnError', false);
fclose(fileID);

Ct = dataArray_TAC{:,7};

end




