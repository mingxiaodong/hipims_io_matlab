function WriteWinscpCMD(filename,loacalDir,remoteDir,mainStatements)
%%
%********************write file******************
fileID = fopen(filename,'w');
cmdFileHeadstr = {'cd C:\Program Files (x86)\WinSCP','@echo off'...
    'winscp.com /command ^',...
    '  "open mysession" ^',...
    ['  "cd ' remoteDir ' " ^'],...
    ['  "lcd ' loacalDir ' " ^']};
for i=1:length(cmdFileHeadstr)
    fprintf(fileID,'%s\n',cmdFileHeadstr{i});
end
for i=1:numel(mainStatements)
    fprintf(fileID,'  %s\n',['"' mainStatements{i} '"^']);
end
fprintf(fileID,'  %s\n','"exit"');
fclose(fileID);
end
%*************************************************