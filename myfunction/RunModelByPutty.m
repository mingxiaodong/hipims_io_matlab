function RunModelByPutty(remoteDir)
    currentFolder = pwd;
    %create a txt script to run model via putty
    fileID = fopen([currentFolder '\\PuttyRemoteScript.txt'],'w');
    fprintf(fileID,['cd ' remoteDir '\n']);
    % model name and remote path
    fprintf(fileID,'../MultiGPUsUrbanFloodSimulator');
    fclose(fileID);
    % create a cmd file to execute PuttyRemoteScript.txt 
    fileID = fopen([currentFolder '\\RunPuttyModel.cmd'],'w');
    fprintf(fileID,'C: \n');
    % local path of putty
    fprintf(fileID,'cd "C:\\Program Files (x86)\\PuTTY" \n');
    currentFolderC = strrep(currentFolder,'\','\\');
    % IP, username and password of remote server
    fprintf(fileID,...
        ['plink.exe -ssh ceg-gpu01.ncl.ac.uk -l b1055010 -pw xiaxilingpu -m ',...
        currentFolderC '\\PuttyRemoteScript.txt']);
    fclose(fileID);
    command = [currentFolder '\\RunPuttyModel.cmd'];
    disp(['Model at ' remoteDir ' is running...'])
    system(command,'-echo');
    cd(currentFolder)
end