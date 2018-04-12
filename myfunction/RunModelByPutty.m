function RunModelByPutty(remoteDir)
    currentFolder = pwd;
    fileID = fopen([currentFolder '\\PuttyRemoteScript.txt'],'w');
    fprintf(fileID,['cd ' remoteDir '\n']);
    fprintf(fileID,'../MultiGPUsUrbanFloodSimulator');
    fclose(fileID);
    %require a file
    fileID = fopen([currentFolder '\\RunPuttyModel.cmd'],'w');
    fprintf(fileID,'C: \n');
    fprintf(fileID,'cd "C:\\Program Files (x86)\\PuTTY" \n');
    currentFolderC = strrep(currentFolder,'\','\\');
    fprintf(fileID,...
        ['plink.exe -ssh ceg-gpu01.ncl.ac.uk -l b1055010 -pw xiaxilingpu -m ',...
        currentFolderC '\\PuttyRemoteScript.txt']);
    fclose(fileID);
    command = [currentFolder '\\RunPuttyModel.cmd'];
    disp(['Model at ' remoteDir ' is running...'])
    system(command,'-echo');
    cd(currentFolder)
end