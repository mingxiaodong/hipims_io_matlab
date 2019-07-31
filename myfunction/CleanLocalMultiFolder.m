function CleanLocalMultiFolder(localMultiDir)
% CleanLocalMultiFolder(localMultiDir,numGPU) clean the all output folder of multi-GPU results
for i=0:7
    dirName = [localMultiDir '/' num2str(i) '/output'];
    if exist(dirName,'dir')==7
        delete([dirName '/*'])
    end
end