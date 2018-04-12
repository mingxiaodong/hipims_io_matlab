function [ status ] = NIMROD_wintest1( dirname )
%NIMROD_WINTEST1 Demo function for decompressing, reading and displaying a 
%set of BADC - UK MetOffice NIMROD composite uk-1km radar data files (V1.7)
% This program is suitable only for running under MS Windows versions of
% MATLAB; However, adaptation to unix platforms should be relatively
% simple.
% TYPICAL MATLAB COMMAND TO CALL THIS FUNCTION:
% status = NIMROD_wintest1( 'NIMRODdata\2007-08-19to2007-08-21\20' );
%
% ASSUMPTIONS:  
%   1) Files have been downloaded from BADC website in *.gz compressed 
%   format; 
%   2) Program 7-zip is used to unzip each *.gz compressed file.
%   NB: If you do not have 7-zip (7z.exe) installed on your PC, you will
%   need to modify the 'lcommand' and 'uzcommand' entries below to use a
%   suitable alternative file decompression program, or download and
%   install 7-zip in C:\Program Files\7-Zip directory.
%   3) This program has been tested on MATLAB V 2010a (64-bit) -it may work
%   on other versions!
%   4) Directory 'dirname' contains the NIMROD uk-1km composite data files  
%   of interest and is on the current MATLAB path (I used a subdirectory
%   of the current directory where both *.m files were. 'dirname' then just
%   contained a relative path from the current directory).
%   5) For NIMROD composite uk-1km data files, there is just one file in 
%   each *.gz archive file. This contains a snapshot of the entire UK 
%   est Total Precipitation Rate, 1 pixel per 1km square, corresponding 
%   to Ordnance Survey National Grid squares, for a period of 5-minutes. 
%   There is a header of 523 bytes, followed by a rectangular integer array
%   of estimated Precipitation Rates; a value of 1 in an array location 
%   corresponds with 1/32 of a mm of rain per hour (valid for the 5-minute 
%   reference period).
%
% DETAILED DESCRIPTION:
%   A listing of the *.gz files in directory passed in 'dirname' is taken;
%   A MATLAB figure is created for the resultant radar data plots;
%   Each compressed *.gz file is unzipped in turn (into the current 
%   directory) prior to reading in the unzipped file data; 
%   Each unzipped file is deleted to save space;
%   Each radar data 5-min image is displayed in the above figure, using 
%   axes corresponding to the Ordnance Survey National Grid. Note that the
%   y-data have to be plotted in reverse order to ensure image is displayed
%   correct way up and y-axis scale is displayed correctly using axis xy.
%   A yellow square (20x20 pixels) is included in the image in the top 
%   right hand corner, to help ensure it is being displayed the right way 
%   round. NB: Remove the code for this, when no longer needed.
%   The clims parameter used for the imagesc command limits the colours
%   diplayed in the plot to 0 (no rain) -(although some cells will actually
%   contain negative values) <= 1000, which corresponds with intensity of
%   1000/32 = 31.25mm/hr. The maximum possible value from the data array is
%   +32767/32 = 1024mm/hr, so clims can be changed if studying very intense
%   rainstorms, but less intense ones will not be visible in this case.
%   It is possible to overlay an outline of Great Britain over the top of
%   the display by using 'plot' statements to plot a vectorised outline.
%   The outline data are available for download from OS Edina website:
%   http://edina.ac.uk/ukborders/
%   They have not been included here due to copyright considerations.
%   An example location is marked with a white 'x', for illustration.
%   For this, OS 6-digit easting and 6-digit northing national grid values 
%   should be used. More information can be found at: 
%   http://en.wikipedia.org/wiki/British_national_grid_reference_system
%

% National Grid co-ordinates of example location of interest 
% (define as required)
loc_easting = 406790;
loc_northing = 442082;

if (~isdir(dirname))
    error('Specified directory not found in NIMROD_wintest1');
end

% read compressed NIMROD filenames from directory
F = dir(sprintf('%s%s*.gz', dirname, filesep));
No_files = length(F);

% create MATLAB figure to fill left-hand half of screen (minus taskbar at
% bottom)
scrsz = get(0,'ScreenSize');
fh1 = figure('OuterPosition',[1 scrsz(4)*.05 scrsz(3)*0.5 scrsz(4)* 0.95]);

% uncompress each .gz file using 7-zip, read it and display it on National
% Grid co-ordinate axes
for f = 1:No_files
    Fname = F(f).name;
    % list contents of this .gz file (no paths,assume yes to all questions)
    lcommand =['C:\PROGRA~1\7-Zip\7z.exe l -y ' dirname filesep Fname ];
    [status,result] = system(lcommand,'-echo'); 
    if (status~=0)
        error('*.gz file contents not listed correctly in NIMROD_wintest1')
    end
    % retrieve name of the single file it contains
    Zfname = strtrim(result(347:414));
    % uncompress contents of the .gz file 
    % (no paths, assume yes to all questions)
    uzcommand =['C:\PROGRA~1\7-Zip\7z.exe e -y ' dirname filesep Fname ];
    [status,result] = system(uzcommand,'-echo'); 
    if (status~=0)
        error([ '*.gz file contents did not uncompress in NIMROD_wintest1' ...
            result ]);
    end
    
    % read the uncompressed NIMROD composite uk-1km rain radar data file
    [int_gen_hd, rl_gen_hd, rl_datsp_hd, char_hd, int_datsp_hd, ...
    rr_dat_mat] = rdnim1km( Zfname );
    
    % delete uncompressed file to save space
    delete( Zfname );
    % Check orientation is correct by making yellow box in top right hand 
    % corner
    xmax = int_gen_hd(19);
    rr_dat_mat(1:20,xmax-19:xmax) = 600;
    
    % Plot rain radar image on National Grid (NG) scale
    % -------------------------------------------------
    % Set colour limits for image and colourbar (1000=31.25mm/hr - you may
    % need to increase this for extreme weather events).
    clims = [ 0 1000 ];
    % create x and y vectors one value for each pixel along edge of plot
    %NB: yg is deliberately in descending order to allow correct display of
    % both array and NG co-ordinates on plot y-axis
    xg = linspace(rl_datsp_hd(2),rl_datsp_hd(4),int_gen_hd(19));
    yg = linspace(rl_datsp_hd(1),rl_datsp_hd(5),int_gen_hd(18));
    % display the rain radar data image
    imagesc(xg,yg,rr_dat_mat,clims);
    % use rainbow colour scheme
    colormap(jet);
    % display axis scales corresponding to...
    % NG co-ordinates of area covered by image: left right bottom top
    axis([rl_datsp_hd(2) rl_datsp_hd(4) rl_datsp_hd(5) rl_datsp_hd(1) ]);
    axis xy
    % ensure square grid i.e 1km easting is same as 1km northing
    axis equal
    % display uncompressed filename (which includes timestamp) as title
    title(Zfname);
    % display key to colours on rhs of plot
    colorbar;
    % allow successive plots to use same figure and plot area
    hold on
    
    % mark location of interest with white x
    plot(loc_easting,loc_northing,'xw');
    % force display to update
    drawnow
    % allow brief time for interrupts to be serviced
    pause(0.1);
            
end

close(fh1);
fprintf( 'Directory NIMROD file read complete: %d files processed\n',...
    No_files);
close all
clear all


end
