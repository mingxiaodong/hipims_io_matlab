function [int_gen_hd, rl_gen_hd, rl_datsp_hd, char_hd, int_datsp_hd, ...
    rr_dat_mat] = rdnim1km(fl_name)

%AUTHOR: Chethana Nagaraja (adapted by Andrew Duncan)
%INSTITUTION: University of Glamorgan (University of Exeter)
%DATE: 14-12-2007 (2010-11-05)
%TITLE: Reading binary data file from Nimrod database on BADC (British
%Atmospheric Data Centre)
% Downloaded from:
% http://badc.nerc.ac.uk/browse/badc/ukmo-nimrod/software/Matlab/rdnim.m
% Modified by Andrew Duncan
% University of Exeter, Centre for Water Systems
% http://student.secam.ex.ac.uk/staff-directory.dhtml
% 2010-10-29 This version has been tested to work with 
% http://badc.nerc.ac.uk/browse/badc/ukmo-nimrod/data/composite/uk-1km
% data files:
% Matlab fread statement fills array columnwise, whereas NIMROD uk-1km
% file stores data row-wise, so [data_c, data_r] have been reversed, then
% matrix transposed afterward to correct this difference.
%
%function 'rdnim.m' reads rain rate from nimrod data files stored on
%BADC. input argument, 'fl_name' is name of data file which
%should be a '.dat' file. function returns following parameters:
%1. 'int_gen_hd': This is general header section of type integer
%2. 'rl_gen_hd': This is general header section of type real
%3. 'rl_datsp_hd': This is data specific header section of type real
%4. 'char_hd': This is header section of type character
%5. 'int_datsp_hd': This is data specific header section of type integer
%6. 'rr_dat_mat': This is rain rate stored in matrix
%All above parameters except 'int_datsp_hd' are in Big Endian format
%'int_datsp_hd' is in Little endian format
%
%more detailed description of above header entries are available
%from Nimrod Dataset webpage on BADC. 
%These header entries must be comprehended considering following:
%1.first 2 elements of 'int_gen_hd' are '0' & '512' respectively. 
%2. last 2 elements of 'int_datsp_hd' are also '512' & '0' respectively.
%
%ALWAYS CHECK FORMAT OF MACHINE, i.e., LITTLE / BIG ENDIAN FORMATS.
%This can be done by checking to see if first 2 elements of the
%'int_gen_hd' parameter are returned as '0' & '512' respectively. If there
%is a mismatch, then machine format can be appropriately changed as
%explained below:
%byte swapping between above 2 formats can be performed by specifying
%'ieee-le' for little endian format & 'ieee-be' for big endian
%format in 'fread' function

ele_int_gen = 33;%number of elements in 'int_gen_hd'

%byte number where 'rl_gen_hd' begins, from beginning of data file
st_rl_gen = 66;
%number of elements in 'rl_gen_hd'
ele_rl_gen = 28;

%byte number where 'rl_datsp_hd' begins, from beginning of data file
st_rl_datsp = 178;
%number of elements in 'rl_datsp_hd'
ele_rl_datsp = 45;

%byte number where 'char_hd' begins, from beginning of data file
st_char = 358;
%number of characters in 'char_hd'
ele_char = 56;

%byte number where 'int_datsp_hd' begins, from beginning of data file
st_int_datsp = 415;
%number of elements in 'int_dat_sp'
ele_int_datsp = 53;

%byte number where data matrix begins, from beginning of data file
st_rr_dat = 524;

fid = fopen(fl_name,'r');
if fid == -1
    error('File does not exist');
end

int_gen_hd = fread(fid,ele_int_gen,'int16','ieee-be');
data_r = int_gen_hd(18);%row dimension of data matrix
data_c = int_gen_hd(19);%column dimension of data matrix

fseek(fid,st_rl_gen,'bof');
rl_gen_hd = fread(fid,ele_rl_gen,'float32','ieee-be');

fseek(fid,st_rl_datsp,'bof');
rl_datsp_hd = fread(fid,ele_rl_datsp,'float32','ieee-be');

fseek(fid,st_char,'bof');
char_hd = fread(fid,ele_char,'*char');

fseek(fid,st_int_datsp,'bof');
int_datsp_hd = fread(fid,ele_int_datsp,'int16','ieee-le');

% read main body array of rainfall intensity values
fseek(fid,st_rr_dat,'bof');
rr_dat_mat = fread(fid,[data_c, data_r],'int16','ieee-be');
% data is stored in file row-wise, whereas fread fills array columnwise, so
% transpose array
rr_dat_mat = rr_dat_mat';

fclose(fid);