function [ svStruct ] = readSupervoxelFromBinary (fid )
%Reads supervoxel from binary file with identifier fid written by C++ file supervoxel:: writeToBinary

dimsImage = 3;

svStruct.TM = fread(fid, 1, 'int32');
svStruct.dataSizeInBytes = fread(fid, 1, 'uint64');
svStruct.dataDims = fread(fid, dimsImage, 'uint64');

ll = fread(fid, 1, 'uint32');
if ( ll > 0)
    svStruct.PixelIdxList = fread(fid, ll, 'uint64');
else
    svStruct.PixelIdxList = [];

end

