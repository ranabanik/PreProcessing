disp('Enter 1 for CT, 2 for MR and 3 for MRCT');
prompt = 'Insert number: '
x = input(prompt);
%
if x == 1
    mod = 'VCT';
    VCT = niftiread('CT-MRI-5634_CT_norm.nii');
elseif x ==2
    mod = 'VMR';
%     VMR = niftiread('Smith_137398.03.01.10-47-08.WIP_T2W_DRIVE_CLEAR.01.nii');
    VMR = mristack;
elseif x == 3
    mod = 'VMRCT';
    VMRCT = niftiread('sMT0160-0013-00001-000001-02.nii');
%     VMRCT = niftiread('CT-MRI-5634_MR_norm_lab_2_.nii');
else 
    mod = 'MRmask'
    MRmask = BW;

end

if x == 1
    [sag, cor, axial] = size(VCT);  %row, column and number of slices
elseif x ==2
    [sag, cor, axial] = size(VMR);
elseif x == 3
    [sag, cor, axial] = size(VMRCT);
else
    [sag, cor, axial] = size(MRmask);
end

t = 1;
while t<10
    t=t+1;
    for i = 1:axial
        if x == 1            
            imagesc(squeeze(VCT(:,:,i)));
        elseif x == 2
            imagesc(squeeze(VMR(:,:,i)));
        elseif x ==3
            imagesc(squeeze(VMRCT(:,:,i)));
        else
            imagesc(squeeze(MRmask(:,:,i)));
        end
    if (i == 1)
        cLim_v = get(gca, 'CLim');
    else
        set(gca, 'CLim', cLim_v)
    end
    axis image
    axis on
    colormap(gray)
if x == 1    
    title(['CT data: ', num2str(i)]);
elseif x ==2
    title(['MR data: ', num2str(i)]);
else
    title(['MR/CT data: ', num2str(i)]);
end
    drawnow
    mov(i) = getframe;
  
    end
end
%% CT image voxel
voxCT = V(117,144,100)
%% MR image voxel
voxMR = V(151,145,91)
%%
i = 100;
% for nSlice = 1:axial 
  for nCols = 1:cor
    for nRows = 1:sag
        plane(nRows,nCols) = V(nRows,nCols,i); 
    end
  end
%% Converting .gz to .nii
imds = imageDatastore(cd,'IncludeSubfolders',true,'LabelSource','FolderNames','FileExtensions','.gz');
nGZ = numel(imds.Files);
for i = 1:nGZ
    gunzip(imds.Files{i});
    clear imds.Files{i};
end
%% Create Mask(Main code)
% Create mask from 2D MR images and then implement that on CT 
% runs with the Tiff.m file in matlab. 
% The folders are arranged as below
% (base) D:\MRCTdata_2D_tif_256_256_single_int32>
% ????Test
% ?   ????CT
% ?   ????Masked_CT
% ?   ????Masked_MR
% ?   ????MR
% ????Train
% ?   ????CT
% ?   ????Masked_CT
% ?   ????Masked_MR
% ?   ????MR
% ????Valid
%     ????CT
%     ????Masked_CT
%     ????Masked_MR
%     ????MR
clear all
close all
Dir = uigetdir; %select the folder where MR and CT branch off
mod1 = '\MR';
mod2 = '\CT';
masked_MR_dir = '\Masked_MR';
masked_CT_dir = '\Masked_CT';
mrDir = strcat(Dir,mod1); %'D:\MRCTdata_2D_tif_256_256_single_int32\trial\MR'
ctDir = strcat(Dir,mod2);
mrFiles = dir(fullfile(mrDir,'*.tif'));
ctFiles = dir(fullfile(ctDir,'*.tif'));

mask = zeros(256,256);

for i = 1:length(mrFiles)
    baseFileName = mrFiles(i).name; %'2809.tif'
    MRfullFileName = fullfile(mrDir,baseFileName); %'D:\MRCTdata_2D_tif_256_256_single_int32\trial\MR\2809.tif'
    CTfullFileName = fullfile(ctDir,baseFileName);
    [folder, name, extension] = fileparts(MRfullFileName);
    MR = imread(MRfullFileName);
    CT = imread(CTfullFileName);
    
    for c = 1:size(MR,2)
        for r = 1: size(MR,1)
            if MR(r,c) > 4.9
                mask(r,c) = 1;
            else
            mask(r,c) = 0;
            end
        end
    end
    
    filled = imfill(mask,'holes');
    w = 3;
    SE = strel('square',w);
    fill_er = imerode(filled,SE);
    masked_MR = MR.*fill_er;
    masked_CT = int32(zeros(size(CT)));
    for c = 1:256
        for r = 1:256
            if fill_er(r,c) == 1
                masked_CT(r,c) = CT(r,c);
            else
                masked_CT(r,c) = -1010;
            end
        end
    end
    
    string = strcat('Masked_CT/',name,'.tif');
    finalCT = Tiff(string,'w');
    
    tagstruct.ImageLength = size(masked_CT,1);
    tagstruct.ImageWidth = size(masked_MR,2);
    tagstruct.SampleFormat = 2; % uint
    tagstruct.Photometric = Tiff.Photometric.MinIsBlack;
    tagstruct.BitsPerSample = 32 ;
    tagstruct.SamplesPerPixel = 1;
    tagstruct.YCbCrSubSampling = [1,1];
    tagstruct.Compression = Tiff.Compression.None;  
    tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky; 
    tagstruct.Software = 'MATLAB'; 
    
    setTag(finalCT,tagstruct);
    write(finalCT,masked_CT);
    close(finalCT);
    
    string = strcat('Masked_MR/',name,'.tif');
    finalMR = Tiff(string,'w');
    
    tagstruct.ImageLength = size(masked_CT,1);
    tagstruct.ImageWidth = size(masked_CT,2);
    tagstruct.SampleFormat = 3; % saves image in single datatype
    tagstruct.Photometric = Tiff.Photometric.MinIsBlack;
    tagstruct.BitsPerSample = 32;
    tagstruct.SamplesPerPixel = 1;
    tagstruct.YCbCrSubSampling = [1,1];
    tagstruct.Compression = Tiff.Compression.None;  
    tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky; 
    tagstruct.Software = 'MATLAB'; 
    
    setTag(finalMR,tagstruct);
    write(finalMR,masked_MR);
    close(finalMR);
end

