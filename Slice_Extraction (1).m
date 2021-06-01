clear;

%% prof functions

addpath(genpath('/Users/alessandro/Desktop/Medical Imaging/thirdparty-libraries'));

%% Load all images

files = dir('*.dcm');
size(dir('/Users/alessandro/Desktop/untitled folder/*.dcm'))
info = dicominfo(files(1).name);
info.PixelSpacing
info.SliceThickness

for n = 1:size(files,1)
    pathtemp = files(n).name;
    info__temp = dicominfo(pathtemp);
    dcm__temp = dicomread(pathtemp);
    dcm__temp = double(dcm__temp.*info__temp.RescaleSlope);
    stackimg(n,:,:) = dcm__temp;
    disp(['The dimension of stackimg after ' num2str(n)...
        ' iterations is ' num2str(size(stackimg))]);
end
nvoxel = size(stackimg,2) * size(stackimg,3);
disp(['Each one of the ' num2str(size(stackimg,1))...
    ' loaded images is made of ' num2str(nvoxel) ' voxels']);



% Normalize the image intensity between 0 and 1 (just for visualizing
% purposes!)
minimum = min(min(min(stackimg)));
normalized_stackimg = stackimg - minimum;
maximum = max(max(max(stackimg)));
normalized_stackimg = normalized_stackimg/maximum;
normalized_stackimg = normalized_stackimg*1;

% The new calculated min/max values are:
minimum = min(min(min(stackimg))); disp(num2str(minimum));
maximum = max(max(max(stackimg))); disp(num2str(maximum));

%% Tumor search: visualize all slices of the entire volume

% From the top (with mirror effect)

figure();
for n = 150:size(stackimg,1)
    imshow(squeeze(normalized_stackimg(n,:,:)), 'InitialMagnification', 600);
    title(n,'FontSize', 16);
    pause(1);
end

% Full frontal

figure();
for n = 320:size(stackimg,2)
    imshow(squeeze(normalized_stackimg(:,n,:)), 'InitialMagnification',600); 
    title(n,'FontSize', 18);
    pause(0.5);
end

% Lato
figure();
for n = 1:size(stackimg,3)
    imshow(squeeze(normalized_stackimg(:,:,n)), 'InitialMagnification',800);
    title(n,'FontSize', 18);
    pause(0.2);
end



%% Slice of Interest

img = squeeze(stackimg(147,:,:));

figure();
imshow(img, 'InitialMagnification',600);
caxis('auto');

%% Write Slice of Interest

fullpath = '/Users/alessandro/Desktop/Medical Imaging/progetto_new/Lesioni Segmentate/';

class(dicomread(files(1).name))

dicomwrite(int16(img), append(fullpath, info.PatientName.FamilyName, '-integra', '.dcm'), info);
