clear;

%% prof functions

addpath(genpath('/Users/alessandro/Desktop/Medical Imaging/thirdparty-libraries'));

%% Load all images

files = dir('*.dcm');

info = dicominfo(files(1).name)

for n = 80:90%size(files,1)
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
normalized_stackimg = stackimg/maximum;
normalized_stackimg = normalized_stackimg*1;

% The new calculated min/max values are:
minimum = min(min(min(stackimg))); disp(num2str(minimum));
maximum = max(max(max(stackimg))); disp(num2str(maximum));

%% Tumor search: visualize all slices of the entire volume

% From the top (with mirror effect)

figure();
for n = 80:90%size(stackimg,1)
    imshow(squeeze(normalized_stackimg(n,:,:)), 'InitialMagnification', 600);
    title(n,'FontSize', 16);
    pause(0.5);
end

% Full frontal

figure();
for n = 190:size(stackimg,2)
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

img = squeeze(stackimg(85,:,:));

figure();
imshow(img, 'InitialMagnification',600);
caxis('auto');

%% Write Slice of Interest

fullpath = '/Users/alessandro/Desktop/Medical Imaging/progetto_new/Lesioni Segmentate/';

class(dicomread(files(1).name))

dicomwrite(int16(img), append(fullpath, info.PatientName.FamilyName, '---', sprintf('%03d', 1), 'integra', '.dcm'), info);

%% ROI extraction

clearvars -except fullpath;



h = imfreehand;

pos = getPosition(h); 

mask = createMask(h);
cutted_img = img.*mask;

figure();
imshow(cutted_img,[]);
caxis('auto');

%%%%% Filter

% Calculate the threshold
max_cutted = max(max(cutted_img));
threshold = 0.45* max_cutted;

% Apply the threshold to the ROI

cutted_img_roi = cutted_img;
cutted_img_roi(cutted_img_roi < threshold) = 0;

figure();
imshow(cutted_img_roi);
caxis('auto');

%% Export lesions

fullpath = '/Users/alessandro/Desktop/Medical Imaging/progetto_new/Lesioni Segmentate/';

class(dicomread(files(1).name))

dicomwrite(int16(cutted_img_roi), append(fullpath, info.PatientName.FamilyName, '---', sprintf('%03d', 1), '.dcm'), info);
dicomwrite(int16(img), append(fullpath, info.PatientName.FamilyName, '---', sprintf('%03d', 1), 'integra', '.dcm'), info);

%% Haralick feature extraction

clear;

fullpath = '/Users/alessandro/Desktop/Medical Imaging/progetto_new/Lesioni Segmentate';

addpath(genpath('/Users/alessandro/Desktop/Medical Imaging/thirdparty-libraries'));

cd(fullpath);

files = dir('*1.dcm');


vars = []

for n = 1:size(files,1)
    img = dicomread(files(n).name);
    m = graycomatrix(img);
    single_vars = getGLCMtextures(m);
    vars = [vars single_vars];
end

vars = reshape(struct2array(vars), [numel(fieldnames(vars)), size(files,1)])'

labels = [1 0 1 1 0 1 1 0 0 0 1 1]

%vars(:,(size(vars,2)+1)) = labels;

idx = kmeans(vars, 2,'Replicates',1000); % che ci facciamo?
% vars(:,(size(vars,2)+1)) = idx

SVMModel = fitcsvm(vars, labels, 'Standardize', false);
CVSVMModel = crossval(SVMModel,'Leaveout','on');
supervised__accuracy = 1 - kfoldLoss(CVSVMModel);













%% K-means

% Use k-means clustering to perform automatic segmentation and clean images

% Non viene fatto slice per slice ma tutto insieme, non dovrebbe cambiare
% molto, solo che facendo di slice in slice ogni volta cambiava le
% etichette e quindi in alcune slice la lesione era bianca e in altre era
% nera... un casino! Così assegna le etichette in modo coerente.

reshaped_cutted_img = reshape(cutted_img,[size(cutted_img,1)*size(cutted_img,2),1,1]);

idx = [];
idx_reshaped = [];
idx(:) = kmeans(reshaped_cutted_img(:,:), 2,'replicate', 10);
idx_reshaped(:,:) = reshape(idx(:),[size(cutted_img,1),size(cutted_img,2)]);

minimum = min(min(idx_reshaped)); 
idx_reshaped = idx_reshaped - minimum;

maximum = max(max(idx_reshaped)); 
idx_reshaped = idx_reshaped/maximum;
idx_reshaped = idx_reshaped*1;

figure();
for n = 1:size(idx_reshaped,3)
    imshow(squeeze(idx_reshaped(:,:,n)), 'InitialMagnification',600);
    pause(0.5);
end



%% morphological features
lesion__volume = sum(sum(sum(stackimg))) * dim__voxel
lesion__radius = (lesion__volume *3/(4*pi))^(1/3);
lesion__diameter = lesion__radius * 2 % il json dice che il tumor_size_in_cm è 4.2 !!!!!!! :) :)
[lesion__area, surf_mat] = compute__surface(stackimg, double([10 * dimx, 10 * dimy, 10 * dimz]));
lesion__area = double(lesion__area)% cm^2
spherical__disproportion = lesion__area/(4*pi*(lesion__radius)^2)
sphericity = ((pi^(1/3))*((6*lesion__volume)^(2/3)))/lesion__area
surfacevolume__ratio = lesion__area/lesion__volume
morphological__features = [lesion__volume lesion__area spherical__disproportion sphericity surfacevolume__ratio];


