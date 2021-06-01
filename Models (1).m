%% Haralick feature extraction

clear;

cd('/Users/alessandro/Desktop/Medical Imaging/progetto_new/Lesioni Segmentate/');
addpath(genpath('/Users/alessandro/Desktop/Medical Imaging/thirdparty-libraries'));
addpath(genpath('/Users/alessandro/Downloads/ADASYN_upd2'));

files = dir('/Users/alessandro/Desktop/Medical Imaging/progetto_new/Lesioni Segmentate/*1.dcm');

vars = [];

for n = 1:size(files,1)
    filtered_img = dicomread(files(n).name);
    filtered_img = double(filtered_img);
    PatientID = dicominfo(files(n).name).PatientName.FamilyName;
    m = graycomatrix(filtered_img);
    single_vars = getGLCMtextures(m);
    PatientID = str2num(extractAfter(PatientID,"-"));
    %single_vars.PatientID = PatientID;
    vars = [vars single_vars];
end

labels = [0 0 1 0 1 1 0 0 0 1 0 0]; % solido 0 --- sub-solido 1

table_vars = struct2table(vars);
table_vars.labels = labels';


%% Export for knime

writetable(table_vars, '/Users/alessandro/Desktop/original_data.xlsx');


%% Models

vars = rmfield(vars, 'PatientID');

vars_matrix = squeeze(cell2mat(struct2cell(vars)))';



[oversampled_vars, oversampled_labels] = ADASYN(vars_matrix, labels, 1, 2, 3);

oversampled_data = vars_matrix;
oversampled_data(13:17,:) = oversampled_vars;

labels = [labels oversampled_labels'];

% Export

export_data = oversampled_data;
export_data(:,10) = labels;
writematrix(export_data, '/Users/alessandro/Desktop/data.csv'); 

%% PCA

[COEFF, SCORE, LATENT, TSQUARED, EXPLAINED] = pca(oversampled_data);

EXPLAINED

new_data = SCORE(:,1);

% SVM, RF, NN, DT

% Estrarre il modello 

SVMModel = fitcsvm(new_data, labels, 'Standardize', false);

CVSVMModel = crossval(SVMModel, 'Leaveout', 'on');

supervised__accuracy = 1 - kfoldLoss(CVSVMModel, 'Mode', 'average');




vars







