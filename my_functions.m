
function [img, info] = load_patient(PatientName)
    
    img_fullpath = append(dir_fullpath, PatientName, '-integra.dcm');
    img = dicomread(dir(img_fullpath).name);
    img = double(img);
    info = dicominfo(dir(img_fullpath).name);
    
end
