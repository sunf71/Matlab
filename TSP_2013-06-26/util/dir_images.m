function [files] = dir_images(input_folder)

files = dir([input_folder '*.png']);
if (numel(files)==0)
    files = dir([input_folder '*.bmp']);
    if (numel(files)==0)
        files = dir([input_folder '*.jpg']);
        if (numel(files)==0)
            files = dir([input_folder '*.JPG']);
            if (numel(files)==0)
                files = dir([input_folder '*.ppm']);
            end
        end
    end
end
