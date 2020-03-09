clearvars
close all

image_index = ["05", "08", "12", "15", "19"];

fileID = fopen('out\res.txt', 'w');

for i=1:length(image_index)
    % load input image
    filename = strcat("dataset/kodim", image_index(i), ".png");

%     check if image exists
    if(isfile(filename))
        input = imread(filename);
    else
        error("Wrong file specified.");
    end
    
    [rows, cols, ~] = size(input);

    % get image dimension
    bits = numel(dec2bin(input));
    dim = bits/8000; % [kB]
    n_pixels = rows*cols;
    bpp = bits/n_pixels;
    
    disp(['bpp: ' num2str(bpp)])

    fprintf(fileID, 'IMAGE %2.0f\n', image_index(i));
    
    tic
    
    quality_factor = [5 15 30 50 65 80 100];
    for j=1:length(quality_factor)
        
        fprintf(fileID, 'QUALITY %3.0f -------------------------------------\n', quality_factor(j));
        
        filename = ['jpeg\' num2str(image_index(i)) 'q' num2str(quality_factor(j)) '.jpeg'];
        imwrite(input, filename, 'jpeg', 'Quality', quality_factor(j));
        
        toc
        time = toc;
        fprintf(fileID, 'Elapsed time: %4.2f seconds\n', time);
        
        image_rgb_rec = imread(['jpeg\' num2str(image_index(i)) 'q' num2str(quality_factor(j)) '.jpeg']);
        
        enc_bits = numel(dec2bin(image_rgb_rec));
        enc_dim = bits/8000; % [kB]
        
        disp(['Encoded image dimension: ' num2str(enc_dim) ' kB'])
        fprintf(fileID, 'Encoded image dimension: %4.2f kB\n\n', enc_dim);
        
        % compute the PSNR between the original and compressed image
        MSE_R = (double(input(:,:,1)) - double(image_rgb_rec(:,:,1))).^2;
        MSE_G = (double(input(:,:,2)) - double(image_rgb_rec(:,:,2))).^2;
        MSE_B = (double(input(:,:,3)) - double(image_rgb_rec(:,:,3))).^2;

        MSE_R = sum(MSE_R(:)) / n_pixels;
        MSE_G = sum(MSE_G(:)) / n_pixels;
        MSE_B = sum(MSE_B(:)) / n_pixels;

        MSE = (MSE_R+MSE_G+MSE_B)/3;

        PSNR = 10*log10( 255^2 / MSE);

        disp(['PSNR: ' num2str(PSNR) ' dB'])
        
        fprintf(fileID, 'PSNR: %4.2f dB\n\n', PSNR);
        
    end
end