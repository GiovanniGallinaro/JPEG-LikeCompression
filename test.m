close all
clearvars

prompt = "Select an image (type the index):\n";
dir dataset/*.png
%image_index = input(prompt, 's');

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

    disp(' ')

    % get image size
    [rows, cols, ~] = size(input);

    % get image dimension
    bits = numel(dec2bin(input));
    dim = bits/8000; % [kB]
    n_pixels = rows*cols;
    bpp = bits/n_pixels;

    disp(['Original image dimension: ' num2str(dim) ' kB'])
    disp(['bpp: ' num2str(bpp)])
    fprintf(fileID, 'IMAGE %2.0f\n', image_index(i));
    fprintf(fileID, 'Original image dimension: %4.2f kB\n\n', dim);

    padding_flag = false;
    % pad image if dimensions are not multiple of 8
    if(mod(rows, 8)~=0)
        input = [input; zeros(8 - mod(size(input, 1), 8), size(input, 2), 3)];
        padding_flag = true;
    end
    if(mod(cols, 8)~=0)
        input = [input; zeros(size(input, 1), 8 - mod(size(input, 2), 8), 3)];
        padding_flag = true;
    end

    tic

    % convert to YCbCr
    ycbcr = rgb2ycbcr(input);

    % extract components and center them to zero
    Y = double(ycbcr(:, :, 1)) - 128;
    Cb = double(ycbcr(:, :, 2)) - 128;
    Cr = double(ycbcr(:, :, 3)) - 128;

    % chroma subsampling (4:2:0)
    Cb = Cb(1:2:end, 1:2:end);
    Cr = Cr(1:2:end, 1:2:end);

    % compute 8x8 2D dct
    dct_2D = @(block_struct) dct2(block_struct.data);
    % compute 8x8 2D dct
    Y_dct = blockproc(Y, [8, 8], dct_2D);
    Cb_dct = blockproc(Cb, [8, 8], dct_2D);
    Cr_dct = blockproc(Cr, [8, 8], dct_2D);

    quality_factor = [15 30 50 65 80 100];
    for j=1:length(quality_factor)
        [Q_Y, Q_C] = compute_tables(quality_factor(j));

        quant_Y = @(block_struct) floor((block_struct.data ./ Q_Y) + 0.5);
        quant_C = @(block_struct) floor((block_struct.data ./ Q_C) + 0.5);

        Y_dct_quant = blockproc(Y_dct, [8, 8], quant_Y);
        Cb_dct_quant = blockproc(Cb_dct, [8, 8], quant_C);
        Cr_dct_quant = blockproc(Cr_dct, [8, 8], quant_C);

        % encoding
        [Y_encoded, dict_Y] = encode(Y_dct_quant);
        [Cb_encoded, dict_Cb] = encode(Cb_dct_quant);
        [Cr_encoded, dict_Cr] = encode(Cr_dct_quant);

        % compute the dimension of the encoded representation
        enc_dim = length(Y_encoded) + length(Cb_encoded) + length(Cr_encoded);
        enc_dim = enc_dim/8000;  % [kB]
        disp(['Encoded image dimension: ' num2str(enc_dim) ' kB'])
        
        fprintf(fileID, 'QUALITY %3.0f -------------------------------------\n', quality_factor(j));
        fprintf(fileID, 'Encoded image dimension: %4.2f kB\n', enc_dim);

        % computation of the bpp for the compressed image
        bpp_rec = (enc_dim*8000)/n_pixels;
        disp(['bpp: ' num2str(bpp_rec)])
        
        fprintf(fileID, 'bpp: %4.2f\n', bpp_rec);

        % decoding
        Y_decoded = decode(Y_encoded, dict_Y, rows, cols);
        Cb_decoded = decode(Cb_encoded, dict_Cb, rows/2, cols/2);
        Cr_decoded = decode(Cr_encoded, dict_Cr, rows/2, cols/2);

        % reconstruction step
        recon_Y = @(block_struct) block_struct.data .* Q_Y;
        recon_C = @(block_struct) block_struct.data .* Q_C;

        Y_dct_rec = blockproc(Y_decoded, [8, 8], recon_Y);
        Cb_dct_rec = blockproc(Cb_decoded, [8, 8], recon_C);
        Cr_dct_rec = blockproc(Cr_decoded, [8, 8], recon_C);

        % compute 8x8 2D idct
        idct_2D = @(block_struct) idct2(block_struct.data);

        Y_rec = blockproc(Y_dct_rec, [8,8], idct_2D);
        Cb_rec = blockproc(Cb_dct_rec, [8,8], idct_2D);
        Cr_rec = blockproc(Cr_dct_rec, [8,8], idct_2D);

        % chroma luminance interpolation (original size)
        Cb_rec = imresize(Cb_rec, 2, 'nearest');
        Cr_rec = imresize(Cr_rec, 2, 'nearest');

        % create YCbCr image and convert it to RGB
        image_ycbcr_rec = zeros(rows, cols, 3);
        image_ycbcr_rec(:, :, 1) = Y_rec;
        image_ycbcr_rec(:, :, 2) = Cb_rec;
        image_ycbcr_rec(:, :, 3) = Cr_rec;
        image_ycbcr_rec = uint8(image_ycbcr_rec + 128); % translate values by 128
        image_rgb_rec = ycbcr2rgb(image_ycbcr_rec);

        toc
       
        time = toc;
        fprintf(fileID, 'Elapsed time: %4.2f seconds\n', time);

        % if there was padding, cut the image
        if (padding_flag)
            image_rgb_rec = image_rgb_rec(1:rows_original, 1:cols_original, :);
        end

        compression_ratio = dim/enc_dim;
        disp(['Compression ratio: ' num2str(compression_ratio)])
        
        fprintf(fileID, 'Compression ratio: %4.2f\n', compression_ratio);

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

        % imshowpair(input, image_rgb_rec, 'montage');

        %imshow(image_rgb_rec);
        filename = ['out\' num2str(image_index(i)) 'q' num2str(quality_factor(j)) '.jpeg'];
        imwrite(image_rgb_rec, filename);
    end
end

fclose(fileID);

function [Q_Y, Q_C] = compute_tables(quality_factor)
    Q_Y = [...
        16 11 10 16 24  40  51  61;
        12 12 14 19 26  58  60  55;
        14 13 16 24 40  57  69  56;
        14 17 22 29 51  87  80  62;
        18 22 37 56 68  109 103 77;
        24 36 55 64 81  104 113 92;
        49 64 78 87 103 121 120 101;
        72 92 95 98 112 100 103 99];
    Q_C = [...
        17 18 24 47 99 99 99 99;
        18 21 26 66 99 99 99 99;
        24 26 56 99 99 99 99 99;
        47 66 99 99 99 99 99 99;
        99 99 99 99 99 99 99 99;
        99 99 99 99 99 99 99 99;
        99 99 99 99 99 99 99 99;
        99 99 99 99 99 99 99 99];
    % if not, make it integer
    quality_factor = round(quality_factor);
    % check range
    if quality_factor < 1
        quality_factor = 1;
    else
        if quality_factor > 100
            quality_factor = 100;
        end
    end
    Q_Y = round(Q_Y * 1/quality_factor * 50);
    Q_C = round(Q_C * 1/quality_factor * 50);
    % turn zero-valued entries into 1s
    Q_Y(Q_Y == 0) = 1;
    Q_C(Q_C == 0) = 1;
end

function [symbols, prob] = generate_dictionary(quantized_coefficients)
    image_flattened = single(quantized_coefficients(:));
    symbols = unique(image_flattened);
    counts = hist(image_flattened, symbols);
    prob = counts / sum(counts);
    if (length(symbols) == 1) % to cope with all equals elements
        prob = 1;
    end
end

function [encoded, dict] = encode(image_grayscale)
    [~, ~, number_of_channels] = size(image_grayscale);
    if (number_of_channels > 1)
        error('The input image must have only 1 channel!')
    end
    [symbols, prob] = generate_dictionary(image_grayscale);
    if (length(symbols) > 1)
        dict = huffmandict(symbols, prob); % created manually
    else
        dict = {symbols, 0};
    end
    image_flattened = image_grayscale(:);
    encoded = huffmanenco(image_flattened, dict);
end

function decoded = decode(encoded_vector, dict, rows, cols)
    if (~isvector(encoded_vector))
        error('The input must be a vector!')        
    end
    decoded = huffmandeco(encoded_vector, dict);
    decoded = reshape(decoded, [rows, cols]);
end






