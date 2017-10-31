% Date:29/10/2017
% This code refers to the code written by Changkai Zhao
% Email:changkaizhao1006@gmail.com
% completed and corrected in 20/06/2013 and in 23/06/2013
%
% This algorithm is described in details in
%
% "Single Image Haze Removal Using Dark Channel Prior",
% by Kaiming He Jian Sun Xiaoou Tang,
% In: CVPR 2009

% details about guilded image filter in
% "Guilded image filtering"
% by Kaiming He Jian Sun Xiaoou Tang

% OUTPUT:
% J is the obtained clear image after visibility restoration.
% tmap is the raw transmission map.
% tmap_ref is the refined transmission map.

function [J, tmap, tmap_ref] = darkChannel(I, px, w)

    % calculate running time
    tic;

    % by default, constant parameter w is set to 0.95.
    if (nargin < 3)
       w = 0.95;
    end

    % by default, the patch size is set to 15.
    if (nargin < 2)
       px = 15;
    end

    % error, no argument.
    if (nargin < 1)
    	msg1 = sprintf('%s: Not input.', upper(mfilename));
        eid = sprintf('%s:NoInputArgument',mfilename);
        error(eid,'%s %s',msg1);
    end

    % error, parameter out of bound.
    if ((w >= 1.0) || (w<=0.0))
    	msg1 = sprintf('%s: w is out of bound.', upper(mfilename));
        msg2 = 'It must be an float between 0.0 and 1.0';
        eid = sprintf('%s:outOfRangeP',mfilename);
        error(eid,'%s %s',msg1,msg2);
    end

    % error, px out of bound.
    if (px < 1)
    	msg1 = sprintf('%s: px is out of bound.', upper(mfilename));
        msg2 = 'It must be an integer higher or equal to 1.';
        eid = sprintf('%s:outOfRangeSV',mfilename);
        error(eid,'%s %s',msg1,msg2);
    end

    % Pick the top 0.1% brightest pixels in the dark channel.
    Image = im2double(I);
    [dimr, dimc, col] = size(I);
    dx = floor(px / 2);

    % L0Smoothing
    % Image = L0Smoothing(Image, 0.003);

    % Initial four matrices
    J          = zeros(dimr, dimc, col);
    t_map      = zeros(dimr, dimc);
    J_darktemp = zeros(dimr, dimc);
    tmap_ref   = zeros(dimr, dimc);

    if(col == 3)
        A_r = 0; A_g = 0; A_b = 0;

        % Estimate the atmospheric light
        J_darkchannel = min(Image, [], 3);  % find the minimum value in RGB channel
        for i = (1 : dimr)
            for j = (1 : dimc)
                winLeft = i - dx; winRight = i + dx;
                winUp   = j - dx; winDown  = j + dx;

                % check the windows range
                if(i - dx < 1)
                    winLeft = 1;
                end
                if(i + dx > dimr)
                    winRight = dimr;
                end
                if(j - dx < 1)
                    winUp = 1;
                end
                if(j+dx>dimc)
                      winDown=dimc;
                end

                % find the minimum value in the patch windows
                J_darktemp(i,j) = min(min(J_darkchannel(winLeft : winRight, winUp : winDown)));
            end
        end
        J_darkchannel = J_darktemp;

        % get the 0.1% most brightest pixels in the dark channel
        lightLimit = quantile(J_darkchannel(:), [.999]);
        [lightRow, lightCol] = find(J_darkchannel >= lightLimit);
        [enum, ~] = size(lightRow);
        for i = (1 : enum)
            A_r = Image(lightRow(i), lightCol(i), 1) + A_r;
            A_g = Image(lightRow(i), lightCol(i), 2) + A_g;
            A_b = Image(lightRow(i), lightCol(i), 3) + A_b;
        end

        % get the Airlight
        Airlight = [A_r / enum, A_g / enum, A_b / enum];

        % Estimating the raw transmission map(color)
        Im_n(:, :, 1) = Image(:, :, 1) ./ Airlight(1);
        Im_n(:, :, 2) = Image(:, :, 2) ./ Airlight(2);
        Im_n(:, :, 3) = Image(:, :, 3) ./ Airlight(3);
        tmap = min(Im_n, [], 3);

        for i = (1 : dimr)
            for j = (1 : dimc)
                winLeft = i - dx; winRight = i + dx;
                winUp   = j - dx; winDown  = j + dx;

                 % check the windows range
                if(i - dx < 1)
                    winLeft = 1;
                end
                if(i + dx > dimr)
                    winRight = dimr;
                end
                if(j - dx < 1)
                    winUp = 1;
                end
                if(j + dx > dimc)
                    winDown = dimc;
                end

                % find the minimum value in the patch windows
                % get the image transmittance
                t_map(i, j) = 1 - w * min(min(tmap(winLeft : winRight, winUp : winDown)));
            end
        end
        
        % Refine the raw transmission map(color)
        % using softmatting
        tmap_ref = softmatting(Image, t_map);

        % using guidedfilter_color
        % tmap_ref = guidedfilter_color(Image, t_map, 40, 0.001);

        % Getting the clear image(color)
        % set the lowest t0
        [lightRow, lightCol] = find(tmap_ref < 0.1);
        [enum, ~] = size(lightRow);
        for i = (1 : enum)
            tmap_ref(lightRow(i), lightCol(i)) = 0.1;
        end

        J(:, :, 1) = (Image(:, :, 1) - Airlight(1)) ./ tmap_ref + Airlight(1);
        J(:, :, 2) = (Image(:, :, 2) - Airlight(2)) ./ tmap_ref + Airlight(2);
        J(:, :, 3) = (Image(:, :, 3) - Airlight(3)) ./ tmap_ref + Airlight(3);

        % show the result
        figure,imshow(Image),title('Input Image');
        figure,imshow(t_map),title('Raw t map');
        figure,imshow(tmap_ref),title('Refined t map');
        figure,imshow(J),title('Output Image');
    end

    % stop calculate running time
    toc;
end