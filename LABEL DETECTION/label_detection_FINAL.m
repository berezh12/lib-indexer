clear all;
img = imread('input2.jpg');
img_double = im2double(img);
img_grey = rgb2gray(img_double);

%test_img_threshold = (img_grey > 0.85) & (img_grey < 0.95);

%{
imshow(img_threshold);
figure(2);
se = strel('rectangle', [5,5]);
imshow(imdilate(img_threshold, se));
%}

%% CORNER DETECTION
img_grey_gaus = imgaussfilt(img_grey, 4);
corners = harris_corners(img_grey);
disp(size(corners(1)));
disp(size(corners(2)));
img_corner_highlight = img_grey;


%% INTEGRAL IMAGING
bright_th = img_grey;
bright_th = bright_th .^ 2.5 .* 3;
bright_th(bright_th > 0.7) = 1;
bright_th(bright_th < 0.3) = 0;
imshow(bright_th);
img_integral_brights = generate_integral_image(bright_th);

dark_th = 1-img_grey;
dark_th = dark_th .^ 2.75 .* 2.95;
dark_th(dark_th > 0.3) = 1;
dark_th(dark_th < 0.3) = 0;
imshow(dark_th);
img_integral_darks = generate_integral_image(dark_th);

%% FIND LABELS
brights_darks_ratio_range = [7, 60];   % ratio between II brightness sum/II darkness sum [min acceptable, max acceptable]
x_diff_range = [20, 125];               % absolute distance that two corners need to be apart in x direction [min, max]
y_diff_range = [100, 200];              % absolute distance that two corners need to be apart in y direction [min, max]
labels = [];
corners = sortrows(corners);

for c = 1:size(corners, 1)          % for every corner, find neighbors within x_diff_range and y_diff_range and check bright/dark ratio
    % quadtree to find neaarest neighbor with minimum distance of
    % [x_diff_range(1), y_diff_range(1)] and maximum distance of [x_diff_range(2), y_diff_range(2)]
    % if no other corner is in range, ignore the current corner and move on
    % if another corner is found in range, do integral imaging and
    % calculate the brightness/darkness ratio

    neighbors = find_neighbors(corners, c, x_diff_range, y_diff_range);
    for n = 1:size(neighbors, 1)
        if check_if_label(img_integral_brights, img_integral_darks, corners(c, :), neighbors(n, :), brights_darks_ratio_range)
            labels = [labels; corners(c, :), neighbors(n, :)];
        end
    end
    
end

% labels contains final result

%% DEBUG
imshow(img_grey);
hold on;
for t = 1:size(labels, 1)
    plot([labels(t, 2), labels(t, 4)], [labels(t, 1), labels(t, 3)], 'LineWidth', 2);
end
cornerst =[corners(:,2),corners(:,1)];
for c = 1:size(cornerst, 1)
    th = 0:pi/2:2*pi;
    r = 5;
    x = r * cos(th) + cornerst(c, 1);
    y = r * sin(th) + cornerst(c, 2);
    plot(x, y, 'r');
end
hold off;


%% FUNCTIONS
function result = find_neighbors(coords, coords_target_index, x_diff_range, y_diff_range)
    % find nearest neighbors in range
    
    result = [];
    for j = coords_target_index + 1:size(coords, 1)
        neighbor = coords(j, :);
        if check_for_range(coords(coords_target_index, :), neighbor, x_diff_range, y_diff_range)
            %disp(corners_x(i,:));
            %disp(next);
            %plot([corners_x(i,1), next(1,1)], [corners_x(i,2), next(1,2)]); 
            result = [result; neighbor];
        end
    end
    
end

%check if the two points are not too far from each other
function result = check_for_range(coordA, coordB, x_diff_range, y_diff_range)

    diff_x = abs(coordA(2) - coordB(2));
    diff_y = abs(coordA(1) - coordB(1));
    
    result = diff_x >= x_diff_range(1) && ...
             diff_x <= x_diff_range(2) && ...
             diff_y >= y_diff_range(1) && ...
             diff_y <= y_diff_range(2);
   
end

%{
    Returns true if the two target vectors could describe a label (based on
    a heuristic), or false if it is no label that can be detected.
    Sources:
        http://delivery.acm.org/10.1145/810000/808600/p207-crow.pdf
        accessed on 2019/11/12
        https://en.wikipedia.org/wiki/Summed-area_table#/media/File:Summed_area_table.png
        used only for variable naming
        accessed on 2019/11/12
    Author:
        Laurenz Edmund Fiala (11807869)
%}
function result = check_if_label(integral_image_brights, integral_image_darks, corner, neighbor, brights_darks_ratio_range)
    
    ii_brightness = integral_image_result(integral_image_brights, corner, neighbor);
    ii_darkness = integral_image_result(integral_image_darks, corner, neighbor);
        
    %target_diff   = abs(diff([corner(1), neighbor(1), corner(2), neighbor(2)]));
    %ii_max_value  = target_diff(1) * target_diff(3);
    %ii_darkness   = ii_max_value - ii_brightness;
    
    ii_bright_dark_ratio = ii_brightness / max(ii_darkness, 1e-15);
    
    result = ii_bright_dark_ratio >= brights_darks_ratio_range(1) && ...
             ii_bright_dark_ratio <= brights_darks_ratio_range(2);

end

function result = integral_image_result(integral_image, corner, neighbor)
    %{
disp("neighb");
    disp(neighbor);
    disp("corner");
    disp(corner);
%}
    ii_coords_1 = integral_image(corner(1), corner(2));
    ii_coords_12 = integral_image(corner(1), neighbor(2));
    ii_coords_2 = integral_image(neighbor(1), neighbor(2));
    ii_coords_21 = integral_image(neighbor(1), corner(2));
    
    % neighbor lies to the bottom left
    if corner(1) < neighbor(1) && corner(2) > neighbor(2)
        A = ii_coords_12;
        B = ii_coords_1;
        C = ii_coords_2;
        D = ii_coords_21;
        
    % neighbor lies to the bottom right
    elseif corner(1) < neighbor(1) && corner(2) < neighbor(2)
        A = ii_coords_1;
        B = ii_coords_12;
        C = ii_coords_21;
        D = ii_coords_2;
        
    else
        error('We assume the corners are sorted by Y-axis and neighbors have a Y greater than and X not equal to corner.');
    end
    
    result = D + A - B - C;
    
end

%{
    Returns an integral image (also called summed area table) from the
    given grayscale image (2-dimensional, [0, 1]-image).
    Sources:
        http://delivery.acm.org/10.1145/810000/808600/p207-crow.pdf
        accessed on 2019/11/12
    Author:
        Laurenz Edmund Fiala (11807869)
%}
function result = generate_integral_image(greyscale_image)
    
    result = greyscale_image;
    for x = 1:size(greyscale_image, 1)
       for y = 1:size(greyscale_image, 2)
           
           r = 0;
           if x > 1
               r = r + result(x-1, y);
           end
           
           if y > 1 
               r = r + result(x, y-1);
           end
           
           if x > 1 && y > 1
               r = r - result(x-1, y-1);
           end
           
           result(x, y) = r + greyscale_image(x, y);
           
       end
    end

end

function result = harris_corners(greyscale_image)
    % 2nd Step:
    %{    
     The sum of squared differences between 2 patches will be given 
     with following equation:
     (1) F(u,v) = sum([I(x+u,y+v) - I(x,y)]^2).
     Now using Taylor expansion we can approximate I(x+u,y+v) as following:
     I(x+u,y+v) = I(x,y) + dx * d(I(x,y))/dx + dy * d(I(x,y))/dy

     We can than rewrite the equation (1)  in the following way :
        F(u,v)*w(x,y) = sum(w(x,y)*[dx * d(I(x,y))/dx + dy * d(I(x,y))/dy]^2)
     Let call d(I(x,y))/dx = Ix and d(I(x,y))/dy = Iy
    (w(x,y) is window function, which can be gaussian or just a rectangular) 
    or in matrix form: F(u,v)*w(x,y)~=   [dx, dy] * M * [dx, dy]^T, where
    M = sum(w(x,y) [Ix^2 IxIy;IxIy Iy^2]) 
    From the PSA method we know, that the eigenvalues of this matrix would give us 
    the directions in which the data are mostly spread. 
    Thats why by analizing the eigenvalues of this matrix in each window, we 
    destinglish between corner, edge or none of both. 
%}
   
    %Spatial derivative calculation (Ix,Iy)
    % Create sobel operator in horizontal direction: 
    fx = [-1 0 1; -2 0 2; -1 0 1];
    % Apply it to the image
    Ix = filter2(fx,greyscale_image);
    % Create sobel operator in vertical direction:
    fy = [1 2 1; 0 0 0; -1 -2 -1];
    % Apply it to the image
    Iy = filter2(fy,greyscale_image);
    % We need to calculate also Ix^2, Iy^2, Ixy
    Ix2 = Ix.^2;
    Iy2 = Iy.^2;
    Ixy = Ix.*Iy;

    %3rd Step:  
    % as window function we choose gaussian. In the next step, we apply it to
    % the result
    h= fspecial('gaussian',[5 5],2); 
    Ix2 = filter2(h,Ix2);
    Iy2 = filter2(h,Iy2);
    Ixy = filter2(h,Ixy);

    % Eigenvalues are very expensive to calculate, thats why we calculate a
    % parameter R(i,j) = det(M)-0.06*(trace(M))^2;
    % where 0,06 is a chosen coefficient between [0.04, 0.06]
    % Now we can destinglish between different R
    % create a matrix R of a size of the Image
    R = zeros(size(greyscale_image,1),size(greyscale_image,2));

    % set Rmax to 0 first
    Rmax = 0; 
    for i = 1:size(greyscale_image,1)
        for j = 1:size(greyscale_image,2)
            M = [Ix2(i,j) Ixy(i,j);Ixy(i,j) Iy2(i,j)]; 
            R(i,j) = det(M)-0.04*(trace(M))^2;
            if R(i,j) > Rmax
                 Rmax = R(i,j);
            end
        end
    end
    % searching for local maxima as corner within the window of 3/3 
    result = zeros(size(greyscale_image,1),size(greyscale_image,2)); 
    for i = 3:size(greyscale_image,1)-2 % height 
         for j = 3:size(greyscale_image,2)-2 %width
                if R(i,j) > 0.0006*Rmax && ...
                    R(i,j)< 0.5*Rmax && ...
                    R(i,j) > R(i-1,j-1) && ...
                    R(i,j) > R(i-1,j) && ...
                    R(i,j) > R(i-1,j+1) &&...
                    R(i,j) > R(i,j-1) &&...
                    R(i,j) > R(i,j+1) && ...
                    R(i,j) > R(i+1,j-1) &&...
                    R(i,j) > R(i+1,j) && ...
                    R(i,j) > R(i+1,j+1) && ...
                    R(i,j) > R(i+2,j+2)&& ...
                    R(i,j) > R(i+2,j-2) &&...
                    R(i,j) > R(i+2,j) && ...
                    R(i,j) > R(i,j+2) && ...
                    R(i,j) > R(i,j-2) &&...
                    R(i,j) > R(i-2,j+2) &&...
                    R(i,j) > R(i-2,j-2) && ...
                    R(i,j) > R(i-2,j) 
                    result(i,j) = 1;
                end
         end
    end
    [posc, posr] = find(result == 1);
    % transporieren gleich hier
    result = [posc, posr];
end
