%spec_struct

clear; clc;
% Initialize fusome/ring system

x_voxel = 0.1515;
z_voxel = 0.2098;
voxelSize = [x_voxel x_voxel z_voxel];

%input sample number and number of cells connected by fusome piece
image_num = '4';
num_cells = '13';
sample = '1';
smooth = 5;
uiopen(strcat('C:\Users\rockyd\Desktop\Manuscripts in Progress\Male_Fusome_Structure_Science\Male_Fusome_sqh_add\Trained\Male_Fusome_Test_',image_num,'_',num_cells,'cell_',sample,'_Probabilities.h5'),1);

% preparing probabilities:
image_germ = double(I.img);
image_germ = permute(image_germ, [1 2 4 3]);
image_germ = image_germ/max(image_germ(:));
ring = image_germ(:,:,:,2);
fusome = image_germ(:,:,:,1);

% Clean up rings and segment fusome
bw_fusome = fusome > 0.5 | ring > 0.5;
frame = true(size(bw_fusome));
frame(2:end-1, 2:end-1, :) = false;
labels = bwlabeln(bw_fusome);
to_remove = integer_unique(labels(frame));
labels(ismember(labels, to_remove)) = 0;
labels(~isolate_lcc(labels)) = 0;
newfus = isolate_lcc(labels > 0 & fusome > 0.5);
newring = bwlabeln(labels > 0 & ring > 0.5);

%remove small objects that are not rings
s = bwskel(newring > 0);
change_occur = true;
while change_occur
    new_s = s;
    new_s(bwmorph3(s,'endpoints')) = false;
    change_occur = ~isequal(new_s, s);
    s = new_s;
end

% finish removing things that are not rings...
s = bwmorph3(s, 'clean');
rings_to_keep = integer_unique(nonzeros(newring(s)));
newring(~ismember(newring, rings_to_keep)) = 0;
newring = bwlabeln(newring > 0);
num_rings = length(unique(newring))-1;
fus_parts = length(unique(newfus))-1;
[x,y,z] = ind2sub(size(s), find(s));
x = x * voxelSize(1);
y = y * voxelSize(2);
z = z * voxelSize(3);

%expand out size to make isotropic directions
in_x = (0:size(newfus,1)-1)*voxelSize(1);
in_y = (0:size(newfus,2)-1)*voxelSize(2);
in_z = (0:size(newfus,3)-1)*voxelSize(3);

out_x = (0:size(newfus,1)-1)*voxelSize(1);
out_y = (0:size(newfus,2)-1)*voxelSize(1);
out_z = (0:size(newfus,3)-1)*voxelSize(1);

interp_newfus = medfilt3(interp3(single(newfus), 1:size(newfus,2), (1:size(newfus,1))', 1 : voxelSize(1)/voxelSize(3) : size(newfus,3), 'linear') >= 0.5);
interp_newring = medfilt3(bwlabeln(interp3(single(newring > 0), 1:size(newring,2), (1:size(newring,1))', 1 : voxelSize(1)/voxelSize(3) : size(newfus,3), 'linear') >= 0.5));
%% visualization before segmenting

for i = 1:num_rings
    p(i) = isosurface(interp_newring == i);
    p(i).vertices(:,1:2) = p(i).vertices(:,1:2)*voxelSize(1);
    p(i).vertices(:,3) = p(i).vertices(:,3)*voxelSize(3)*voxelSize(1)/voxelSize(3);
end

figure;
for j = 1:length(p)
    patch(p(j),'FaceColor','r','EdgeColor','none'); axis image; alpha(1)
end

p1 = isosurface(interp_newfus == 1);
p1.vertices(:,1:2) = p1.vertices(:,1:2)*voxelSize(1);
p1.vertices(:,3) = p1.vertices(:,3)*voxelSize(3)*voxelSize(1)/voxelSize(3);
patch(p1,'FaceColor',[0.7 0.7 0.7],'EdgeColor','none'); axis image; alpha(1)
box off; grid off; axis off;
view([-38 -59])
%% split fusome and rings into component parts
clc;
newfus_labeled = split_fusome(0, interp_newfus, interp_newring, smooth);

%% Find adjacencies between objects
bw = imdilate(interp_newring, create_ball(1));
newimage = newfus_labeled.*(bw > 0);

if exist('list','var')
    clear list
end

for i = 1:(num_rings)
    check_obj = bwconncomp(bw == i,26);
    props = regionprops(check_obj,'Area','BoundingBox','Centroid','PixelIdxList','PixelList');
    list{i} = props.PixelIdxList;
end

matchup = zeros(length(list),2);

for j = 1:length(list)
    output = newimage(list{j});
    
    number = zeros(1,num_rings);
    for k = 1:num_rings+1
        number(k) = sum(output == k);
    end
    
    if nnz(number) < 2
        matchup(j,1) = 0;
        matchup(j,2) = 0;
    elseif nnz(find(number == max(number))) == 2
        [indices] = find(number == max(number));
        matchup(j,1) = indices(1);
        matchup(j,2) = indices(2);
    elseif max(number) == max(number(number < max(number)))
        matchup(j,1) = find(number == max(number));
        matchup(j,2) = find(number == max(number));
    else
        matchup(j,1) = find(number == max(number));
        matchup(j,2) = find(number == max(number(number<max(number))));
    end
end

% Remove zero rows and rows with one element
rule = matchup(:,2) == 0;
matchup(rule,:) = [];

% Remove zero columns
matchup(:,all(~matchup,1)) = [];

matchupswitch = [matchup(:,2) matchup(:,1)];
matchup = unique([matchup; matchupswitch], 'rows');

%Remove duplicate pairs
pair = matchup(:,1) < matchup(:,2);
matchup = matchup.*pair;
rule = matchup(:,2) == 0;
matchup(rule,:) = [];

total_seg = max(matchup(:));

%get fusome volumes for each portion
fus_vol_out = zeros(total_seg,1);
for i = 1:total_seg
    fus_vol_out(i) = nnz(newfus_labeled(:) == i).*voxelSize(1).*voxelSize(2).*voxelSize(3).*voxelSize(1)./voxelSize(3);
end

matchupmatrix = zeros(total_seg);
for i = 1:size(matchup,1)
    matchupmatrix(matchup(i,1),matchup(i,2)) = 1;
    matchupmatrix(matchup(i,2),matchup(i,1)) = 1;
end

%relabel fusome pieces
output_fus = zeros(size(newfus_labeled));
fusome_vol = zeros(size(fus_vol_out));
for m = 1:(num_rings+1)
    output_fus(newfus_labeled == m) = m;
    fusome_vol(m) = fus_vol_out(m);
end