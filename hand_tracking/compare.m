

if 1
clear
num_img = 6;
figure(4)
hold off
for img_ind = 1:num_img
    img_name = sprintf('%02d.jpg',img_ind);
    img = imread(img_name);
    r = img(:,:,1);
    g = img(:,:,2);
    b = img(:,:,3);
    rgbvals = [r(:),g(:),b(:)]; hsv = rgb2hsv(img);
    v = hsv(:,:,3); v = v(r>5);
    s = hsv(:,:,2); s = s(r>5);
    h = hsv(:,:,1); h = h(r>5); 
    sample = randsample(numel(h), 10000);
    hall(1+(img_ind-1)*numel(sample):img_ind*numel(sample)) = h(sample);
    sall(1+(img_ind-1)*numel(sample):img_ind*numel(sample)) = s(sample);
    vall(1+(img_ind-1)*numel(sample):img_ind*numel(sample)) = v(sample);
    figure(i+20)
    %scatter3(h(sample),s(sample),v(sample),ones(size(sample))*0.01, double(rgbvals(sample,:))/255.0);
    %scatter3(h(sample),s(sample),v(sample),ones(size(sample))*0.01, repmat(double([i, i, 0.0])/6.0, 10000, 1));
    hold on
end
end
%H = [hall', sall', vall'];
H = [hall', sall'];
centH = H - repmat(mean(H), size(H,1), 1);
HTH = centH'*centH;
[U, S, V] = svd(HTH)
Hinv = inv(HTH');
Hmean = mean(H)

if 1
for img_ind = 1:num_img
    img_ind
    img_name = sprintf('t%02d.jpg',img_ind);
    img = imread(img_name);
    %img = imresize(img, 0.1);
    r = img(:,:,1);
    g = img(:,:,2);
    b = img(:,:,3);
    rgbvals = [r(:),g(:),b(:)]; hsv = rgb2hsv(img);
    h = hsv(:,:,1);
    s = hsv(:,:,2);
    v = hsv(:,:,3);
    hsvvals = [h(:),s(:)];
    repmean = repmat(Hmean,numel(h), 1);
    diff = (hsvvals-repmean)';
    temp = (diff'*Hinv);
    dists = exp(-sqrt(sum(temp.*diff', 2))/0.02);
    dists = reshape(dists,size(h));
    figure(img_ind)
    imagesc(dists);
end
end
