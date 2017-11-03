x = magic(8)
w = hamming(8)
fft2(x.*(w*w'))
y1 = fft2(x.*(w*w'));
y2 = fft(bsxfun(@times,fft(bsxfun(@times,x,w)),w'),[],2);
imagesc(x)
max(abs(y1(:)-y2(:)))

figure
imagesc(ifft2(y1))