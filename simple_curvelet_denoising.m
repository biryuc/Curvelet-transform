% [filename,pathname] = uigetfile('*.*','Select the input file');
% 
% filewithpath = strcat(pathname,filename);
% 
% img = imread(filewithpath);
% 
% imgd = double(img);

[imgd,dx,dt,NR,NT] = SYNTH_REFL_SPARS();

full_size = size(imgd,1)*size(imgd,2);

n = size(imgd,1);
m =  size(imgd,2);


%Добавляем шум 
a =  abs(min( imgd ,[], "all" ));
b =  max( imgd ,[], "all" );
r = a + (b-a).*rand(1,1);

sigma = 0.1;
noisy_img = zeros(n,m);
Rand_noise = randn(n,m);
Rand_noise_sigma = sigma*randn(n,m);
noisy_img = imgd + Rand_noise_sigma;
% for s = 1:n
%     for w = 1:m
%        r = a + (b-a).*rand(1,1);
%        noisy_img(s,w) =  imgd(s,w) + sigma*r; 
%     end
% end


F = ones(n,m);

X = fftshift(ifft2(F)) * sqrt(numel(F));

Cn = fdct_wrapping(X,0,1,5);

E = cell(size(Cn));


for s = 1:length(Cn)
    E{s} = cell(size(Cn{s}));
    for w = 1:length(Cn{s})
        A = Cn{s}{w}; 
        E{s}{w} = sqrt(sum(sum(A.*conj(A))) / numel(A));
    end
end

C = fdct_wrapping(noisy_img,1,1,5);
Ct = C;

% t01 = max(C{1,1}{1,1},[],'all');
% t02_max = zeros(16,1);
% for s = 1:16
%     t02_max(s) = max(C{1,2}{1,s},[],'all');
% end
% t02 = max(t02_max);
% t03 = max(C{1,3}{1,1},[],'all');
% 
% t0 = [t01,t02,t03];
for s = 1:length(C)  
    tresh = 3*sigma + sigma*(s == length(C));
    for w = 1:length(C{s})
        Ct{s}{w} = C{s}{w}.* (abs(C{s}{w}) > tresh*E{s}{w});
    end
end


restored_image = real(ifdct_wrapping(Ct,1,n,m));

subplot(131);
wigb(imgd,60,dx:dx:dx*NR,dt:dt:dt*NT);


subplot(132);
wigb(noisy_img,60,dx:dx:dx*NR,dt:dt:dt*NT);

subplot(133);
wigb(restored_image,60,dx:dx:dx*NR,dt:dt:dt*NT);


orig_vs_Noisy = 20*log10(norm(imgd(:))/norm(imgd(:)-noisy_img(:)))
orig_vs_Denoised = 20*log10(norm(imgd(:))/norm(imgd(:)-restored_image(:)))
    
    
