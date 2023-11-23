close all;
clear all;

dx = 50;
dt = 0.002;


agc_data = zeros(400,50);
NR = size(agc_data,2);
NT = size(agc_data,1);

a_(1) =0.1;
a_(2) = -0.15;
a_(3) = 0.05;
a_(4) = 0.3;

b_(1) = 0.2;
b_(2) = 0.3;
b_(3) = 0.15;
b_(4) = 0.4;


for nn = 1:length(a_)
for i = 1:NR
% ---------------------------------------------------------    
    delta_t = i*a_(nn)/NR + b_(nn); 
    for j = 1:NT
        t = j*dt;
        r=pi*40.0*(t - delta_t);
        agc_data(j,i) = agc_data(j,i)+(1-2*r.*r).*exp(-r.*r); 
    end
end
end
subplot(1,3,1)

 wigb(agc_data,40,dx:dx:dx*NR,dt:dt:dt*NT);
  

N_DFT = 10*NR;
N_steps = 100;

Ks = 1/dx;
X_max = NR*dx;
dk = Ks/N_DFT;
K = dk:dk:Ks;


 for j = 1:NR
     data_fft_0(:,j) = fft(agc_data(:,j));
 end

 n = 1;
k_j(1) = 1;
data_sparse(:,1) =  agc_data(:,1);
 for j = 2:(NR-1)
   
    nn = rand;
    if nn>0.7 
        n = n+1;
        k_j(n) = j;
        data_sparse(:,n) = agc_data(:,j);
    end
 end

n = n+1;
k_j(n) = NR;
data_sparse(:,n) = agc_data(:,NR);


NR_sparse = length(k_j);

for k = 2:(NR_sparse-1)
    delta_x(k) = (k_j(k+1)-k_j(k-1))/2*dx;
end
 delta_x(1) = ((k_j(2) - k_j(1))*0.5+1)*dx;
 delta_x(NR_sparse) = ((k_j(NR_sparse) - k_j(NR_sparse-1))*0.5+1)*dx;


 x_offset = k_j.*dx;
subplot(1,3,2)
 wigb(data_sparse,40,x_offset,dt:dt:dt*NT);


 for j = 1:NR_sparse
     data_fft(:,j) = fft(data_sparse(:,j));
 end

delta_X = NR*dx;
G = zeros(NT,N_DFT );
data_fft_reg  = zeros(NT,NR);
x_j = k_j*dx;
tic;

for i = 1:400
    data_fft_u = data_fft(i,:);
    res = 1.0;
    nn = 1;
    while (res>0.2) & (nn<100)
           nn = nn+1;        
 X=1/delta_X*delta_x.*data_fft_u;
  F = nufft(X,x_j,K);
  [argvalue, argmax] = max(abs(F));    

      G(i,argmax) =  G(i,argmax) + F(argmax);
      A_ls  = F(argmax);


            for j = 1:NR_sparse
             data_fft_u(j) = data_fft_u(j)-F(argmax)*exp(2*pi*1i*K(argmax)*k_j(j)*dx);

            end      
 res = norm(data_fft_u)/norm(data_fft(i,:));
  hold off; 

    end

end


toc;
data_fft_reg =zeros(NT,NR);
for i = 1:NT
for j = 1:NR
    for k = 1:N_DFT
    data_fft_reg(i,j) = data_fft_reg(i,j) + G(i,k)*exp(2*pi*1i*K(k)*j*dx);
    end
end
end

 for j = 1:NR
     data_reg(:,j) = ifft(data_fft_reg(:,j));
 end


subplot(1,3,3)
 wigb(real(data_reg),40,dx:dx:dx*NR,dt:dt:dt*NT);
 
k0 = 100;
figure;
subplot(3,1,1);
plot(real(data_fft_reg(k0,:)))
subplot(3,1,2);
plot(real(data_fft_0(k0,:)))
hold on;
plot(x_offset/50,real(data_fft(k0,:)),'b*')
subplot(3,1,3);
plot(real(data_fft_0(k0,:))-real(data_fft_reg(k0,:)))
imagesc(abs(fft2(data_reg)))


% imagesc(abs(data_fft));
