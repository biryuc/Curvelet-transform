close all;
clear all;

dx = 50; % Шаг по х
dt = 0.002; % Шаг по времени

agc_data = zeros(400,50);
NT = size(agc_data,1); %Измерения по времени
NR = size(agc_data,2); %Измерений по расстоянию


%Задаем 4 прямые y = ax + b
a_(1) =0.1;
a_(2) = -0.15;
a_(3) = 0.05;
a_(4) = 0.3;

b_(1) = 0.2;
b_(2) = 0.3;
b_(3) = 0.15;
b_(4) = 0.4;

%Формируем сейсмограмму
for nn = 1:length(a_)
    for i = 1:NR  
        delta_t = i*a_(nn)/NR + b_(nn); 
        for j = 1:NT
            t = j*dt;
            r=pi*40.0*(t - delta_t);
            agc_data(j,i) = agc_data(j,i)+(1-2*r.*r).*exp(-r.*r); 
        end
    end
end


%Прореживаем
[ S,data_sparse,x_offset ] = Sparse(agc_data,NT,NR,dx);
% Полноразмерная sparse_data
sparse_data_full_size = S.*agc_data;
% data_sparse(:,1) =  agc_data(:,1);
% n = 1;
% k_j(1) = 1;
% %Sampling matrix
% S = zeros(NT,NR);
% for j = 2:(NR-1) 
%     k = mod(j,2);
%     if k == 0 
%         n = n+1;
%         k_j(n) = j;
%         data_sparse(:,n) = agc_data(:,j);
%         S(:,n) = ones(NT,1);
%     end
% end
% n = n+1;
% k_j(n) = NR;
% data_sparse(:,n) = agc_data(:,NR);
% 
% x_offset = k_j.*dx;



%Отрисовка 
subplot(1,5,1)
title('Original data');
wigb(agc_data,40,dx:dx:dx*NR,dt:dt:dt*NT);
subplot(1,5,2)
title('Sampling data');
wigb(S,40,dx:dx:dx*NR,dt:dt:dt*NT);
subplot(1,5,3)
title('Sparse data');
wigb(data_sparse,40,x_offset,dt:dt:dt*NT);
% wigb(sparse_data_full_size,40,dx:dx:dx*NR,dt:dt:dt*NT);


%Добавление шума
NT_sparse = size(data_sparse,1);
NR_sparse = size(data_sparse,2);
sigma = 0.1;
% noisy_data = zeros(NT_sparse,NR_sparse);
% Rand_noise_sigma = sigma*randn(NT_sparse,NR_sparse);
% 
% noisy_data = data_sparse + Rand_noise_sigma;

noisy_data = zeros(NT,NR);
Rand_noise_sigma = sigma*randn(NT,NR);

noisy_data_full_size = sparse_data_full_size + Rand_noise_sigma;

%Отрисовка
subplot(1,5,4)
title('Sparse + noise data');
% wigb(noisy_data,40,x_offset,dt:dt:dt*NT);
wigb(noisy_data_full_size,40,dx:dx:dx*NR,dt:dt:dt*NT);


%%%%%%%%%%%%%% Алгоритм %%%%%%%%%%%%%%%
d = zeros(NT,NR);%initial solution d_0
% d_obs = noisy_data;
d_obs = noisy_data_full_size;
N = 10; %iterative number


% Curvelet преобразование подготовленных данных
N_cells = 4;
% Curve_data = fdct_wrapping(noisy_data,1,1,N_cells);
Curve_data = fdct_wrapping(noisy_data_full_size,1,1,N_cells);
%Поиск tau_0
%Нормировка
L_norm = 2; %Какую надо?
for i = 1: N_cells
    for j = 1: length(Curve_data{1, i})
        norm_cell_data(i,j) =  norm(Curve_data{1, N_cells}{1,j},L_norm); 
    end
end

% Точно ли так?
tau_0 = max(norm_cell_data,[],"all");
tau_min = min(norm_cell_data,[],"all");

k = 0;
%treshhold
trh = sqrt(tau_0);

%Алгоритм ALFT по шагам
L_2_norm = 2;
% y_gamma = 
for k = 0: N
   %step 1 calc prediction
   r = d_obs - S.*d;
   %доп условие выхода
%    r_2norm = norm(r,2);
%    if r_2norm< y_gamma
%        break;
%    end
   
   %step2 
   k_ = rand();
   arg_data = d + r.*k_;
   Curve_data_step = fdct_wrapping(arg_data,1,1,N_cells);
  
   %soft_tresh
   for n = 1: N_cells
        for j = 1: length(Curve_data{1, n})
%             Th = ones(size(Curve_data_step{1, N_cells}{1,j},1), size(Curve_data_step{1, N_cells}{1,j},2))
%             Th = Th*trh;
%             Curve_data_soft_tresh(n,j) =  Curve_data_step{1, N_cells}{1,j} - Th; 
            %Отнимаем trh из внутренней матрицы
            for l = 1:size(Curve_data_step{1, N_cells}{1,j},1)
                 for p = 1:size(Curve_data_step{1, N_cells}{1,j},2)
                        if abs(Curve_data_step{1, N_cells}{1,j}(l,p)) < trh
                            Curve_data_step{1, N_cells}{1,j}(l,p) = 0;
                        elseif Curve_data_step{1, N_cells}{1,j}(l,p) > trh
                            Curve_data_step{1, N_cells}{1,j}(l,p) = Curve_data_step{1, N_cells}{1,j}(l,p) - trh;
                        elseif Curve_data_step{1, N_cells}{1,j}(l,p) < -trh
                            Curve_data_step{1, N_cells}{1,j}(l,p) = Curve_data_step{1, N_cells}{1,j}(l,p) + trh; 
                        end
                 end
             end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        end
   end%soft tresh
   Alpha = Curve_data_step;
   
   %Reverse transform
   d = real(ifdct_wrapping(Alpha,1,NT,NR));
   
   %step3
   %update treshhold
   c = log(tau_0/tau_min);
   trh = tau_0*exp(c*(k-1)*(N-1));
   


end


%Отрисовка
subplot(1,5,5)
title('ALFT');
wigb(d,40,dx:dx:dx*NR,dt:dt:dt*NT);

orig_vs_Noisy = 20*log10(norm(agc_data(:))/norm(abs(agc_data(:)-noisy_data_full_size(:))))
orig_vs_Denoised = 20*log10(norm(agc_data(:))/norm(abs(agc_data(:)-d(:))))


function [ S,data_sparse,x_offset ] = Sparse(agc_data,NT,NR,dx)
    data_sparse(:,1) =  agc_data(:,1);
    n = 1;
    k_j(1) = 1;
    %Sampling matrix
    S = zeros(NT,NR);
    for j = 2:(NR-1) 
        k = mod(j,2);
        if k == 0 
            n = n+1;
            k_j(n) = j;
            data_sparse(:,n) = agc_data(:,j);
            S(:,j) = ones(NT,1);
        end
    end
    n = n+1;
    k_j(n) = NR;
    data_sparse(:,n) = agc_data(:,NR);

    x_offset = k_j.*dx;

end













