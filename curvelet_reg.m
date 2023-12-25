close all;
clear all;

% [filename,pathname] = uigetfile('*.*','Select the input file');
% 
% filewithpath = strcat(pathname,filename);
% M = csvread(filewithpath);
% 
% imgd = double(M);


dx = 50; % Шаг по х
dt = 0.002; % Шаг по времени
Ampl = 40.0;

agc_data = zeros(400,50);
NT = size(agc_data,1); %Измерения по времени
NR = size(agc_data,2); %Измерений по расстоянию


% NT_ = size(M,1); %Измерения по времени
% NR_ = size(M,2); %Измерений по расстоянию
% dx_ = dt;
% dt_ = dx;
% figure
% title('Синтетические данные');
% xlabel('Смещение,м') 
% ylabel('Время,с') 
% wigb(M,Ampl,dx_:dx_:dx_*NR_,dt_:dt_:dt_*NT_);

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
            r=pi*Ampl*(t - delta_t);
            agc_data(j,i) = agc_data(j,i)+(1-2*r.*r).*exp(-r.*r); 
        end
    end
end



% [imgd,dx,dt,NR,NT] = SYNTH_REFL_SPARS();
% 
% agc_data = imgd;

agc_max = max(agc_data,[],"all");
agc_min = min(agc_data,[],"all");

Stdev = std(agc_data,0,"all");
Mean = mean(agc_data,"all");

for i = 1: NT
    for j = 1: NR
       agc_data(i,j) =    (agc_data(i,j) - Mean)/Stdev;
    end
end
agc_max_norm = max(agc_data,[],"all");
agc_min_norm = min(agc_data,[],"all");

coef = agc_min_norm/agc_max_norm;

%Отрисовка 
subplot(1,5,1)
title('Синтетические данные');
xlabel('Смещение,м') 
ylabel('Время,с') 
wigb(agc_data,Ampl,dx:dx:dx*NR,dt:dt:dt*NT);

%Добавление шума
noisy_data_full_size = agc_data; % + 2*coef*randn(NT,NR);

subplot(1,5,2)
title('noise data');
xlabel('Смещение,м') 
ylabel('Время,с')
% wigb(noisy_data,40,x_offset,dt:dt:dt*NT);
wigb(noisy_data_full_size,Ampl,dx:dx:dx*NR,dt:dt:dt*NT);

%Прореживаем
[ S,data_sparse,x_offset ] = Sparse(noisy_data_full_size,NT,NR,dx);
NT_sparse = size(data_sparse,1);
NR_sparse = size(data_sparse,2);

% Полноразмерная sparse_data + noise
sparse_data_full_size = S.*noisy_data_full_size;
subplot(1,5,3)
title('Sparse + noise data');
xlabel('Смещение,м') 
ylabel('Время,с')
%wigb(sparse_data_full_size,40,dx:dx:dx*NR,dt:dt:dt*NT);
 wigb(data_sparse,Ampl,x_offset,dt:dt:dt*NT);

%%%%%%%%%%%%%% Алгоритм %%%%%%%%%%%%%%%
d = zeros(NT,NR);%initial solution d_0
d_obs = sparse_data_full_size;
N = 150; %iterative number

% Curvelet преобразование подготовленных данных
N_cells = 4;
is_real = 1;
type_core_function = 1; %//2-wavelet 1-curvelet
Curve_data = fdct_wrapping(d_obs,is_real,type_core_function,N_cells);
%Поиск tau_0
%Нормировка
for i = 1: N_cells
    size_cell = length(Curve_data{1, i});
    for j = 1: length(Curve_data{1, i})
        single_matrix = Curve_data{1, i}{1,j};
        max_cell_data(i,j) =  max(Curve_data{1, i}{1,j},[],"all"); 
        min_cell_data(i,j) =  min(Curve_data{1, i}{1,j},[],"all"); 
    end
end

tau_0 = max(max_cell_data,[],"all"); 
tau_min = min(min_cell_data,[],"all");

k = 0;
%treshhold
trh = sqrt(tau_0)/3;


%Алгоритм  по шагам
L_2_norm = 2;
% y_gamma = 
h = 0;
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
   Curve_data_step = fdct_wrapping(arg_data,is_real,type_core_function,N_cells);
  
   %soft_tresh
   for n = 1: N_cells
        for j = 1: length(Curve_data_step{1, n}) 
            %Зануляем элементы
            for l = 1:size(Curve_data_step{1, n}{1,j},1)
                 for p = 1:size(Curve_data_step{1, n}{1,j},2)
                        if abs(Curve_data_step{1, n}{1,j}(l,p)) < trh
                            Curve_data_step{1, n}{1,j}(l,p) = 0;
                          
                        end
                 end
             end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        end
   end%soft tresh
  
   Alpha = Curve_data_step;
   
   %Reverse transform
   d = ifdct_wrapping(Alpha,is_real,NT,NR);
   
   %step3   
end

d = real(d);
%Отрисовка
subplot(1,5,4)
title('Curvelet');
xlabel('Смещение,м') 
ylabel('Время,с')
wigb(d,Ampl,dx:dx:dx*NR,dt:dt:dt*NT);



orig_vs_Noisy1 = 20*log10(mean(abs(agc_data),"all")/mean(abs(agc_data-noisy_data_full_size),"all"));
orig_vs_reg1 = 20*log10(mean(abs(agc_data),"all")/mean(abs(agc_data-d),"all"));

orig_vs_Noisy = 20*log10(norm(abs(agc_data))/norm(abs(agc_data-noisy_data_full_size)))
orig_vs_reg = 20*log10(norm(abs(agc_data))/norm(abs(agc_data-d)))

%%%%%%%%%%%%%%%%%%%%%% Curvelet Denoising %%%%%%%%%%%%%%%%%%%%%%
% [ d_den ] = Curvelet_Denoising(d,NT,NR,N_cells,trh);
% orig_vs_denois = 20*log10(norm(agc_data)/norm(abs(agc_data-d_den)))
% %Отрисовка
% figure
% subplot(1,1,1)
% title('Curvelet Denoising');
% xlabel('Смещение,м') 
% ylabel('Время,с')
% wigb(d_den,Ampl,dx:dx:dx*NR,dt:dt:dt*NT);


function [ S,data_sparse,x_offset ] = Sparse(agc_data,NT,NR,dx)
    data_sparse(:,1) =  agc_data(:,1);
    n = 1;
    k_j(1) = 1;
    counter = 0;
    %Sampling matrix
    S = zeros(NT,NR);
    for j = 2:(NR-1) 
       % k = mod(j,2);
%         k = 0;
        nn = rand;
        if nn>0.5 
        %if k == 0 
            counter = counter + 1;
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
    counter

end


%%%%%%%%%%%%%%%%%%%%%% Curvelet Denoising %%%%%%%%%%%%%%%%%%%%%%
function [ d ] = Curvelet_Denoising(data,NT,NR,N_cells,tresh)
        F = ones(NR,NT);
       
        X = fftshift(ifft2(F)) * sqrt(numel(F));

        Cn = fdct_wrapping(X,0,1,N_cells);

        E = cell(size(Cn));

        for s = 1:length(Cn)
            E{s} = cell(size(Cn{s}));
            for w = 1:length(Cn{s})
                A = Cn{s}{w}; 
                E{s}{w} = sqrt(sum(sum(A.*conj(A))) / numel(A));
            end
        end

        C = fdct_wrapping(data,1,1,N_cells);
        Ct = C;

        for s = 1:length(C)  
            tresh = tresh + tresh*(s == length(C));
            for w = 1:length(C{s})
                Ct{s}{w} = C{s}{w}.* (abs(C{s}{w}) > tresh*E{s}{w});
            end
        end
        d = real(ifdct_wrapping(Ct,1,NT,NR));

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%













