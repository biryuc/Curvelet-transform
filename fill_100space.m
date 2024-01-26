close all;
clear all;
load('data_curvelet_after_reg.mat');
load('curveletcolormap.mat');
dx = 10; % Шаг по х
dt = 0.002; % Шаг по времени
Ampl = 40.0;
N_steps = 200;
N_cells = 5;

 NR_full = size(data_curvelet_reg,2);
 NT_full = size(data_curvelet_reg,1);
SPace = zeros(NT_full,101);



data_curvelet_reg = AGCgain(data_curvelet_reg,dt,0.5,1);
data_curvelet_reg(1:NT_full,900:1000) = SPace(:,:);
figure(1)
imagesc(data_curvelet_reg);
colormap(CustomColormap);
NR_size = 25;


init_NT = 1;
init_NR = 900;
init_NR_left_bound = 900 - NR_size ;
init_NR_right_bound = 1000;
NT_size = 200;

for j = 1:7
    for i = 1:2
        

        if i == 1
            clear one_part;
            one_part(:,:) = zeros(NT_size,3*NR_size);
            one_part(1:NT_size,1:2*NR_size) = data_curvelet_reg(init_NT: init_NT+ NT_size -1,init_NR_left_bound:init_NR_left_bound + 2*NR_size - 1);
            one_part(1:NT_size,2*NR_size:3*NR_size) = data_curvelet_reg(init_NT:init_NT+ NT_size - 1,init_NR_right_bound:init_NR_right_bound + NR_size);

        else
        clear one_part;
        one_part(:,:) = zeros(NT_size,3*NR_size + NR_size*(i-1));
        one_part(1:NT_size,1:i*NR_size) = data_reg(1: NT_size,1:i*NR_size); 
        one_part(1:NT_size,i*NR_size:(i+1)*NR_size-1) = SPace(1:NT_size,1:NR_size);
        one_part(1:NT_size,(i+1)*NR_size:(i+2)*NR_size) = data_curvelet_reg(init_NT:init_NT+ NT_size - 1,init_NR_right_bound:init_NR_right_bound + NR_size);
        end
        NR = size(one_part,2);
        NT = size(one_part,1);
        [one_part,coef,Stdev,Mean] = preprocess_data(one_part,NR,NT);

        [ S ] = Sparse_Curvelet(one_part,NT,NR);

        sparse_data_full_size = S.*one_part;

         Start_curvelet = cputime;
           [data_reg] = Curvelet_function(sparse_data_full_size,NR,NT,N_steps,N_cells,S);
    %      [data_reg_curvelet] = ALFT_function(data_sparse,NR,NT,NR_sparse,dx,N_steps,k_j);
         Elapsed_curvelet = cputime - Start_curvelet


         [data_reg] = reverse_ampl(data_reg,Mean,Stdev,NR,NT);

         data_curvelet_reg(init_NT: init_NT+NT_size -1,init_NR +(i-1)*NR_size:init_NR +(i)*NR_size) = data_reg(1: NT_size,i*NR_size:(i+1)*NR_size);  


         figure(2)
         imagesc(data_reg);
         colormap(CustomColormap);

          figure(3)
         imagesc(data_curvelet_reg);
         colormap(CustomColormap);

         init_NR_left_bound = init_NR_left_bound + NR_size;
%          init_NR_right_bound = init_NR_right_bound + NR_size;
    end     
    init_NR_left_bound = 900 - 50 ;
     init_NT = init_NT + 200;
end


figure(4)
imagesc(data_curvelet_reg);
colormap(CustomColormap);
