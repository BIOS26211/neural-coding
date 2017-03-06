close all; 

% Creates a raster of one set of data, to visualize it. 

cell1 = load('cell_11l10.mat');

data1 = cell1.data;
dirs = cell1.dirs;

s=size(data1);

figure 
hold on 
xlabel('Time(ms)');ylabel('Direction');
title('Neural Spikes by Stimulus Direction');
for i = 1:s(1)
    for j= 1:s(2)-1
        for k = 1:s(3)
            if data1(i,j,k) == 1
                spkx = [i*2, i*2];
                spky = [dirs(j) - 5,dirs(j) + 5];
                line(spkx,spky,'color','k','LineWidth',1);
            else 
            end
        end
    end
end
hold off

% Load 4 more cells in order to make a 'word'

cell2 = load('cell_11l18.mat');
data2 = cell2.data;
cell3 = load('cell_12r08.mat');
data3 = cell3.data;
cell4 = load('cell_14l07.mat');
data4 = cell4.data;
cell5 = load('cell_14l13.mat');
data5 = cell5.data;

% The data needs to be in 8ms time bins, currently it is in 2
% Since we are looking at it binarily, any point must either be 0 or 1. Any
% value 2 or higher is removed. 

s2=size(data2);s3=size(data3);s4=size(data4);s5=size(data5);
n = min([s(3);s2(3);s3(3);s4(3);s5(3)]);

% It can only create words for the number of trials in the data set with
% the smallest number of trials (otherwise you lose that letter of the
% word)

data1n = zeros(s(1)/4, s(2));
data2n = zeros(s(1)/4, s(2));
data3n = zeros(s(1)/4, s(2));
data4n = zeros(s(1)/4, s(2));
data5n = zeros(s(1)/4, s(2));

for i = 1:s(1)/4
    for j = 1:s(2)
        for k = 1:n
            data1n(i,j,k) = data1((4*i),j,k) + data1(((4*i)-1),j,k) + data1(((4*i)-2),j,k);
            if data1n(i,j,k)>=2
                data1n(i,j,k)=1;
            end
            data2n(i,j,k) = data2((4*i),j,k) + data2(((4*i)-1),j,k) + data2(((4*i)-2),j,k);
            if data2n(i,j,k)>=2
                data2n(i,j,k)=1;
            end
            data3n(i,j,k) = data3((4*i),j,k) + data3(((4*i)-1),j,k) + data3(((4*i)-2),j,k);
            if data3n(i,j,k)>=2
                data3n(i,j,k)=1;
            end
            data4n(i,j,k) = data4((4*i),j,k) + data4(((4*i)-1),j,k) + data4(((4*i)-2),j,k);
            if data4n(i,j,k)>=2
                data4n(i,j,k)=1;
            end
            data5n(i,j,k) = data5((4*i),j,k) + data5(((4*i)-1),k) + data5(((4*i)-2),j,k);
            if data5n(i,j,k)>=2
                data5n(i,j,k)=1;
            end
        end
    end
end 

delta_I = I_sigAB - I_sigA - I_sigB
