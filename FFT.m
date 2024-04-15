%%%%%%%%%%%%%%%%%%%%%%%%%% BREAK %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
clear all
%%%%%%%%%%%
%binsize = 4.0478;
binsize = 10;
T = 300.0;

%input = sprintf('MgAl-Screw%.1f-%.2f',T,binsize);   % without smooth
%input = sprintf('MgAl-Screw%.1f-%.2f-SMOOTH',T,binsize); % smooth 1
input = sprintf('MgAl-Screw%.1f-%.2f-2',T,binsize);   % smooth 2

loadname = sprintf('%s.mat',input);
savename = sprintf('%s-FFT-Smooth.mat',input);
fileID = fopen('MgAl-Screw-MIN-K-Cn.txt','w');

load (loadname)
clip = 0;
num_interp = 800;

kb = 1.38*10^-23; % boltzmann's constant in J/K
kb = kb*1000; % mJ/K
%kb = 0.00008617;
N = length(ix);
% doing interpolation
interp_x = linspace(min(ix), max(ix), num_interp);

%partials
for n = 1:length(i11(:,1))
    ii1(n,:) = interp1(ix,i11(n,:),interp_x,'pchip');   %'pchip' stands for Piecewise Cubic Hermite Interpolation
    ii2(n,:) = interp1(ix,i12(n,:),interp_x,'pchip');
    ii3(n,:) = interp1(ix,i21(n,:),interp_x,'pchip');
    ii4(n,:) = interp1(ix,i22(n,:),interp_x,'pchip');   
end
%full
for n = 1:length(dPy1(:,1))
    dPy1ii1(n,:) = interp1(ix,dPy1(n,:),interp_x,'pchip');   %'pchip' stands for Piecewise Cubic Hermite Interpolation
    dPy2ii2(n,:) = interp1(ix,dPy2(n,:),interp_x,'pchip');  
end

clear x i1 i2 i3 i4 i11 i12 i21 i22 dPy1 dPy2 ix
x = interp_x;
i1 = ii1;
i2 = ii2;
i3 = ii3;
i4 = ii4;
i5 = dPy1ii1;
i6 = dPy2ii2;
clear interp_x ii1 ii2 ii3 ii4 dPy1ii1 dPy2ii2

%%%%%%%
figure(1)
subplot(2,2,1)
plot(x,i1)
axis equal
subplot(2,2,2)
plot(x,i2)
axis equal
subplot(2,2,3)
plot(x,i3)
axis equal
subplot(2,2,4)
plot(x,i4)
axis equal
%%%%%%
    
lx = max(x) - min(x)
%lx = 495.5174
%lx = 400
runs = length(i1(:,1));
fruns = length(i5(:,1));
%%%%%%%%%%% partial
for n = 1:runs
    sum1 = 0;
    sum2 = 0;
    sum3 = 0;
    sum4 = 0;
    for nn = 1:length(x)-1
        dist1 = sqrt((x(nn)-x(nn+1))^2+(i1(n,nn)-i1(n,nn+1))^2);    % distance
        sum1 = sum1 + dist1;
        
        dist2 = sqrt((x(nn)-x(nn+1))^2+(i2(n,nn)-i2(n,nn+1))^2);
        sum2 = sum2 + dist2;
        
        dist3 = sqrt((x(nn)-x(nn+1))^2+(i3(n,nn)-i3(n,nn+1))^2);
        sum3 = sum3 + dist3;
        
        dist4 = sqrt((x(nn)-x(nn+1))^2+(i4(n,nn)-i4(n,nn+1))^2);
        sum4 = sum4 + dist4;
    end
    lx1(n,1) = sum1;
    lx2(n,1) = sum2;
    lx3(n,1) = sum3;
    lx4(n,1) = sum4;
end
%%%%%%%%%%% full
for n = 1:fruns
    sum5 = 0;
    sum6 = 0;
    for nn = 1:length(x)-1
        dist5 = sqrt((x(nn)-x(nn+1))^2+(i5(n,nn)-i5(n,nn+1))^2);    % distance
        sum5 = sum5 + dist5;
        
        dist6 = sqrt((x(nn)-x(nn+1))^2+(i6(n,nn)-i6(n,nn+1))^2);
        sum6 = sum6 + dist6;
    end
    lx5(n,1) = sum5;
    lx6(n,1) = sum6;
end

%partial
l1 = mean(lx1);
l2 = mean(lx2);
l3 = mean(lx3);
l4 = mean(lx4);
%full
l5 = mean(lx5);
l6 = mean(lx6);

lx = mean([l1 l2 l3 l4]);    % mean of the length scales partial
lxf = mean([l5 l6]);    % mean of the length scales full

K = [1:num_interp] * 2 * pi / lx;   % wave number calculation partial
Kf = [1:num_interp] * 2 * pi / lxf;   % wave number calculation full

runs = length(i1(:,1));
fruns = length(i5(:,1));
for run = 1:runs
    zi1(run,:) = i1(run,:) - mean(i1(run,:));
    zi2(run,:) = i2(run,:) - mean(i2(run,:));
    zi3(run,:) = i3(run,:) - mean(i3(run,:));
    zi4(run,:) = i4(run,:) - mean(i4(run,:));
    
    zi1((isnan(zi1))) = 0 ;
    zi2((isnan(zi2))) = 0 ;
    zi3((isnan(zi3))) = 0 ;
    zi4((isnan(zi4))) = 0 ;
    
    fft_zi1(run,:) = (ifft(zi1(run,:)));
    fft_zi2(run,:) = (ifft(zi2(run,:)));
    fft_zi3(run,:) = (ifft(zi3(run,:)));
    fft_zi4(run,:) = (ifft(zi4(run,:)));
end

for run = 1:fruns
    zi5(run,:) = i5(run,:) - mean(i5(run,:));
    zi6(run,:) = i6(run,:) - mean(i6(run,:));
    
    zi5((isnan(zi5))) = 0 ;
    zi6((isnan(zi6))) = 0 ;

    fft_zi5(run,:) = (ifft(zi5(run,:)));
    fft_zi6(run,:) = (ifft(zi6(run,:)));
end

figure(2)
hold on
plot(x,zi1)
plot(x,zi2)
plot(x,zi3)
plot(x,zi4)
hold off
axis equal

for run =1:runs
    Ak1S(run,:) = abs(fft_zi1(run,:)).^2;
    Ak2S(run,:) = abs(fft_zi2(run,:)).^2;
    Ak3S(run,:) = abs(fft_zi3(run,:)).^2;
    Ak4S(run,:) = abs(fft_zi4(run,:)).^2;
end

for run =1:fruns
    Ak5S(run,:) = abs(fft_zi5(run,:)).^2;
    Ak6S(run,:) = abs(fft_zi6(run,:)).^2;
end

%Averaging Power Spectra
for run=1:num_interp
    Ak1(run) = mean(Ak1S(:,run));
    Ak2(run) = mean(Ak2S(:,run));
    Ak3(run) = mean(Ak3S(:,run));
    Ak4(run) = mean(Ak4S(:,run));
    Ak5(run) = mean(Ak5S(:,run));
    Ak6(run) = mean(Ak6S(:,run));
end

Ak1 = Ak1(1:floor(num_interp/2));
Ak2 = Ak2(1:floor(num_interp/2));
Ak3 = Ak3(1:floor(num_interp/2));
Ak4 = Ak4(1:floor(num_interp/2));
%Ak5 = Ak5(1:floor(num_interp/2));
%Ak6 = Ak6(1:floor(num_interp/2));

K = K(1:floor(num_interp/2));
Kf = Kf(1:floor(num_interp));
 %kaxis = K.^2;
 %y1axis = (kb.*T)./(lx.*Ak1);
 %y2axis = (kb.*T)./(lx.*Ak2);
 %y3axis = (kb.*T)./(lx.*Ak3);
 %y4axis = (kb.*T)./(lx.*Ak4);

kaxis = K;
kaxisf = Kf;

y1axis = (Ak1);
y2axis = (Ak2);
y3axis = (Ak3);
y4axis = (Ak4);
y5axis = (Ak5);
y6axis = (Ak6);

kaxis = kaxis(2:end);
kaxisf = kaxisf(2:end);

y1axis = y1axis(2:end);
y2axis = y2axis(2:end);
y3axis = y3axis(2:end);
y4axis = y4axis(2:end);
y5axis = y5axis(2:end);
y6axis = y6axis(2:end);

myaxis = (y1axis+y2axis+y3axis+y4axis)/4;
myaxisf = (y5axis+y6axis)/2;

% Define the data points from french paper for screw dislocation
x_screw = [0.012479914425313888, 0.025577888, 0.037642450053088757, 0.049997632193936846, ...
    0.062842666, 0.075928241, 0.088897718, 0.10084655145164367, 0.11350439153143622, ...
    0.1257493827117733, 0.1393153783658891, 0.15193871285895555, 0.1643934071579327, ...
    0.1879695515765981, 0.1879695515765981, 0.20019760798462424, 0.2166185277479304, ...
    0.2556036764787999, 0.2556036764787999, 0.274400555, 0.3040080631649444, ...
    0.3341585299313714, 0.36731092194962767, 0.39123691388552095, 0.43347165531460036, ...
    0.48406135380797327, 0.5279472808152246, 0.5849401567480439, 0.6481061734951398, ...
    0.7294421431249125, 0.801836798, 0.8608443503013042, 1.0237669744962554, ...
    1.0237669744962554, 1.098966365355237, 0.9388003776656697];

y_screw = [0.36965846769284505, 0.1216102227944868, 0.06683199, 0.044832229, 0.032297696, ...
    0.026445272125948684, 0.021051244, 0.018250932923396752, 0.015823897117818238, ...
    0.014317942725018795, 0.012955309450671169, 0.011395307038857102, 0.010761011, ...
    0.009067931, 0.009067931, 0.008686462, 0.00786092, 0.006717532, 0.006717532, ...
    0.006079407, 0.00542329, 0.004977561, 0.004440576, 0.003962291, 0.003387282, ...
    0.002895578, 0.002440714, 0.002086517, 0.001733787, 0.001524719, 0.001322152, ...
    0.001146664, 0.001008051, 0.001008051, 0.000979404, 0.001052473];
figure(3)
loglog(x_screw, y_screw, 'o-', 'DisplayName','Geslin & Rodney 2018');
hold on;
%plot(kaxis,myaxis,'o')
loglog(kaxis,myaxis,'s', 'DisplayName','NVE bin=2.4 int=400 NY=2 fr=0.01');
loglog(kaxisf,myaxisf,'s', 'DisplayName','NVE bin=2.4 int=400 NY=2 fr=0.01 f');
%loglog(myaxis,kaxis,'s')
%hold off
xlim([0.01 2.0])
ylim([0.0007 1.5])
xlabel('k_n [A^{-1}]') 
%ylabel('<|C_n|^2> [A^2]') 
legend('show', 'Location', 'Northeast');

for kk =1:length(kaxis) 
    fprintf(fileID,'%f %e\n',kaxis(kk),myaxis(kk));
end

save(savename,'kaxis','y1axis', 'y2axis', 'y3axis','y4axis', 'y5axis', 'y5axis',  'myaxis', 'myaxisf')



