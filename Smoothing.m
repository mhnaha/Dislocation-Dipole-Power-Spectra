%%%%%%%%%%%%%%%%%%%%%%%%%% BREAK %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
clear all
%%%%%%%%%%%
%binsize = 4.0478;
binsize = 2.5;
int = 400; % interpolation 
mbin = binsize/2; % moving average cutoff in angstrom
T = 300;
%%%%%%%%%%%
%input = sprintf('MgAl-Screw300.0-%.2f',binsize);
input = sprintf('MgAl-Screw%.1f-%.2f',T,binsize);

%input = sprintf('MgAl-Screw%.1f-%.2f',TEMP,binsize);

loadname = sprintf('%s.mat',input);
savename = sprintf('%s-2.mat',input);
load (loadname)
% plotting four partials with respect to x direction
figure(1)
subplot(2,2,1)
plot(ix,i11)
xlabel('ix') 
ylabel('i11')
axis equal

subplot(2,2,2)
plot(ix,i12)
xlabel('ix') 
ylabel('i12')
axis equal

subplot(2,2,3)
plot(ix,i21)
xlabel('ix') 
ylabel('i21')
axis equal

subplot(2,2,4)
plot(ix,i22)
xlabel('ix') 
ylabel('i22')
axis equal
%%%%%%%%%%%
runs = length(i11(:,1));
%oN = length(i11(1,:));
%dix = ix(2)-ix(1);


N = int;
iix = linspace(min(ix),max(ix),N);
diix = iix(2) - iix(1);   % diix is calculated as the difference between the second 
                          % and first elements of the array iix
MA = floor(mbin/diix);    % mbin=binsize/2 it seems to be a "moving average" smoothing method

%%%%%%%%%%%
for run = 1:runs
    ii11(run,:) = interp1(ix,i11(run,:),iix,'pchip'); % This line performs 1-D interpolation using the 'pchip' method 
                                                      % for the data in i11(run,:). 
    ii12(run,:) = interp1(ix,i12(run,:),iix,'pchip');
    ii21(run,:) = interp1(ix,i21(run,:),iix,'pchip');
    ii22(run,:) = interp1(ix,i22(run,:),iix,'pchip');
end
%%%%%%%%%%%

for run = 1:runs
    for n = 1:N
        if n <= MA
            si11(run,n) = mean(ii11(run,1:n+MA));
            si12(run,n) = mean(ii12(run,1:n+MA));
            si21(run,n) = mean(ii21(run,1:n+MA));
            si22(run,n) = mean(ii22(run,1:n+MA));
        else if n >= N-MA
                si11(run,n) = mean(ii11(run,n-MA:end));
                si12(run,n) = mean(ii12(run,n-MA:end));
                si21(run,n) = mean(ii21(run,n-MA:end));
                si22(run,n) = mean(ii22(run,n-MA:end));
            else
                si11(run,n) = mean(ii11(run,n-MA:n+MA));
                si12(run,n) = mean(ii12(run,n-MA:n+MA));
                si21(run,n) = mean(ii21(run,n-MA:n+MA));
                si22(run,n) = mean(ii22(run,n-MA:n+MA));
            end
        end
    end
end

%%%%%%%%%%% uninterpolate
for run = 1:runs
    sii11(run,:) = interp1(iix,si11(run,:),ix,'pchip');
    sii12(run,:) = interp1(iix,si12(run,:),ix,'pchip');
    sii21(run,:) = interp1(iix,si21(run,:),ix,'pchip');
    sii22(run,:) = interp1(iix,si22(run,:),ix,'pchip');
end

      
%%%%%%%%%%%

figure(2)
subplot(2,2,1)
plot(ix,sii11)
xlabel('ix') 
ylabel('sii11')
axis equal
title('After Smoothing')

subplot(2,2,2)
plot(ix,sii12)
xlabel('ix') 
ylabel('sii12')
axis equal
title('After Smoothing')

subplot(2,2,3)
plot(ix,sii21)
xlabel('ix') 
ylabel('sii21')
axis equal
title('After Smoothing')

subplot(2,2,4)
plot(ix,sii22)
xlabel('ix') 
ylabel('sii22')
axis equal
title('After Smoothing')
clearvars -except ix sii11 sii12 sii21 sii22 savename
i11 = sii11;
i12 = sii12;
i21 = sii21;
i22 = sii22;

clearvars sii11 sii12 sii21 sii22

save(savename,'ix','i11','i12','i21','i22')
