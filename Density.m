clc
clear all
close all

TEMP = 300.0;   % Temperature
output = sprintf('MgAl-Screw%.1f',TEMP);    % output file generation

%pcutoff = 6;
interp = 800; % number of cells along the y direction 

latparam = 4.7; % Lattice parameter (change if you want). Should not necessarily be the same as lattice paprameter

binsize = latparam;   
savename = sprintf('%s-%.2f.mat',output,binsize);
plotname = sprintf('%s-%.2f.tiff',output,binsize);

%write_video = 0;
%make_video = 0;

dump_start =1;% first input dump file frame number to be read  
dump_end = 8780;%  last input dump file frame number to be read 
%dump_end = 6044;%  last input dump file frame number to be read 

ddump = 1;  % read each dump file, or every two dump, so on and so forth 

%dev = latparam;
%peak = latparam;
iter = 5;

heightcut1 = -10;% change  for Y min
heightcut2 = 500;% change  for Y max
thickcut1 = -10;% change   for X min
thickcut2 = 1000;% change   for X max


dump_list = dump_start:ddump:dump_end;
runs = length(dump_list);
fruns = length(dump_list);
for run = 1:runs
    clearvars -except output binsize savename plotname ...
        dump_start dump_end ddump startx ...
        heightcut1 heightcut2 thickcut1 thickcut2 dump_list runs run ix ...
        lx nx ix i11 i12 i21 i22 LENGTH THICK HEIGHT interp pcutoff TEMP
    
    status = sprintf('%i of %i configurations',run,runs);
    disp(status)
    
    
    input = sprintf('O_T300F1S1.%i',dump_list(run)); % change
    %# ATOMS id type x y z c_csys c_cna vx vy vz c_3[1] c_3[2] c_3[3]

    %M = dlmread(filename,delimiter,R1,C1) starts reading at row offset R1 and column offset C1. 
    % For example, the offsets R1=0, C1=0 specify the first value in the file.
    %To specify row and column offsets without specifying a delimiter, 
    % use an empty character as a placeholder, for example, M = dlmread(filename,'',2,1).
    dump = dlmread(input,'',9,0);
    
    N = length(dump(:,1));   % length of the dump file
    
    % Reading the orientation from the OVITO
    sZ = dump(:,5); 
    sY = dump(:,4); 
    sX = dump(:,3); 
    
    % Switching the orientation of the X and Z 
    dump(:,3) = sZ; %Length % change
    dump(:,4) = sY; %Height % change
    dump(:,5) = sX; %Thick % change
    
    %reorient structure
    
    
    if run ==1
        %plot raw dump file
        figure(1)
        subplot(1,3,1)        
        plot(dump(:,4),dump(:,5),'.')
        
        xlabel('y') 
        ylabel('x')
        title('YX - raw dump file')
        ylim([60 ,180]);
               

        subplot(1,3,2)        
        plot(dump(:,3),dump(:,5),'.')
        
        xlabel('z') 
        ylabel('x') 
        title('ZX - raw dump file')
        ylim([60 ,180]);
        %axis equal
        
        
        
        subplot(1,3,3)       
        plot(dump(:,3),dump(:,4),'.')
        
        xlabel('z') 
        ylabel('y') 
        title('ZY - raw dump file')
        ylim([70 ,400]);
        
        %axis equal
        
    else end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %trim dump file we basically do not need this part because we are
    %already trimming in OVITO using Cluster analysis
    count = 1;
    for n = 1:N     % loops over dump files to apply the following condition  (heightcut for y and thickcut for x)
        if dump(n,4) > heightcut1 && dump(n,4) < heightcut2 && ...
                dump(n,5) > thickcut1 && dump(n,5) < thickcut2
            tdump(count,:) = dump(n,:);
            count = count + 1;
        else if dump(n,4) > heightcut1 && dump(n,4) < heightcut2 && ...
                    dump(n,5) < -thickcut1 && dump(n,5) > -thickcut2
                tdump(count,:) = dump(n,:);
                count = count + 1;
           else end
        end
    end
    
    indices = find(tdump(:,7)==1);
    tdump(indices,:) = [];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if run ==1
        %plot trim dump file
        figure(2)
        subplot(1,3,1)
        plot(tdump(:,5),tdump(:,4),'.')
        title('XY - trimed dump file')
        xlabel('x') 
        ylabel('y') 

        axis equal
        
        subplot(1,3,2)
        plot(tdump(:,5),tdump(:,3),'.')
        title('XZ - trimed dump file')
        xlabel('x') 
        ylabel('z') 

        axis equal
        
        subplot(1,3,3)
        plot(tdump(:,4),tdump(:,3),'.')
        title('YZ - trimed dump file')
        xlabel('y') 
        ylabel('z') 
        axis equal
    else end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    z = tdump(:,5);
    y = tdump(:,4);
    x = tdump(:,3);
    
    
    % Separate dislocations
    
    c1 = 1;
    c2 = 1;
    for n = 1:length(x)
        if y(n) > 119 && y(n) < 128;     %here we can change the values based on fluctuation amplitude 
            dis1(c1,1) = x(n);
            dis1(c1,2) = y(n);
            dis1(c1,3) = z(n);
            c1 = c1+ 1;
        else
            dis2(c2,1) = x(n);
            dis2(c2,2) = y(n);
            dis2(c2,3) = z(n);
            c2 = c2+1;
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if run ==1
        %plot dislocations separately
        figure(3)
        subplot(1,2,1)
        plot(dis1(:,1),dis1(:,2),'.');
        xlabel('x')
        ylabel('y')
        title('Dislocation1')
        axis equal

        subplot(1,2,2)
        plot(dis2(:,1),dis2(:,2),'.');
        xlabel('x')
        ylabel('y')
        title('Dislocation2')
        axis equal
    else end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if run == 1
        lx = max(x(:,1)) - min(x(:,1));       
        nx = floor(lx/binsize);
        ix = min(x)+binsize/2:binsize: max(x)-binsize/2;  % this line forms ix based on the bin size and x value
        startx = min(x);
    else end
    
    v1 = startx;
    % in the following loop the makelist1 and makelist2 are generating the same output values!! 
    for n = 1:nx
        count = 1;
        v2 = v1 + binsize;
        for nn = 1:length(dis1(:,1))
            if dis1(nn,1) >= v1 && dis1(nn,1) < v2;
                makelist1(count,n) = dis1(nn,3);    %counting the number of atoms in each bin for dis1
                count = count+1;
            else end
        end
        v1 = v2;
    end
    
    v1 = startx;
    
    for n = 1:nx
        count = 1;
        v2 = v1 + binsize;
        for nn = 1:length(dis2(:,1))
            if dis2(nn,1) >= v1 && dis2(nn,1) < v2;  
                makelist2(count,n) = dis2(nn,3);     %counting the number of atoms in each bin for dis2
                count = count+1;
            else end
        end
        v1 = v2;
    end
    %After populating makelist1 and makelist2, it replaces any zero entries with nan
    makelist1(~makelist1) = nan;
    makelist2(~makelist2) = nan;
    
    %it sorts each column of makelist1 and makelist2 based on the values in the first column
    for n = 1:length(makelist1(1,:))
        makelist1(:,n) = sortrows(makelist1(:,n),1);
    end
    
    for n = 1:length(makelist2(1,:))
        makelist2(:,n) = sortrows(makelist2(:,n),1);
    end
    %disp(size(makelist1));
    %disp(size(makelist2));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % MAKE DENSITY PLOTS %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % in the following three lines it interpolates between 
    % min and max in y direction (namely z) 
    iy = linspace(thickcut1,thickcut2,interp); 
    %This line generates a linearly 
    % spaced vector iy with interp points between the values thickcut1 and thickcut2. 
    % The function linspace creates equally spaced values over the specified range.

    ny = length(iy);
    %This line calculates the length of the vector iy and assigns it 
    % to the variable ny. It represents the number of elements in the vector iy.


    %NY = floor(ny/55);
    NY = 2;
    %This line calculates the value of NY by dividing the length of the vector iy (ny) by 55 
    % and then rounding down to the nearest integer using the floor function. The result is 
    % stored in the variable NY.

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % in the following loop Py1 and Py2 are generating zero double files
    % For each bin, it counts the number of atoms in the corresponding 
    % makelist1 and makelist2 columns that fall within the range [ix(n)-binsize, ix(n)+binsize]
    for nn = 1:nx
        for n = 1:ny
            Py1(n,nn) = sum(makelist1(:,nn) >= iy(n)-binsize & makelist1(:,nn) <= iy(n)+binsize);
            Py2(n,nn) = sum(makelist2(:,nn) >= iy(n)-binsize & makelist2(:,nn) <= iy(n)+binsize);
        end
    end
    
    dPy1 = zeros(ny,nx);
    dPy2 = zeros(ny,nx);

    % both of the following dPy1 and dPy2 are generating the same values
    % Calculates the difference between the sum of values in a sliding 
    % window of size NY for both Py1 and Py2.
    for nn = 1:nx
        for n = NY+1:ny-NY
            dPy1(n,nn) = sum(Py1(n:n+NY,nn)) - sum(Py1(n-NY:n,nn));
            dPy2(n,nn) = sum(Py2(n:n+NY,nn)) - sum(Py2(n-NY:n,nn));
        end
        
        MIN1 = find(dPy1(:,nn) == min(dPy1(:,nn)));   % This line finds the indices of the minimum values in the nn-th column 
                                                      % of the matrix dPy1 and stores them in the variable MIN1.
        MAX1 = find(dPy1(:,nn) == max(dPy1(:,nn)));
        MIN2 = find(dPy2(:,nn) == min(dPy2(:,nn)));
        MAX2 = find(dPy2(:,nn) == max(dPy2(:,nn)));


         %Extracts corresponding iy values for the identified minimum and maximum values.
         %Updates the matrices i11, i12, i21, and i22 with the interpolated values.
        i11(run,nn) = iy(min(MIN1));
        i12(run,nn) = iy(max(MAX1));
        i21(run,nn) = iy(min(MIN2));
        i22(run,nn) = iy(max(MAX2));
        %i1ave(run,nn) = (i11(run,nn) + i12(run,nn))/2; % Copy i11 values
        %i2ave(run,nn) = (i21(run,nn) + i22(run,nn))/2; % Copy i11 values
    end    
    % Calculate the size of the matrices
    %[num_runs, num_nn] = size(i11);

    % Initialize i1ave matrix
    %i1ave = zeros(num_runs, num_nn);
    %i2ave = zeros(num_runs, num_nn);
    for nn = 1:nx
    % Concatenate columns of i11 and i12 into i111
        i1ave(:,nn) = (i11(:,nn) + i12(:,nn))/2; % Copy i11 values
        i2ave(:,nn) = (i21(:,nn) + i22(:,nn))/2; % Copy i11 values
    end   

    subplot(2,2,1)
    plot(Py1(:,1),iy)    
    title('Py1,iy')

    subplot(2,2,2)
    plot(Py2(:,1),iy)
    title('Py2,iy')

    subplot(2,2,3)
    plot(dPy1(:,1),iy)
    title('dPy1,iy')

    subplot(2,2,4)
    plot(dPy2(:,1),iy)
    title('dPy2,iy')

    [iX,iY] = meshgrid(ix,iy);
    
    if run == 1
    figure(4)
    
    subplot(2,1,1)
    hold on
    plot(ix,i11(run,:),'k-','linewidth',2)  % black line plot
    plot(ix,i12(run,:),'k-','linewidth',2)  %fluctuation plot
    contour(iX,iY,Py1)
    title('i11,i12')

    axis equal
    hold off
    
    subplot(2,1,2)
    hold on
    plot(ix,i21(run,:),'k-','linewidth',2)  %black line plot
    plot(ix,i22(run,:),'k-','linewidth',2)  %fluctuation plot
    contour(iX,iY,Py2)
    title('i21,i22')

    axis equal
    hold off
    else end
  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if run == 1
        
        figure(5)
        subplot(2,1,1)
        hold on
        f1=plot(dis1(:,1),dis1(:,3),'k.','DisplayName','full dislocation #1');
        f2=plot(ix,i11(run,:),'g-',ix,i12(run,:),'g-','DisplayName','two partials');
        legend('show', 'Location', 'Northeast');

        axis equal
        hold off
        title('Dislocation #1')
        xlim([thickcut1-50 ,thickcut2+50])
        ylim([heightcut1, heightcut2])
        
        subplot(2,1,2)
        hold on
        f3=plot(dis2(:,1),dis2(:,3),'k.','DisplayName','full dislocation #2');
        f4=plot(ix,i21(run,:),'r-',ix,i22(run,:),'r-','DisplayName','two partials');
        legend('show', 'Location', 'Northeast');

        axis equal
        title('Dislocation #2')
        xlim([thickcut1-50 ,thickcut2+50])
        ylim([heightcut1 ,heightcut2])
        hold off
        
        print(plotname,'-dtiff','-r300')
        
    else 
        
    end
    
end

save(savename,'ix','i11','i12','i21','i22','i1ave','i2ave')
