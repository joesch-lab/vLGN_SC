function [lineOut, fillOut] = stdshade_mean_std_linecolor(m,s,alpha,acolor,linecolor,F,smth)
% usage: stdshade_mean_std(m.s,alpha,acolor,F,smth)
% m is mean, s is std
% plot mean and sem/std coming from a matrix of data, at which each row is an
% observation. sem/std is shown as shading.
% - acolor defines the used color (default is red) 
% - F assignes the used x axis (default is steps of 1).
% - alpha defines transparency of the shading (default is no shading and black mean line)
% - smth defines the smoothing factor (default is no smooth)
% smusall 2010/4/23

%The code was downloaded from:
%  Simon Musall (2024). stdshade (https://www.mathworks.com/matlabcentral/fileexchange/29534-stdshade), MATLAB Central File Exchange. Abgerufen9. Dezember 2024. 

% modified to provide the mean (to plot as a line) and SD or SEM (to plot as shaded area)

if exist('acolor','var')==0 || isempty(acolor)
    acolor='r'; 
end

if exist('linecolor','var')==0 || isempty(linecolor)
    linecolor='k'; 
end

if exist('F','var')==0 || isempty(F)
    F=1:length(m);
end

if exist('smth','var'); if isempty(smth); smth=1; end
else smth=1; %no smoothing by default
end  

if ne(size(F,1),1)
    F=F';
end

amean = m; %get man over first dimension
if smth > 1
    amean = boxFilter(m,smth); %use boxfilter to smooth data
end
astd = s; % to get std shading
% astd = nanstd(amatrix,[],1)/sqrt(size(amatrix,1)); % to get sem shading

if exist('alpha','var')==0 || isempty(alpha) 
    fillOut = fill([F fliplr(F)],[amean+astd fliplr(amean-astd)],acolor,'linestyle','none','HandleVisibility','off');
    acolor='k';
else
    fillOut = fill([F fliplr(F)],[amean+astd fliplr(amean-astd)],acolor, 'FaceAlpha', alpha,'linestyle','none','HandleVisibility','off');
end

if ishold==0
    check=true; else check=false;
end

hold on;
lineOut = plot(F,amean, 'color', linecolor,'linewidth',1); %% change color or linewidth to adjust mean line

if check
    hold off;
end

end


function dataOut = boxFilter(dataIn, fWidth)
% apply 1-D boxcar filter for smoothing

fWidth = fWidth - 1 + mod(fWidth,2); %make sure filter length is odd
dataStart = cumsum(dataIn(1:fWidth-2),2);
dataStart = dataStart(1:2:end) ./ (1:2:(fWidth-2));
dataEnd = cumsum(dataIn(length(dataIn):-1:length(dataIn)-fWidth+3),2);
dataEnd = dataEnd(end:-2:1) ./ (fWidth-2:-2:1);
dataOut = conv(dataIn,ones(fWidth,1)/fWidth,'full');
dataOut = [dataStart,dataOut(fWidth:end-fWidth+1),dataEnd];

end

