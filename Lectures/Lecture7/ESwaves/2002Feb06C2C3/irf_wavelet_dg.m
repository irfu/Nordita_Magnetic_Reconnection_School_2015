function specrec=irf_wavelet_dg(varargin)
%IRF_WAVELET  compute fast wavelet spectrogram
%  Calculate wavevelet spectrogram based on fast FFT algorithm
%  Written by A. Tjulin
%
% specrec=IRF_WAVELET(data) data are column vectors, first column assumed time
%
% specrec=IRF_WAVELET(data,parameter1,parameter1_value,...)
%   input parameters to the routine can be
%       'Fs' - sampling rate, if not specified calculated from time series, if no time series assumed 1S/s
%       'f'  - vector [fmin fmax], calculate spectra between frequencies fmin and fmax
%       'nf' - number of frequency bins
%       'wavelet_width' - width of Morlet wavelet, default 5.36. (same as 'w')
%
%	SPECREC is a structure:
%		SPECREC.T - time
%		SPECREC.P - spectrum
%		SPECREC.F - frequency
%
% if no output specified, plot the spectrogram
%
% See also IRF_SPECTROGRAM

% With modification by D. B. Graham

%% Check the input
if nargin == 0,
    help irf_wavelet;
    return;
end
data=varargin{1};
args=varargin(2:end);
if numel(args)>0,
    flag_have_options=1;
else
    flag_have_options=0;
end
%% Default values
% Fs
if size(data,1)==1, % no time axis specified and row vector input
    Fs=1;
    t=1:numel(data);
    t=t(:);
else                % assume time axis is first column
    if size(data,2) == 1 % no time axis specified, column vector input
        Fs=1;
        t=1:numel(data);
        t=t(:);
    else
        Fs=1/(data(2,1)-data(1,1));
        t=data(:,1);
        data(:,1)=[];
    end
end
% f
amin=0.01; % The highest frequency to consider is 0.5*sampl/10^amin
fmax=0.5*Fs/10^amin;
amax=2; % The lowest frequency to consider is 0.5*sampl/10^amax
fmin=0.5*Fs/10^amax;
% nf
nf=200;
% wavelet_width
wavelet_width=5.36;
%% Check for NaNs
inan=isnan(data);
data(inan)=0;

%% Check the options
while flag_have_options
    l = 2;
    switch(lower(args{1}))
        case 'Fs'
            if numel(args)>1 && isnumeric(args{2})
                Fs = args{2};
            else
                irf_log('fcal','parameter ''Fs'' without parameter value')
            end
        case 'nf'
            if numel(args)>1 && isnumeric(args{2})
                nf = args{2};
            else
                irf_log('fcal','parameter ''nf'' without parameter value')
            end
        case {'wavelet_width','w'}
            if numel(args)>1 && isnumeric(args{2})
                wavelet_width = args{2};
            else
                irf_log('fcal','parameter ''wavelet_width'' without parameter value')
            end
        case 'f'
            if numel(args)>1 && isnumeric(args{2})
                if numel(args{2})== 2,
                    fmin=max(fmin,args{2}(1)); 
                    fmax=min(fmax,args{2}(2));
                else
                    irf_log('fcal','parameter ''f'' should have vector with 2 elements as parameter value.')
                end
            else
                irf_log('fcal','parameter ''f'' without parameter value')
            end
        otherwise
            irf_log('fcal',['Unknown flag: ' args{1}])
            l=1;
            break
    end
    args = args(l+1:end);
    if isempty(args), flag_have_options=0; end
end

sampl=Fs;
w0=sampl/2; % The maximum frequency
anumber=nf; % The number of frequencies
sigma=wavelet_width/(Fs/2); % The width of the Morlet wavelet
amin=log10(0.5*Fs/fmax); % The highest frequency to consider is 0.5*sampl/10^amin
amax=log10(0.5*Fs/fmin); % The lowest frequency to consider is 0.5*sampl/10^amax
a=logspace(amin,amax,anumber);

%% Remove the last sample if the total number of samples is odd

if size(data,1)/2 ~= floor(size(data,1)/2)
    data(end,:)=[];
    t(end)=[];
end

%% Find the frequencies for an FFT of all data

nd2=size(data,1)/2;
nyq=1/2;
freq=sampl*(1:nd2)/(nd2)*nyq;
w=[0,freq,-freq(end-1:-1:1)];% The frequencies corresponding to FFT


%% construct specrec
specrec.p=cell(1,size(data,2));
specrec.t=t;
% Get the correct frequencies for the wavelet transform
newfreq=w0./a;
specrec.f=newfreq(:);
[newfreqmat,temp]=meshgrid(newfreq,w);
[temp,ww]=meshgrid(a,w); % Matrix form
for i=1:size(data,2), % go through all the datacolumns
    %% Make the FFT of all data
    datacol=data(:,i);
    Sw=fft(datacol);
    %  Sw=fft(detrend(Ex));
    [aa,Sww]=meshgrid(a,Sw); % Matrix form
    
    %% Calculate the FFT of the wavelet transform
    
    Ww=sqrt(1).*Sww.*exp(-sigma*sigma*((aa.*ww-w0).^2)/2);
    
    %% Get the wavelet transform by IFFT of the FFT
    
    W=ifft(Ww);
    
    %% Calculate the power spectrum
    
    %power=abs((2*pi)*conj(W).*W./newfreqmat);
    
    %% Remove data possibly influenced by edge effects
    
    %censur=floor(2*a);
    %power2=power;
    %for j=1:anumber;
    %    power2(1:censur(j),j)=NaN;
	%	power2(numel(datacol)-censur(j):numel(datacol),j)=NaN;
    %end
    specrec.p{i}=W;
	
	%% remove NaNs
	
	specrec.p{i}(inan(:,i),:)=NaN;
end

if nargout==0,
    irf_spectrogram(specrec)
end
