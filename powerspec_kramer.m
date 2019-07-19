function powerspec_kramer(time, csc)
%plot power spectrum

%calulate time info from time
T = time(end) - time(1); %duration
dt = time(2) - time(1); %samp rate

%cut off first sample, which we used to calculate T
time = time(2:end);
csc = csc(2:end);

xf = fft(csc); %fourier transform x

Sxx = 2*dt^2/T * xf.*conj(xf); %compute power spectrum

Sxx = Sxx(1:floor(length(csc)/2)+1); %ignore negatives
Sxx = real(Sxx); %ignore imaginary component

df = 1/max(T); %frequency resolution

fNQ=1/dt/2; %Nyquist frequency

faxis = (0:df:fNQ);

plot(faxis, Sxx)

xlim([0 100]) %alternatively: [0 fNQ]
xlabel('Frequency [Hz]')
ylabel('Power')

end