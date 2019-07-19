function cohere_kramer(time, csc1, csc2)
%from kramer 2013

%Define the number of indices per trial.
N = size(csc1,2); 

%Define the sampling interval.
dt = time(2)-time(1);

%Define the duration of data.
T = time(end);

%Create variables to save the spectra.
K = size(csc1,1);
Sxx = zeros(K,N);
Syy = zeros(K,N);
Sxy = zeros(K,N);

%Compute the spectra for each trial.
for k=1:K
    Sxx(k,:) = 2*dt^2/T * fft(csc1(k,:)) .* conj(fft(csc1(k,:)));
    Syy(k,:) = 2*dt^2/T * fft(csc2(k,:)) .* conj(fft(csc2(k,:)));
    Sxy(k,:) = 2*dt^2/T * fft(csc1(k,:)) .* conj(fft(csc2(k,:)));
end

%Ignore negative frequencies.
Sxx = Sxx(:,1:N/2+1);
Syy = Syy(:,1:N/2+1);
Sxy = Sxy(:,1:N/2+1);

%Average the spectra across trials.
Sxx = mean(Sxx,1);
Syy = mean(Syy,1);
Sxy = mean(Sxy,1);

%Compute the coherence.
cohr = abs(Sxy) ./ (sqrt(Sxx) .* sqrt(Syy));

%Determine the frequency resolution.
df = 1/max(T);
%Determine the Nyquist frequency.
fNQ = 1/ dt / 2;


%Construct frequency axis.
faxis = (0:df:fNQ);

%Plot the results
plot(faxis, real(cohr));

xlim([0 100]); ylim([0 1])
xlabel('Frequency [Hz]')
ylabel('Coherence [ ]')