function [Ug,F,T,U,y] = DispersionCurve(x,fs,alpha,fini,fend,df,dist,Ulim,levels)

% A function to calculate and plot the dispersion curve (group velocity as a function of wave period) of a wave packet contained
% in an input signal. It is calculated numerically by applying the FTAN (Frequency-Time Analysis) method proposed by Levshin et al in 
% a 1972 entitled "On a frequency-time analysis of oscillations" for applications such as the analysis of seismic surface waves.

% Author: Alejandro Silva, Universidad Politécnica de Madrid, Spain
% (alejandro.silva@upm.es), April 2025.

% INPUT PARAMETERS
% x: One-dimensional column array containing the signal with the wave.
% fs: Signal sampling frequency (Hz).
% alpha: A parameter that determines the time-frequency resolution of the filtering.
% Higher values for higher frequency resolution, lower values for higher time resolution.
% fini: The lowest filtering central frequency. It determines the upper bound of the curve's period (horizontal) axis.
% fend: The highest filtering central frequency. It determines the lower bound of the curve's period (horizontal) axis.
% df: The spacing between two consecutive filtering central frequencies.
% dist: The distance (usually, in km) between the wave source and the signal's measurement point
% Ulim: The upper limit of the group velocity for plotting (vertical axis)
% levels: The number of contour level curves in the figure showing the dispersion curve and the FTAN output function

% OUTPUTS
% Ug: One-dimensional array that contains the dispersion curve representing the wave group velocities (in units of distance/s).
% F: Vector containing the wave frequency (in Hz) that corresponds to each value of the dispersion curve array Ug.
% T: Vector containing the wave period (in s) that corresponds to each value of the dispersion curve array Ug.
% U: The vector of wave group velocities for plotting (vertical axis).
% y: Two-dimensional array containing the FTAN output function
% (it represents the log-scale instantaneous amplitude of the filtered signal for each filtering central frequency in time domain).

F = fini:df:fend; % The central filtering frequencies (in Hz).
T = 1./F; % Converts central filtering frequency values to their corresponding periods (in s).
f = linspace(0,fs,length(x)); % Vector of frequencies (in Hz) of the signal spectrum.
t = linspace(1/fs,length(x)/fs,length(x)); % Vector of times (in s) of the input signal.
y = zeros(length(x),length(F)); % Initializes the FTAN output function.
X = fft(x); % Fourier spectrum of the input signal.
for k=1:length(F) % For each central filtering frequency.
    H = exp(-alpha*(f-F(k)).^2/F(k)^2)'; % The filter's transfer function (a real-valued Gaussian in the frequency domain).
    Y = H.*X; % Filters the signal in the frequency domain.
    y(:,k) = ifft(Y); % The FTAN output function array is filled column-wise with the Inverse Fourier Transform of filter output Y.
end
% Outputs the index of the highest log-absolute value of y for each central frequency (i.e. for each column in y).
[~,indicesmax] = max(log(abs(y)));
% Wave arrival times are related with their corresponding group velocities through the distance between the source and the sensor.
Ug = dist./t(indicesmax);
U = dist./t; % Idem.

%% PLOTS THE LOG OF THE MODULE OF THE FTAN FUNCTION AND THE ESTIMATED DISPERSION CURVE
figure
contour(T,U,log(abs(y)),levels)
hold on
scatter(T,Ug)
plot(T,Ug)
xlabel('wave period (s)')
ylabel('Group velocity (units of distance/s)')
ylim([U(end) Ulim])
grid on
colorbar

end