//CTFT EXP 7 rect pulse 
clc
clear all
wid=5;
t=-5:0.01:5;
x=zeros(size(t));
x(t>=-wid/2&t<=wid/2)=1;
subplot(1,2,1);
plot(t,x,'r', 'linewidth',2);

axis([-20 20 0 10]);
xlabel('Time');
ylabel('Amplitude');
title('Continus Rectangular Pulse');
legend ('Eihab 102395005');
wid=4;
t=-5:0.01:5;
s=sinc(t);
subplot(1,2,2);
plot(t,s, 'b', 'linewidth',2);
axis([-10 10 -3 5]);
xlabel ('Time');
ylabel('Amplitude');
title('Sinc Function');
legend ('Eihab 102395005');
a=trapz(t,x);


//dtft 8th
t=-10:0.1:10;
y=rectpuls(t,1);
z=abs(y);
a=fft(y);
plot(t,fftshift(abs(a)));



//9th z trans
clc
clear all;
b=[1];
a=[1,-0.9];
subplot(2,2,1);
zplane(b,a);
w=-6* pi:6*pi;
x=freqz(b,a,w);
T=abs(x);
angle=angle(x);
subplot(2,2,2);
plot(w/pi,T);
xlabel('Normalised Freq');
ylabel('Magnitude (db)');
subplot(2,2,3);
plot(w/pi,angle);
xlabel('Normalised Freq');
ylabel ('phase');




//10th polezero
clc
clear all;
syms n z
x(n) = ((1/2).^n)+(2.^n);
z=ztrans(x,n,z);
subplot(2,2,1);
b=[4,-5];
a=[2,-5,2];
zplane(b,a)





//11th
% DFT and frequency spectrum

clear all,

% Define the signal length
N = 4;

% Imaginary unit
i = sqrt(-1);

% Define the signal sequence
xn = [1, 1, 2, 3];

% Initialize the DFT sequence to zeros
XK = zeros(1, N);

% Loop to compute DFT
for k = 0:1:N-1
  for n = 0:1:N-1
    XK(k+1) = XK(k+1) + xn(n+1) * exp(-j*2*pi*k*n/N);
  end
end

% Display the DFT sequence
disp('The DFT sequence is: ');
disp(XK);

% Compute the magnitude of the DFT
MagXK = abs(XK);

% Define the frequency range
k = 0:1:N-1;

% Plot the magnitude spectrum
subplot(2,1,1);
stem(k, MagXK);
title('Magnitude Spectrum');
xlabel('K');
ylabel('Magnitude');

% Plot the phase spectrum
subplot(2,1,2);
stem(k, angle(XK));
title('Phase Spectrum');
xlabel('K');
ylabel('Phase');

% DFT computation using FFT

clear all,

% Define the signal length
N = 4;

% Define the signal sequence
xn = [1, 1, 2, 3];

% Display message
disp('fft of the sequence');

% Compute DFT using FFT
Xk = fft(xn, N);

% Compute the magnitude of the DFT
MagXk = abs(Xk);

% Compute the phase of the DFT
Phaxk = angle(Xk);

% Display messages
disp('The magnitude sequence is: ');
disp(MagXk);
disp('The phase sequence is: ');
disp(Phaxk);

% Define the frequency range
n = 0:1:N-1;
wk = 2*pi*n/N;

% Plot the magnitude spectrum
subplot(2,2,2);
stem(wk, MagXk);
title('Magnitude Spectrum');
xlabel('Frequency');
ylabel('Magnitude');

% Plot the phase spectrum
subplot(2,2,2);
stem(wk, Phaxk);
title('Phase Spectrum');
xlabel('Frequency');
ylabel('Phase');

