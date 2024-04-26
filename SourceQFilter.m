clear all;
clc;
close;

#new username test change for github
#graphics_toolkit("gnuplot")
#plot([0 1 2],[0 1 2])
#close;
#graphics_toolkit("qt")
#set( 0, "defaultlinelinewidth",3)
#plot([0 1 2],[2 1 0])

Tp=0.4
Tn=0.16
glottPulse = @(x) (0).*((x<0)) + (3.*(x./Tp).^2 - 2.*(x./Tp).^3).*((0<=x) & (x<Tp)) + (1-((x-Tp)./Tn).^2).*((Tp<=x) & (x<Tp+Tn)) + (0).*((Tp+Tn<=x));
glott = @(x) (glottPulse(mod(x,1)) - 0.3066675142915666);

Fsflat = 44100*2;
Xflat100 = glott(100*linspace(0, 1, Fsflat));      #sets Xflat100 to calculate Yflat100 for normalizing spectrum later
Yflat100 = fft(Xflat100);      #sets initial Yflat100
Yflat100 = abs(Yflat100)./max(abs(Yflat100));
plot(linspace(0,Fsflat,Fsflat),abs(Yflat100))
axis([0 Fsflat/2])
axis([0 200])
YflatMatrix = reshape(Yflat100(50:Fsflat/2-mod(Fsflat/2-1,100)+148),100,length(50:Fsflat/2-mod(Fsflat/2-1,100)+148)/100);
[m, ~] = max(YflatMatrix);
line = plot(log2((1:length(m))),20*log10(m))

slopeMatrix = (20*log10(m(2:length(m))))./log2(2:length(m));                      #slope in Db/octave
slope = sum(slopeMatrix(2:length(slopeMatrix))/(length(slopeMatrix)-1));
-11.82548879493844     #avg slope over all points
-11.06519211255925     #avg slope 2:14

clf
curve = plot(1:length(m),m)
axis([0 10])
hold on
ndegree = 24;
p = polyfit(1:150,m(1:150),ndegree)
plot(1:150,polyval(p,1:150))

#This is where to loop program on UI update! Before this line is loading screen

Fs=44100;           #sampling rate
T=1/Fs;             #period
L=Fs*2;             #length of domain in samples
#L=length(X)
t=(0:L-1)*T;        #time in seconds
f0=120;             #fundamental frequency
N=50;               #number of iterations

#X = audioread('source.wav');
#X = horzcat(X',zeros(1,(L-length(X))));

#X=0;
#for n = 1:N
#  X=X+((10^(-6/20))^log2(n))*cos(f0*2*pi*t.*(n));
#endfor
#X=X.*0.95/max(X);

Tp=0.4
Tn=0.16
glottPulse = @(x) (0).*((x<0)) + (3.*(x./Tp).^2 - 2.*(x./Tp).^3).*((0<=x) & (x<Tp)) + (1-((x-Tp)./Tn).^2).*((Tp<=x) & (x<Tp+Tn)) + (0).*((Tp+Tn<=x));
glott = @(x) (glottPulse(mod(x,1)) - 0.3066675142915666);

X = glott(f0*t);

Xsum=0;
for n = 1:L
  Xsum = Xsum + X(n);      #ensuring the average of X is near 0
endfor
Xsum=Xsum/L;

waveform = plot(t,X,"linewidth",3);
axis([0 2/f0 -1 1])
xlabel("Time (s)")
ylabel("Amplitude")
title("Two-period Waveform Glott")
playX=audioplayer(X,Fs);
play(playX)
#audiowrite('GlottTest.wav',X,44100)

function transferFunct = formants(fh,F1,Q1,Db1)
  transferFunct =  10^(Db1/20) * (1+(fh./F1 - F1./fh).^2 * Q1^2).^-0.5;
endfunction

Y=fft(X);
Yh=Y(1:L/2);               #removes reflection of fft
f=Fs/L*(1:length(Y));      #frequency domain for fft
fh=Fs/L*(1:length(Yh));    #frequencies up to nyquist frequency

plot(fh,abs(Yh))
hold on
plot(fh,max(abs(Yh))*((10.^(-11.06519211255925/20)).^log2(fh./f0)))
hold on
plot(f0*(1:length(m)),max(abs(Yh))*m)
axis([0 2000 0 max(abs(Yh))])
xlabel("Frequecy (Hz)")
ylabel("Amplitude")
title("Spectrogram of Glott and 11.8 dB Rolloff")

#splitlength =
#raise Y peaks to same level

transform = formants(fh,500,10,0);      #linear plot of transform
plot(fh,transform)
axis([0 N*f0+f0])
xlabel("Frequency (Hz)")
ylabel("Amplitude")
title("Linear Plot of Transform")

plot(fh,20*log10(transform))
hold on
#plot(fh,log10(2^-0.5))


