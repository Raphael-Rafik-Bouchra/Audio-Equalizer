close all
clear 
clc
%-----------------------------------Inputs----------------------------------------------
fileName = uigetfile('*.wav');                  %opens gui to choose the audio file
[inputSignal,InputFs] = audioread(fileName);         %reads the audio file and store at inputSignal

G1 = [];G2 = [];G3 = [];G4 = [];G5 = [];G6 = [];G7 = [];G8 = [];G9 = []; %Initialization
while isempty(G1)
    G1 = input('Enter the Gain of the Frequency band 0-170HZ in DB:\n');  %Every while is to validate the input
end
while isempty(G2)
    G2 = input('Enter the Gain of the Frequency band 170-310HZ in DB:\n');
end
while isempty(G3)
    G3 = input('Enter the Gain of the Frequency band 310-600HZ in DB:\n');
end
while isempty(G4)
    G4 = input('Enter the Gain of the Frequency band 600-1000HZ in DB:\n');
end
while isempty(G5)
    G5 = input('Enter the Gain of the Frequency band 1-3kHZ in DB:\n');
end
while isempty(G6)
    G6 = input('Enter the Gain of the Frequency band 3-6kHZ in DB:\n');
end
while isempty(G7)
    G7 = input('Enter the Gain of the Frequency band 6-12kHZ in DB:\n');
end
while isempty(G8)
    G8 = input('Enter the Gain of the Frequency band 12-14kHZ in DB:\n');
end
while isempty(G9)
    G9 = input('Enter the Gain of the Frequency band 14-16kHZ in DB:\n');
end

filterType = -1;                                                %Initialization for validation
while (filterType ~= 1) && (filterType ~= 2)
    filterType = input('Enter the type of filters used:\n1)FIR    2)IIR\n');
end

OutputFs = -1;                                                  %Initialization for validation
x = ['Input Fs = ' ,num2str(InputFs)];                               %Prints the input samle rate 
disp(x);
while OutputFs <= 0                                         %Input validation    fs >= 2fm    fm = 16k
    OutputFs = input('Enter The output sample rate:\n');
end
%--------------------------------------------------------------------------------------

%----------------------------------------MAIN------------------------------------------
if OutputFs > 32000                                         %Checks if the output sampling rate < 32K
    Fs = OutputFs;                                          %it takes the Fs = 34K and then resamples
else                                                        %the output signal as the desired Output FS.
    Fs = 34000;
end

if filterType == 1                                               %Type 1 means FIR
    N = 299;                                                     % High order to prevent gain manipulation
    [H1,f1,num1] = constructLowPassFIR(Fs,N,170);          %This function constructs a FIR filter
    [H2,f2,num2] = constructBandPassFIR(Fs,N,170,310);
    [H3,f3,num3] = constructBandPassFIR(Fs,N,310,600);
    [H4,f4,num4] = constructBandPassFIR(Fs,N,600,1000);
    [H5,f5,num5] = constructBandPassFIR(Fs,N,1000,3000);
    [H6,f6,num6] = constructBandPassFIR(Fs,N,3000,6000);
    [H7,f7,num7] = constructBandPassFIR(Fs,N,6000,12000);
    [H8,f8,num8] = constructBandPassFIR(Fs,N,12000,14000);
    [H9,f9,num9] = constructBandPassFIR(Fs,N,14000,16000);
    den1=1;den2=1;den3=1;den4=1;den5=1;den6=1;den7=1;den8=1;den9=1;     %den = 1 because it is a FIR filter
else                                                                    %else filter type is IIR
    N = 3;                                                              %with IIR filters use Low order
    [H1,f1,num1,den1] = constructLowPassIIR(Fs,N,170);            %This function constructs IIR filter
    [H2,f2,num2,den2] = constructBandPassIIR(Fs,N,170,310);
    [H3,f3,num3,den3] = constructBandPassIIR(Fs,N,310,600);
    [H4,f4,num4,den4] = constructBandPassIIR(Fs,N,600,1000);
    [H5,f5,num5,den5] = constructBandPassIIR(Fs,N,1000,3000);
    [H6,f6,num6,den6] = constructBandPassIIR(Fs,N,3000,6000);
    [H7,f7,num7,den7] = constructBandPassIIR(Fs,N,6000,12000);
    [H8,f8,num8,den8] = constructBandPassIIR(Fs,N,12000,14000);
    [H9,f9,num9,den9] = constructBandPassIIR(Fs,N,14000,16000);
end

plotFilter(Fs,H1,f1,num1,den1,filterType,1);                     %This function Analyzes each filter
plotFilter(Fs,H2,f2,num2,den2,filterType,2);
plotFilter(Fs,H3,f3,num3,den3,filterType,3);
plotFilter(Fs,H4,f4,num4,den4,filterType,4);
plotFilter(Fs,H5,f5,num5,den5,filterType,5);
plotFilter(Fs,H6,f6,num6,den6,filterType,6);
plotFilter(Fs,H7,f7,num7,den7,filterType,7);
plotFilter(Fs,H8,f8,num8,den8,filterType,8);
plotFilter(Fs,H9,f9,num9,den9,filterType,9);


output1 = filter(num1,den1,inputSignal);                               %This function gets the output when each filter...
output2 = filter(num2,den2,inputSignal);                               %is applied to the input signal
output3 = filter(num3,den3,inputSignal);
output4 = filter(num4,den4,inputSignal);
output5 = filter(num5,den5,inputSignal);
output6 = filter(num6,den6,inputSignal);
output7 = filter(num7,den7,inputSignal);
output8 = filter(num8,den8,inputSignal);
output9 = filter(num9,den9,inputSignal);

plotOutputSignal(Fs,output1,1);                                  %This function plots the signal in time & Freq domain
plotOutputSignal(Fs,output2,2);
plotOutputSignal(Fs,output3,3);
plotOutputSignal(Fs,output4,4);
plotOutputSignal(Fs,output5,5);
plotOutputSignal(Fs,output6,6);
plotOutputSignal(Fs,output7,7);
plotOutputSignal(Fs,output8,8);
plotOutputSignal(Fs,output9,9);


compositeSignal = output1*convDb(G1) + output2*convDb(G2) + output3*convDb(G3) + output4*convDb(G4) ...  %To create the output signal after applying all filters
                  + output5*convDb(G5) + output6*convDb(G6) + output7*convDb(G7) + output8*convDb(G8) + output9*convDb(G9);
plotOutputSignal(InputFs,inputSignal,10);                                            %Plots signal before applying it to filters
plotOutputSignal(OutputFs,compositeSignal,11);                                  %Plots signal after applying it to filters

sound(compositeSignal,OutputFs);                                                %Plays the edited signal
audiowrite('Edited Sound.wav',compositeSignal,OutputFs);                        %saves the output signal
%--------------------------------------------------------------------------------------
    

%-----------------------------------FUNCTIONS------------------------------------------
function [H,f,num,den] = constructLowPassIIR(Fs,N,fc)                           %Constructs low pass IIR filter (BUTTERWORTH)             
    Wc=(fc)/(Fs/2);                                                             %Fs-->sampling, rate N-->order, fc-->CuttOff freq 
    [num,den]=butter(N,Wc);                           
    [H,w]=freqz(num,den);
    f=(w/2/pi)*Fs;
end
function [H,f,num,den] = constructBandPassIIR(Fs,N,fc1,fc2)                     %Constructs Band pass IIR filter (BUTTERWORTH)            
    Wc=[2*fc1/Fs,2*fc2/Fs];                                                     %fc1-->CuttOff freq1, fc2-->CuttOff freq2
    [num,den]=butter(N,Wc);
    [H,w]=freqz(num,den);
    f=(w/2/pi)*Fs;
end

function [H,f,b] = constructLowPassFIR(Fs,N,wc)                                 %Constructs Low pass FIR filter
    Fn=Fs/2;                                                                    %Fs-->sampling, rate N-->order, wc-->CuttOff freq
    fc=wc/Fn;
    b=fir1(N,fc,hanning(N+1));
    [H,w]=freqz(b,1);
    f=(w/2/pi)*Fs;
end

function [H,f,b] = constructBandPassFIR(Fs,N,wl,wu)                             %Constructs Band pass FIR filter
    Fn = Fs/2;                                                                  %wl-->CuttOff freq1, wu-->CuttOff freq2
    wl1 = wl/Fn;
    wu1 = wu/Fn;
    wc = [wl1,wu1];
    b = fir1(N,wc,'bandpass');
    [H,w] = freqz(b,1);
    f=(w/2/pi)*Fs;
end

function plotFilter(Fs,H,f,num,den,t,n)                                         %This function analyzes the filter
    figure;                                                                     
    
    subplot(3,2,1);
    plot(f,abs(H));grid;                                                        %Magnitude Response
    xlabel('Frequency (Hz)');
    ylabel('Magnitude');
    title('Magnitude Response');
    
    subplot(3,2,2);
    plot(f,angle(H)*180/pi);grid;                                               %Phase Response
    xlabel('Frequency (Hz)');
    ylabel('Phase (deg)');
    title('Phase Response');
    
    
    subplot(3,2,3);
    Input1 = impseq(0,0,250);
    Output1 = filter(num,den,Input1);
    time1 = linspace(0,length(Output1)/Fs,length(Output1));
    plot(time1,Output1);grid;                                                   %Impulse Response
    xlabel('Time (Sec)');
    ylabel('Output');
    title('Impulse Response');
    
    subplot(3,2,4);
    Input2 = stepseq(0,0,250);
    Output2 = filter(num,den,Input2);
    time2 = linspace(0,length(Output1)/Fs,length(Output2));
    plot(time2,Output2);grid;                                                   %Step Response
    xlabel('Time (Sec)');
    ylabel('Output');
    title('Step Response');
    
    subplot(3,2,[5 6]);
    zplane(roots(num),roots(den));grid;                                         %Z-PLane
    title('z-plane');                                           
    
    if t == 1                                                                   %Used to print the suitable title
        switch n
        case 1
            suptitle('FIR FILTER 0-170 HZ')
        case 2
            suptitle('FIR FILTER 170-310 HZ')
        case 3
            suptitle('FIR FILTER 310-600 HZ')
        case 4
            suptitle('FIR FILTER 600-1000 HZ')
        case 5
            suptitle('FIR FILTER 1000-3000 HZ')
        case 6
            suptitle('FIR FILTER 3000-6000 HZ')
        case 7
            suptitle('FIR FILTER 6000-12000 HZ')
        case 8
            suptitle('FIR FILTER 12000-14000 HZ')
        case 9
            suptitle('FIR FILTER 14000-16000 HZ')
        end
    else
        switch n
        case 1
            suptitle('IIR FILTER 0-170 HZ')
        case 2
            suptitle('IIR FILTER 170-310 HZ')
        case 3
            suptitle('IIR FILTER 310-600 HZ')
        case 4
            suptitle('IIR FILTER 600-1000 HZ')
        case 5
            suptitle('IIR FILTER 1000-3000 HZ')
        case 6
            suptitle('IIR FILTER 3000-6000 HZ')
        case 7
            suptitle('IIR FILTER 6000-12000 HZ')
        case 8
            suptitle('IIR FILTER 12000-14000 HZ')
        case 9
            suptitle('IIR FILTER 14000-16000 HZ')
        end
    end
end

function plotOutputSignal(Fs,output,n)                                              %Plots the signal in Timme & Freq domain
    time = linspace(0,length(output)/Fs,length(output));
    f = (0:length(output)-1)*100/length(output);
    
    figure;
    subplot(2,2,[1 2]);
    plot(time,output);grid;                                                         %Time domain
    title('Signal In Time Domain');
    xlabel('Time (Sec)');
    ylabel('Output');
    
    subplot(2,2,3);
    foutput = fftshift(fft(output));
    plot(f,abs(foutput));
    title('Magnitude');                                                             %Magnitude
    xlabel('Frequency');
    ylabel('Magnitude');
    
    subplot(2,2,4);
    plot(f,angle(foutput)*180/pi);                                                  %Phase
    title('Phase');
    xlabel('Frequency');
    ylabel('Phase');
    
     switch n                                                                       %Used to prints the suitable title
         case 1
             suptitle('Output (Time-Freq)  input-->0-170Hz Filter')
         case 2
             suptitle('Output (Time-Freq)  input-->170-310Hz Filter')
         case 3
            suptitle('Output (Time-Freq)  input-->310-600Hz Filter')
         case 4
             suptitle('Output (Time-Freq)  input-->600-1000Hz Filter')
         case 5
             suptitle('Output (Time-Freq)  input-->1000-3000Hz Filter')
         case 6
             suptitle('Output (Time-Freq)  input-->3000-6000Hz Filter')
         case 7
             suptitle('Output (Time-Freq)  input-->6000-12000Hz Filter')
         case 8
             suptitle('Output (Time-Freq)  input-->12000-14000Hz Filter')
         case 9
             suptitle('Output (Time-Freq)  input-->14000-16000Hz Filter')
         case 10
             suptitle('Input Signal (Time-Freq)')
         case 11
             suptitle('Composite Output Signal (Time-Freq)')
      end
end

function x = convDb(num)                    %used to convert the gain from DB
    x = 10^(num/20);
end
%--------------------------------------Copied Functions----------------------
    %These functions supposed to be built-in on matlab but in this version they are NOT.

function [x,n] = impseq(n0,n1,n2)                   %To create an impulse function
% Generates x(n) = delta(n-n0); n1 <= n,n0 <= n2
% ----------------------------------------------
% [x,n] = impseq(n0,n1,n2)
%
if ((n0 < n1) || (n0 > n2) || (n1 > n2))
	error('arguments must satisfy n1 <= n0 <= n2')
end
n = [n1:n2];
%x = [zeros(1,(n0-n1)), 1, zeros(1,(n2-n0))];
x = [(n-n0) == 0];
end


function [x,n] = stepseq(n0,n1,n2)                 %To create a unit step function
% Generates x(n) = u(n-n0); n1 <= n,n0 <= n2
% ------------------------------------------
% [x,n] = stepseq(n0,n1,n2)
%
if ((n0 < n1) || (n0 > n2) || (n1 > n2))
	error('arguments must satisfy n1 <= n0 <= n2')
end
n = [n1:n2];
%x = [zeros(1,(n0-n1)), ones(1,(n2-n0+1))];
x = [(n-n0) >= 0];
end
%--------------------------------------------------------------------------------------