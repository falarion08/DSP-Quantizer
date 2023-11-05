%Read audio file

[audio, audioFs] = audioread('Secret-Message.mp3');
% sound(audio,audioFs *1); 

% Stores quantization step-size or bit-depth
% n = input("Input bit-depth: ");
n = 5; 
time =linspace(0,3,length(audio));

% Part 1
% xq = midTreadQuintizer(n,audio);
% e = quantizationError(xq,audio);
% 
% 
% subplot(3,1,1)
% plot(time,audio);
% title('Input Audio');
% xlabel('Time in seconds')
% 
% subplot(3,1,2)
% plot(time,xq);
% title('Quantized Audio(Bit Depth = 5)');
% xlabel('Time in seconds')
% 
% subplot(3,1,3)
% plot(time,e);
% title('Quantization Error');
% xlabel('Time in seconds');
% 
% SNR1 = signalToNoiseRatio_RMS(audio,n);
% SNR2 = signalToNoiseRatio(e,audio);

% Part 4 - mu Law Compander Implementation
xCompressed = muCompressor(audio);
xq = midTreadQuintizer(n,xCompressed);
recovered = muExpander(xq,audio);
e = quantizationError(recovered, audio);

subplot(4,1,1);
plot(time,audio);
title('Input Audio');
xlabel('time in seconds');

subplot(4,1,2)
plot(time,xq);
title('Quantized Audio (Bit-Depth = 5)');
xlabel('time in seconds');

subplot(4,1,3)
plot(time,recovered);
title('Recovered');
xlabel('time in seconds');

subplot(4,1,4)
plot(time,e);
title('Error');
xlabel('time in seconds');

SNR2 = signalToNoiseRatio(e,audio); 



function xq = midTreadQuintizer(n,audio)
    xMin = min(audio(:,1));
    xMax = max(audio(:,1)); 
    L = (2^n - 1);

    % Get quantization interval
    qStep = round(xMax - xMin,4) /L;

    audioLen = length(audio); 
    
    % Initialize output
    xq = zeros(audioLen,2);

    % Compute for recovered amplitude
    for i=1:audioLen
        xq(i,:) = round(audio(i,1) / qStep) * qStep; 
    end
    
end

function e = quantizationError(quantizedAudio, inputAudio)
audioLen = length(inputAudio); 
e = zeros(audioLen,2);
    for i=1:audioLen
        e(i,:) = quantizedAudio(i) - inputAudio(i); 
    end

end

function SNR = signalToNoiseRatio_RMS(audio,n)
% SNR computation using rms value of signal, x
audioLen = length(audio); 
sum = 0;

xMin = min(audio(:,1));
xMax = max(audio(:,1)); 
L = (2^n - 1);

% Get quantization interval
qStep = round(xMax - xMin,4)/L;

% Computing for RMS value
for i=1:audioLen
    sum = sum + (audio(i,1) * audio(i,1));
end

rms = sqrt((1/audioLen)* sum);
SNR = 10.79 + 20*log10(rms/qStep);

end

function SNR = signalToNoiseRatio(qError, input)
    numerator = 0;
    denominator = 0;
    
    audioLen = length(input);
    for i =1:audioLen
        numerator = numerator + input(i,1)^2;
        denominator = denominator + qError(i,1)*qError(i,1);
    end

    SNR = 10 * log10(numerator / denominator);
end

function muComp = muCompressor(audio)
    xMax = max(audio(:,1)); 
    audioLen = length(audio); 
    muComp = zeros(audioLen,2); 
    mu = 255; 

    for i=1:audioLen
        x = audio(i,1);
        if(x >= 0)
            sign = 1;
        else 
            sign = -1;
        end
        muComp(i,:) = sign * log(1 + mu*(abs(x)/abs(xMax))) / log(i + mu);
    end
end

function muExpan = muExpander(quantizedAudio, inputAudio)
    xMax = max(inputAudio(:,1));
    mu = 255; 
    
    audioLen = length(inputAudio);
    muExpan = zeros(audioLen,2);

    for i = 1:audioLen
        y = quantizedAudio(i,1); 

        if y>= 0
            sign =1;
        else
            sign = -1;
        end

        muExpan(i,:) = abs(xMax) * sign * ((1+mu)^abs(y) - 1)/ mu; 
    end 
end


% Part 5 a
function encoded = myBinaryEncoder(sample)
    sample = min(max(sample, -1), 1);
    if sample >= 0
        sign_bit = '1';
    else
        sign_bit = '0';
        sample = -sample; 
    end
    mantissa = '';
    for i = 1:11 
        sample = sample * 2;
        if sample >= 1
            mantissa = strcat(mantissa, '1');
            sample = sample - 1;
        else
            mantissa = strcat(mantissa, '0');
        end
    end
        encoded = strcat(sign_bit, mantissa);
end


% part 5 b

function compressed8bit = myMuCompressor(linear12bit)
    % condition for the input to be a 12-bit string
    if ~ischar(linear12bit) || length(linear12bit) ~= 12
        error('Input must be a 12-bit string.');
    end
    
    % Extract the sign bit
    signBit = linear12bit(1);
    segmentBits = ['0', '0', '0'];
    remainingBits = ['0', '0', '0', '0'];
    
    if linear12bit(2) == '0' && linear12bit(3) == '0' && linear12bit(4) == '0' && linear12bit(5) == '0' && linear12bit(6) == '0' && linear12bit(7) == '0' && linear12bit(8) == '0'
        segmentBits = ['0', '0', '0'];
        remainingBits = [linear12bit(9), linear12bit(10), linear12bit(11), linear12bit(12)];
    end 

    if linear12bit(2) == '0' && linear12bit(3) == '0' && linear12bit(4) == '0' && linear12bit(5) == '0' && linear12bit(6) == '0' && linear12bit(7) == '0' && linear12bit(8) == '1'
        segmentBits = ['0', '0', '1'];
        remainingBits = [linear12bit(9), linear12bit(10), linear12bit(11), linear12bit(12)];
    end 

    if linear12bit(2) == '0' && linear12bit(3) == '0' && linear12bit(4) == '0' && linear12bit(5) == '0' && linear12bit(6) == '0' && linear12bit(7) == '1'
        segmentBits = ['0', '1', '0'];
        remainingBits = [linear12bit(8), linear12bit(9), linear12bit(10), linear12bit(11)];
    end 
    if linear12bit(2) == '0' && linear12bit(3) == '0' && linear12bit(4) == '0' && linear12bit(5) == '0' && linear12bit(6) == '1' 
        segmentBits = ['0', '1', '1'];
        remainingBits = [linear12bit(7), linear12bit(8), linear12bit(9), linear12bit(10)];
    end 
    if linear12bit(2) == '0' && linear12bit(3) == '0' && linear12bit(4) == '0' && linear12bit(5) == '1' 
        segmentBits = ['1', '0', '0'];
        remainingBits = [linear12bit(6), linear12bit(7), linear12bit(8), linear12bit(9)];
    end 
    if linear12bit(2) == '0' && linear12bit(3) == '0' && linear12bit(4) == '1'
        segmentBits = ['1', '0', '1'];
        remainingBits = [linear12bit(5), linear12bit(6), linear12bit(7), linear12bit(8)];
    end 
    if linear12bit(2) == '0' && linear12bit(3) == '1' 
        segmentBits = ['1', '1', '0'];
        remainingBits = [linear12bit(4), linear12bit(5), linear12bit(6), linear12bit(7)];
    end 
    if linear12bit(2) == '1'
        segmentBits = ['1', '1', '1'];
        remainingBits = [linear12bit(3), linear12bit(4), linear12bit(5), linear12bit(6)];
    end 


    % Combine all parts to get the final 8-bit mu-law compressed code
    compressed8bit = [signBit, segmentBits, remainingBits];
end

%sample inputs
% compressed8bit = myMuCompressor('100000000101') output = 10000101
% compressed8bit = myMuCompressor('000011101010') output = 01001101

