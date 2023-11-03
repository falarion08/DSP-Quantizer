%Read audio file

[audio, audioFs] = audioread('Secret-Message.mp3');
% sound(audio,audioFs *1); 

% Stores quantization step-size or bit-depth
% n = input("Input bit-depth: ");
n = 5; 

% Part 1
xq = midTreadQuintizer(n,audio);
e = quantizationError(xq,audio);

subplot(3,1,1)
plot(audio);
title('Input Audio');

subplot(3,1,2)
plot(xq);
title('Quantized Audio(Bit Depth = 5)');

subplot(3,1,3)
plot(e);
title('Quantization Error');


function xq = midTreadQuintizer(n,audio)
    xMin = min(audio(:,1));
    xMax = max(audio(:,1)); 
    L = (2^n - 1);
    xq = 0;

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

