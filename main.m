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

%Part 5
mBE12bit = myBinaryEncoder(audio);
dCompressed = dCompressor(mBE12bit);
dExpanded = dExpander(dcompressed);
dDecoded = myBinaryDecoder(dExpanded);
audio=audio(:,1);
dError = dDecoded-audio;
dDelta = (input_max-input_min)/((2^n)-1);
SNR_recovered = signalToNoiseRatio(dError, audio);

figure;
subplot(3,1,1)
plot(time,audio)
title('Input Audio');
xlabel('time in seconds');

subplot(3,1,2)
plot(time,audio)
title('Recovered Signal');
xlabel('time in seconds');

subplot(3,1,3)
plot(time,dError)
title('Error');
xlabel('time in seconds');

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

function SNR = signalToNoiseRatio_RMS(audio, qStep)
    % Calculate the SNR using the first formula
    rms_value = rms(audio);
    SNR = 10.79 + 20 * log10(rms_value / qStep);
end

function SNR = signalToNoiseRatio(qError, audio)
    % Calculate the SNR using the second formula
    pSignal = sum(audio.^2);
    pError = sum(qError.^2);
    SNR = 10 * log10(pSignal / pError);
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

function myBin = myBinaryEncoder(input)
    myBin = zeros(size(input, 1), 12);
    for i = 1:size(input, 1)
        curr = input(i);
        temp = zeros(1, 12);
        if curr < 0
            temp(1) = 0;
            curr = -curr; % Make it positive for encoding
        else
            temp(1) = 1;
        end

        for pos = 2:12
            curr = curr * 2;
            temp(pos) = floor(curr);
            curr = curr - temp(pos);
        end

        myBin(i, :) = temp;
    end
end

function deComp = dCompressor(input)
    deComp = zeros(size(input, 1), 8);
    for i = 1:size(input, 1)
        temp = zeros(1, 8);
        curr = input(i, :);
        temp(1) = curr(1);
        temp(2:8) = [0 0 0 curr(2:5)];

        for pos = 2:8
            if curr(pos) == 1
                offset = 9 - pos;
                temp(5:8) = curr(pos + 1:pos + 4);
                binaryStr = dec2bin(offset, 3);
                for di = 1:3
                    temp(1 + di) = str2double(binaryStr(di));
                end
                break
            end
        end
        deComp(i, :) = temp;
    end
end

function dExpander = dExpander(input)
    dExpander = zeros(size(input, 1), 12);
    for i = 1:size(input, 1)
        temp = zeros(1, 12);
        curr = input(i, :);
        temp(1) = curr(1);
        segment = 7 - bin2dec([num2str(curr(2)) num2str(curr(3)) num2str(curr(4))]);

        if segment ~= 7
            temp(segment + 2) = 1;
            temp(segment + 3:segment + 6) = curr(5:8);

            if segment == 5
                temp(segment + 7) = 1;
            elseif segment < 5
                temp(segment + 7) = 1;
                temp(segment + 8:12) = zeros(1, 12 - segment - 7);
            end
        else
            temp(segment + 2:segment + 5) = curr(5:8);
        end

        dExpander(i, :) = temp;
    end
end

function deBin = myBinaryDecoder(input)
    deBin = zeros(size(input, 1), 1);
    for i = 1:size(input, 1)
        curr = input(i, :);
        temp = 0;

        for exp = 2:12
            if curr(exp) == 1
                temp = temp + 2^-(exp - 1);
            end
        end

        if curr(1) == 0
            temp = -temp;
        end

        deBin(i, 1) = temp;
    end
end
