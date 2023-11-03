%Read audio file

[audio, audioFs] = audioread('Secret-Message.mp3');
% sound(audio,audioFs *1); 

% Stores quantization step-size or bit-depth
% n = input("Input bit-depth: ");
n = 5; 

xq = midTreadQuintizer(n,audio);

subplot(2,1,1)
plot(audio);

subplot(2,1,2)
plot(xq);

function xq = midTreadQuintizer(n,audio)
    xMin = min(audio(:,1));
    xMax = max(audio(:,1)); 
    L = (2^n - 1);
    xq = 0;


    qStep = round(xMax - xMin) /L;  
    audioLen = length(audio); 
    
    xq = zeros(audioLen,2);

    for i=1:audioLen
        xq(i,:) = round(audio(i,1) / qStep) * qStep; 
    end
    
end

