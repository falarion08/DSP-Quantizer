%Read audio file

[audio, audioFs] = audioread('Secret-Message.mp3');
% sound(audio,audioFs *1); 

% Stores quantization step-size or bit-depth
% n = input("Input bit-depth: ");
n = 3; 

yq = length(audio);


function yq = midTreadQuintizer(n,audio)
    xMin = min(audio);
    xMax = max(audio); 

    qStep = (xMax(1) - xMin(1)) / (2^n - 1);
    
    audioLen = length(audio);
    yq = zeros(audioLen,2);

    for i=1:audioLen
        xq = audio(i); 
        j = 0; 
        while (j <= n)
            % % upperThresold = j + ();  
            % lowerThreshold
            % if ( <= j )
            % end
            j = j + 1; 
        end
    end
end

