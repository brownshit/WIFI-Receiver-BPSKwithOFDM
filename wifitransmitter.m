%% wifitransmitter: Creates a Wi-Fi packet
% output = WiFi packet
% Inputs: message = text message, snr = signal to noise ratio,
% level = number of stages of encoding
function output = wifitransmitter(message, level, snr)
%% Default values
if(nargin < 2)
    level = 5;
end
if(nargin < 3)
    snr = Inf;  % m. Noise Power == 0
end

% m. currently snr >= -5 works
%% Sanity checks

% check if message length is reasonable
if(length(message) > 10000)
    fprintf(2, 'Error: Message too long\n');
    output=[];
    return;
end

% check if level is between 1 and 5
if(level > 5 || level < 1)
    fprintf(2, 'Error: Invalid level, must be 1-5\n');
    output=[];
    return;
end

%% Some constants

% We will split the data into a cluster of nfft bits
nfft = 64;
% m. nfft indicates the number of samples in FFT, when you do OFDM

% This is the Encoder/decoder trellis used by WiFi's turbo encoder
Trellis = poly2trellis(7,[133 171]);
% m. constraint length K = 7
% m. 133 (octal representaion) (=> 1011011) and 171 (=> 1111001) shows
% m. connections from the
% m. outputs of the registers to the two adders
% m. k = 1, n = 2 (2 output bits for every single bit)
% m. Explain about trellis coding

% Every WiFi packet will start with this exact preamble (m. 64 bit long)
preamble = [1, 1, 1, 1,-1,-1, 1, 1,-1, 1,-1, 1, 1, 1, 1, 1, 1,-1,-1, 1, 1,-1, 1,-1, 1, 1, 1, 1, 1,-1,-1, 1, 1,-1, 1,-1, 1,-1,-1,-1,-1,-1, 1, 1,-1, -1, 1,-1, 1,-1, 1, 1, 1, 1,-1,-1, -1,-1,-1, 1, 1,-1, -1, 1];
% m. Explain about preamble

% Every 64 bits are mixed up like below:
Interleave = reshape(reshape([1:nfft], 4, []).', [], 1);
% m. .': transpose
% m. used in level 2

%% Lets learn about the message
% Length
len = length(message);
% m. a character is counted as a length 1

%% Level #1: First lets do coding, which adds redundancy to the bits
if (level >= 1)
    % This basically converts the message into a sequence of bits
    bits = reshape(dec2bin(double(message).', 8).', 1, [])-'0';
    %bits = reshape(dec2bin(double(message).', 8).', 1, [])
    % m. double('he') returns ascii values of chars; 104 101.
    % m. dec2bin('he',8) converts chars into 2 x 8 char array;
    % '01101000' (corresponding to 104)
    % '01100101' (corresponding to 101)
    % , where each bit is a single character.
    % m. reshape(dec2bin(double('he').', 8).', 1, []) returns 1 x 16 array;
    % '0110100001100101'
    % , where each bit is a single character.
    % m. subtracting '0' (= 48 in ascii) from a char array made of '0' and
    % '1' results in an integer array made of 0s and 1s.
    % m. a single char '0' in subtrahend (B in A - B) subtracts from every
    % character and outputs the result in integer array.
    % m. if the subtrahend is made of multiple chars, the number of chars
    % in subtrahend must match the chars of minuend (A in A - B).
    % m. length: 8 * length(message).
    % m. for example, two letters produce 16 integers (made of 0s and 1s).
    
    % We append as many bits as necessary to make this a multiple of
    % nfft
    bits = [bits, zeros(1, mod(-length(bits), nfft))];
    % m. bit length: ceil(8 * length(message) / nfft) * nfft
    % m. # of bits is approximately 8 X no. of alphabets in message

    % Next, we apply the turbo coder
    output = convenc(bits, Trellis);
    % m. explain about turbo coding
    % m. we obtain the output, the length of which is twice of the previous
    % m. length: 2 * ceil(8 * length(message) / nfft) * nfft
    % because k/n = 1/2.
    % m. output length is approximately 2 X 8 X (no. of alphabets in
    % message)

    % Finally, let's pre-pend the length to the message
    output = [dec2bin(len, nfft)-'0', output];
    % m. why do we convert the length in 64-bit length?
    % m. because the length info will also be processed in FFT
    % m. this increase the message by 64 bits.
    % m. length: 2 * ceil(8 * length(message) / nfft * nfft) + nfft
    % m. output length is approximately 2 X 8 X (no. of alphabets in
    % message) + 64

end


%% Level #2: Next, lets do interleaving, which permutes the bits
if (level >= 2)
    % Number of symbols in message
    nsym = length(output)/nfft;
    % m. then 4 (= 1/(2 x 8 / 64))chars in message becomes a single symbol
    
    for ii = 1:nsym
        % Collect the iith symbol
        symbol = output((ii-1)*nfft+1:ii*nfft);
        % m. from (ii-1)*nfft+1 to ii*nfft
        
        % Interleave the symbol
        output((ii-1)*nfft+1:ii*nfft) = symbol(Interleave);
        % m. output is a integer array (made of 0s and 1s)
    end
    
    % m. length: (2 * ceil(8 * length(message) / nfft) * nfft) + nfft
    % m. the length has not been changed.
end


%% Level #3: Next, lets do modulation, which maps the bits to a modulation (BPSK)
if (level >= 3)
    % Do BPSK modulation
    output = 2*output-1;
    % m. changes integer 0 into -1 and leaves integer 1 unchanged.
    
    % Prepend a preamble
    output = [preamble, output];
    % m. this increase the length by 64 bits.
    % m. length: (2 * ceil(8 * length(message) / nfft) * nfft) + nfft * 2
end


%% Level #4: Next, lets create an OFDM packet
if (level >= 4)
    % Number of symbols in message
    nsym = length(output)/nfft;
    
    for ii = 1:nsym
        % Collect the iith symbol
        symbol = output((ii-1)*nfft+1:ii*nfft);
        % m. same as the one in iteration within level 2.
        
        % Run an IFFT on the symbol
        % m. IFFT formula has a scaling factor 1/nfft.
        % m. it makes the values too small and be truncated.
        % m. so we multiply nfft so it is not truncated as below.
        % m. but we have to put back 1/nfft in FFT in the receiver site.
        output((ii-1)*nfft+1:ii*nfft) = ifft(symbol)*nfft;
        % m. output becomes a complex number array
    end
    % m. FFT does not change the output length
    % m. length: (2 * ceil(8 * length(message) / nfft) * nfft) + nfft * 2
end


%% Level #5: Finally, lets add some random padding and noise
if (level >= 5)
    % Lets add some (random) empty space to the beginning and end
    % m. to see whether the receiver can locate the preamble
    noise_pad_begin = zeros(1, round(rand*1000));
    % m. rand returns a random number between 0 and 1 exclusively.
    % m. zeros(1,N) returns 1-by-N matrix of zeros.

    noise_pad_end = zeros(1, round(rand*1000));
    
    size(noise_pad_begin)
    size(output)
    size(noise_pad_end)
    % m. length(X) returns the length of vector X.
    % m. this is printed as an output of this 'wifitransmitter' function.
    
    output = [noise_pad_begin, output, noise_pad_end];
    % This adds random numbers to the output.
    
    % Let's add additive white gaussian noise
    output = awgn(output, snr);


    % m. output is complex numbers
    % m. length: random length pad_begin
    % + (2 * ceil(8 * length(message) / nfft) * nfft)
    % + nfft * 2 + random length pad_end
end

end
% m. end of function
