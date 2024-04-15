function [message, length_out, start] = wifireceiver(txsignal, level)
coded_msg = txsignal;
K = 7;
Trellis = poly2trellis(K,[133 171]);
nfft = 64;
ASCII_POWER = 8;
% m. using same as trasmitted.
INTERLEAVE_SZ = 4;
Interleave_undo = reshape(reshape([1:nfft], INTERLEAVE_SZ, []).', [], 1);
% ** to undo interleave, iterate interleaving process as n times.
%       e.q. 1 set : 64 => 16 x 4 // n = 2
%       e.q. 1 set : 16 =>  4 x 4 // n = 1
%       'The n' is iteration number.
syms n; 
iter_num = solve(INTERLEAVE_SZ^n == size(reshape([1:nfft], INTERLEAVE_SZ, []),2),n);
clear n;

% m. preamble is well-known information
preamble = [1, 1, 1, 1,-1,-1, 1, 1,-1, 1,-1, 1, 1, 1, 1, 1, 1,-1,-1, 1, 1,-1, 1,-1, 1, 1, 1, 1, 1,-1,-1, 1, 1,-1, 1,-1, 1,-1,-1,-1,-1,-1, 1, 1,-1, -1, 1,-1, 1,-1, 1, 1, 1, 1,-1,-1, -1,-1,-1, 1, 1,-1, -1, 1];
preamble_sz = size(preamble,2);

%{
figure;
t = linspace(1,size(coded_msg,2),size(coded_msg,2));
plot3(t,real(coded_msg),imag(coded_msg),'k-')
title('ofdm signal(Raw)')
ylabel('real')
zlabel('imag')
%}

%% [Signal noise decoding]
if (level >= 5)
    % m. registers for memory
    pream_ind = 1;              %preamble index.
    correct_num_with_pream = 0; % it needs when we do not find preamble on our coded_msg.

    for ind1 = 1:(size(coded_msg,2)-nfft+1)
        % m. find preamble
        section = coded_msg(ind1:ind1+nfft-1);
        section_value = round(real(fft(section)./nfft));
        sum_TrueBit = sum(section_value == preamble);
        % m. occation : found preamble
        if sum_TrueBit == nfft  % m. excatly found
            pream_ind = ind1;
            break;
        elseif sum_TrueBit > correct_num_with_pream
            pream_ind = ind1;
            correct_num_with_pream = sum_TrueBit;
        end
        if ind1 == (size(coded_msg,2)-preamble_sz+1)
            disp('[Cannot Find Preamble]')
        end
    end
    % m. eliminate preamble and prepand.
    coded_msg = coded_msg(1,pream_ind+preamble_sz:end);
    disp(size(coded_msg))
    coded_msg = coded_msg(1,1:end-mod(size(coded_msg,2),nfft));     % m. set by 64 x n
    disp(size(coded_msg))

    disp(pream_ind)
    start = (pream_ind-1);
else
    start = 0;
end

%{
figure;
t = linspace(1,size(coded_msg,2),size(coded_msg,2));
plot3(t,real(coded_msg),imag(coded_msg),'k-')
title('eliminate preamble')
ylabel('real')
zlabel('imag')
%}


%% [OFDM]   using fft but not domain transformation
if (level >= 4)
    % m. coded_msg is defined at complex plane
    if (level == 4)
        correct_num_with_pream = 0;         %reg. of true_values num
        pream_ind = 1;                      %reg. of preamble index.
        
        for ind1 = 1:(size(coded_msg,2)-nfft+1)
            % m. find preamble
            section = coded_msg(ind1:ind1+nfft-1);
            section_value = round(real(fft(section)./nfft));
            sum_TrueBit = sum(section_value == preamble);
            %disp(sum_TrueBit)
            % m. occation : found preamble
            if sum_TrueBit == nfft  % m. excatly found
                pream_ind = ind1;
                break;
            elseif sum_TrueBit > correct_num_with_pream
                pream_ind = ind1;
                correct_num_with_pream = sum_TrueBit;
            end
            if ind1 == (size(coded_msg,2)-preamble_sz+1)
                disp('[Cannot Find Preamble]')
            end
        end

        % m. when upper loops ends, we could find where preamble started.
        %       which means we could eliminate prepand and preamble in our
        %       coded message.
        %m. now we can got msg without pre-pand and preamble.
        coded_msg = coded_msg(1,pream_ind+preamble_sz:end);
        % m. make our msg + others as number of multiplied by 64.(64의 배수로 msg 조정)
        pre_valid_sz = size(coded_msg,2) - mod(size(coded_msg,2),ASCII_POWER^2);
    
        % m. in level 1, we could eliminate padding of backside with msg length
        %       information
        coded_msg = coded_msg(1:pre_valid_sz);
    end

    %coded_msg
    nsym = length(coded_msg)/nfft;
    %for each sub carrier in 1 set of spectrum
    for ind1 = 1:nsym
        % m. get s_tilda in 64 sample in 1 set of signal
        s_tilda = coded_msg((ind1-1)*nfft+1:ind1*nfft);

        %fft to get symbol.
        coded_msg((ind1-1)*nfft+1:ind1*nfft) = fft(s_tilda)./nfft;
        %disp(coded_msg((ind1-1)*nfft+1:ind1*nfft))
        %coded_msg((ind1-1)*nfft+1:ind1*nfft) = round(coded_msg((ind1-1)*nfft+1:ind1*nfft));
        %{
        
        %disp(coded_msg((ind1-1)*nfft+1:ind1*nfft))

        % m. floating point error counter :
        fpec = 0;
        zero_counter = 0;
        for ind2 = 1:nfft
            % m. eliminate ' imag. ' floating point error 
            if imag(coded_msg(ind2)) > 1e-10
                fpec = fpec + 1;
            end
            if real(coded_msg((ind1-1)*nfft + ind2)) == 0
                zero_counter = zero_counter + 1;
            end
        end
        %{
                disp('---')
        disp(zero_counter)
        disp('---')
        %}
        if fpec == 0
            coded_msg = real(coded_msg);
        end
        %disp('======================[ '+ind1+" ]======================")
        disp(coded_msg((ind1-1)*nfft+1:ind1*nfft))
        %}
    end
    coded_msg = real(coded_msg);
    %{
    figure;
    t = linspace(1,size(coded_msg,2),size(coded_msg,2));
    plot3(t,real(coded_msg),imag(coded_msg))
    title('ofdm demodulation')
    ylabel('real')
    zlabel('imag')
    %}

    % m. leftover of sequence level 52496
    if (level >= 5)
        % m. coded_msg is defined at complex plane
        % m. after ofdm demod, we could our received signal projection to
        %       our symbols. // BPSK
        for ind = 1:size(coded_msg,2)
            if coded_msg(1,ind)>0
                coded_msg(1,ind) = 1;
            else 
                coded_msg(1,ind) = -1;
            end
        end
    else    % m. under level 4...
        coded_msg = round(coded_msg);
    end
    %{
    figure;
    t = linspace(1,size(coded_msg,2),size(coded_msg,2));
    plot3(t,real(coded_msg),imag(coded_msg))
    title('ofdm demodulation')
    ylabel('real')
    zlabel('imag')
    %}
end

%{
figure;
t = linspace(1,size(coded_msg,2),size(coded_msg,2));
plot3(t,real(coded_msg),imag(coded_msg))
title('ofdm demodulation')
ylabel('real')
zlabel('imag')

disp(coded_msg')
size(coded_msg)
coded_msg = round(coded_msg);
disp('-------------')
disp(coded_msg')
disp('-------------')
disp(coded_msg(1,4))
%}
%% Now received bits are binary. [0 or 1]
%% [BPSK demodulation]
if (level >= 3)
    % m. Note that, in out coded msg, it could include zero-padding with
    %   random size, should be eliminated.
    % m. divide preamble part and essential parts.

    % m. lower codes are operated in upper codes.
    if (level == 3)
        correct_num_with_pream = 0;         %reg. of true_values num
        pream_ind = 1;                      %reg. of preamble index.
        
        for ind1 = 1:(size(coded_msg,2)-nfft+1)
            % m. find preamble
            % m. 4,5번과 다르게 ofdm과정이 들어가지 않아 fft가 필요가 없다.
            section_value = coded_msg(ind1:ind1+nfft-1);
            sum_TrueBit = sum(section_value == preamble);
            % m. occation : found preamble
            if sum_TrueBit == nfft  % m. excatly found
                pream_ind = ind1;
                break;
            elseif sum_TrueBit > correct_num_with_pream
                pream_ind = ind1;
                correct_num_with_pream = sum_TrueBit;
            end
            if ind1 == (size(coded_msg,2)-preamble_sz+1)
                disp('[Cannot Find Preamble]')
            end
        end

        % m. when upper loops ends, we could find where preamble started.
        %       which means we could eliminate prepand and preamble in our
        %       coded message.
        %m. now we can got msg without pre-pand and preamble.
        coded_msg = coded_msg(1,pream_ind+preamble_sz:end);
        % m. make our msg + others as number of multiplied by 64.(64의 배수로 msg 조정)
        pre_valid_sz = size(coded_msg,2) - mod(size(coded_msg,2),ASCII_POWER^2);
    
        % m. in level 1, we could eliminate padding of backside with msg length
        %       information
        coded_msg = coded_msg(1:pre_valid_sz);
    end

    % m. temp var for register of (index)Preamble_ind_esti
    
    %{
    for ind1 = 1:(size(coded_msg,2)-preamble_sz+1)
        % coded_msg(1,begin_index:begin_index+preamble_sz-1)
        % preamble
        TrueBit_counter = sum(coded_msg(1,ind1:ind1+preamble_sz-1) == preamble);
        
        if TrueBit_counter > correct_num_with_pream
            % m. when we exactly found our preamble
            if TrueBit_counter == 64
                pream_ind = ind1;% 레벨 4에서 얘로 잡아내는거 보면 여긴 문제가 없다.
                break;
            end

            % m. move present truebit num into register
            correct_num_with_pream = TrueBit_counter;
            % m. move present index of guessed index.
            pream_ind = ind1;
        end
        %Preamble_ind_esti=Preamble_ind_esti+1;
    end
    %disp(size(coded_msg))
    %disp(TrueBit_counter)
    %disp(ind_pream)

    %}
    
    
    %   coded_msg = coded_msg(1,preamble_sz+1:end);
    % m. BPSK demodulation
    coded_msg = (coded_msg+1)./2;
end

%% [interleaving]
if (level >= 2)
    % m. nsym : output data divide into 64 multi sub-carrier.
    %       nsym : index of each sub-carrier.
    % m. reset the form of Interleaving
    nsym = length(coded_msg)/nfft;

    % m. Note that, if we need to undo the interleaving, iterate same
    %    procedure as n times, as I mentioned above.
    for iter = 1:iter_num
       for ind = 1:nsym
            % m. (ind-1)*nfft+1:ind*nfft    indexing
            % m. In 1 set of nsym, undo interleaving.
            symbol = coded_msg((ind-1)*nfft+1:ind*nfft);
            
            % m. Interleave_undo
            coded_msg((ind-1)*nfft+1:ind*nfft) = symbol(Interleave_undo);
        end
    end
end
%% [turbo-decoding]
if (level >= 1)
    % m. len : size of message.
    length_out = bin2dec(string(num2str(coded_msg(1:nfft))));

    % m. get the bit size of msg, by length_out
    %       could take only msg.
    msg_bit_size = length_out*8;
    msg_bit_size = 2*(msg_bit_size + mod(-msg_bit_size,nfft));
    % m. remove pre-pand block which include the info. of message len.
    
    % m. cut... 
    %   | preamble | msg; (nfft+1) ~ nfft+msg_bit_size | ...> | padding at the end.
    coded_msg = coded_msg(1,(nfft+1):(nfft)+msg_bit_size);

    % m. constraint length is optimal between 3K ~ 5K
    bin_message = erase(num2str(vitdec(coded_msg, Trellis,3*K,'trunc','hard')),' ');
    % ASCII have 2^8 assigned character...
    % m. should cut binary stream by 8 bits.
    
    bin_message = reshape(bin_message',ASCII_POWER,size(bin_message,2)/ASCII_POWER)';

    % m. If there is some left over blanks; remove blanks at behinds (filled as '00000000')
    while(1)
        if sum(bin_message(end,:) == '00000000') == 8
            bin_message = bin_message(1:end-1,:);
        else
            % m. if bin_message ~= '00000000', escape.
            break;
        end
    end
    
    % m. change into message string type.
    temp_message = '';
    for i = 1:size(bin_message,1)
        % m. transform 8 bits data stream into char with ASCII
        temp_ASCII = char(bin2dec(bin_message(i, :)));
        % m. accumulation char
        temp_message = [temp_message temp_ASCII];
    end
    message = temp_message;
end

end