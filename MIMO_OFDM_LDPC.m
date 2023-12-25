% ͨ��ϵͳ������SOC���� 2023���＾ѧ��
clear;close all;clc;

load('H_3328_8960.mat');

% SNR��λΪdB�������е���
SNR_list = 0 : 2 : 18;

% ÿ��SNR���������OFDM���������������е���
block_num = 500;

% FFT����
N_fft = 1024;     

% ���ƽ�����4(QPSK)��16(16QAM)
Q = 2;

% ÿ�����ų��صı�����
B = log2(Q);

% ����ؼ���
BER_count = zeros(length(SNR_list), block_num);

% ������Ϣ�����ز�����
N_sc = 896;

% ���䣬������������
N_t = 3;
N_r = 3;

% ÿ��OFDM����������ı�����
N_bit = N_sc * B * N_t;

% bit num of LDPC code
[row, column] = size(ldpc_H);
N_ldpc = column;
N_info = column-row;
H_sparse = sparse(logical(ldpc_H));
enccfg = ldpcEncoderConfig(H_sparse);
block = ceil(N_ldpc/N_bit);
zeros_num = block*N_bit-N_ldpc;
iter_num = 50;

% CP����
length_CP = 73;

% 0:ZF 1:MMSE
equal_method = 0;

% ��ȡ��Ч�����ŵ��弤��Ӧhij������i=1,2,3��j=1,2,3
% hij��ʾ�ӵ�j���������ߵ���i���������ߵĵ�Ч�����ŵ��弤��Ӧ
load('h.mat');

for snr_count = 1 : length(SNR_list)
    
    snr = SNR_list(snr_count);
    
    % ���￼�Ƿ�������֮��ĵȹ��ʷ��䣬��ÿ������������ÿ���������ز��Ϸ�����ŵĹ���Ϊ1 / Nt
    % ��������������и��������ߵ�Ƶ������ŵ�ƽ������֮��Ϊ1��SNR�Ķ���Ϊ��N_sc / (N_fft ^ 2 * ��������)
    sigma2 = N_sc / (10 ^ (snr / 10) * N_fft ^ 2);
    
    parfor count = 1 : block_num
        snr_count, count
        
        %������Ϣ���� + LDPC coding
        msg = round(rand(1, N_info));
        ldpc_msg = ldpcEncode(msg.', enccfg).';
        
        %% MIMO-OFDM���Ʋ���
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % �������Ĵ��룬��ɴӱ��ص�3·����������������QPSK��16QAM��������ӳ��,��Ƶ����ɻ���SVD��Ԥ���룬����OFDM���Ʋ����CP
        % ���룺��Ϣ��������msg
        % �����3·OFDM���ƺ����CP��ʱ���ź�����x_ofdm
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % qpsk/16qam modulation
        % if(Q == 4)
        %     msg_mod = (1/sqrt(N_t))*pskmod(bit2int(reshape(msg, B, []), B), Q);
        % elseif(Q == 16)
        %     msg_mod = (1/sqrt(N_t))*(1/sqrt(10))*qammod(bit2int(reshape(msg, B, []), B), Q);
        % end

        % BPSK modulation 0 -> 1 1 -> -1 + zero padding
        msg_mod = -(ldpc_msg*2-1);
        msg_mod = [msg_mod, zeros(1, zeros_num)];
        
        % divide data into 3 data streams
        msg_3_tr = pagetranspose(reshape(msg_mod, [], 3, ceil(N_ldpc/N_bit)));
        msg_3_tr_zeros = [msg_3_tr, zeros(N_t, N_fft-N_sc, block)];
        
        % channel matrix
        H = zeros(N_r, N_t, N_fft);
        H(1, 1, :) = fft(h11, N_fft);
        H(1, 2, :) = fft(h12, N_fft);
        H(1, 3, :) = fft(h13, N_fft);
        H(2, 1, :) = fft(h21, N_fft);
        H(2, 2, :) = fft(h22, N_fft);
        H(2, 3, :) = fft(h23, N_fft);
        H(3, 1, :) = fft(h31, N_fft);
        H(3, 2, :) = fft(h32, N_fft);
        H(3, 3, :) = fft(h33, N_fft);

        % precoding
        U = zeros(N_r, N_r, N_fft);
        S = zeros(N_r, N_r, N_fft);
        msg_prec = zeros(N_r, N_fft, block);

        for i = 1:N_fft
            [U(:, :, i), S(:, :, i), V] = svd(H(:, :, i));
            msg_prec(:, i, :) = pagemtimes(V, msg_3_tr_zeros(:, i, :));
        end
        
        % OFDM + D to S
        msg_ofdm = ifft(msg_prec, [], 2);
        
        % cyclic prefix
        x_cyc = [msg_ofdm(:, length(msg_ofdm)-length_CP+1:end, :), msg_ofdm]; 
        x_ofdm = [];
        for idx = 1:block
            x_ofdm = [x_ofdm, x_cyc(:, :, idx)];
        end

        %% �ŵ����䲿��
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % �������Ĵ��룬�������x_ofdm�����ྶ�ŵ�hij������նˣ������AWGN��ȥ��CP��ע��ʵ�鲿�������ʸ�Ϊ�������ʵ�һ��
        % �����Ѿ������CP����������һ��OFDM���ŶԱ�OFDM���ŵ�Ӱ��
        % ���룺3·���CP���ʱ��OFDM�������x_ofdm
        % �����3·����ΪN_fft�Ľ���OFDM����r_ofdm
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % r_1 = conv(x_ofdm(1, :), h11, 'same') + conv(x_ofdm(2, :), h12, 'same') + conv(x_ofdm(3, :), h13, 'same');
        % r_2 = conv(x_ofdm(1, :), h21, 'same') + conv(x_ofdm(2, :), h22, 'same') + conv(x_ofdm(3, :), h23, 'same');
        % r_3 = conv(x_ofdm(1, :), h31, 'same') + conv(x_ofdm(2, :), h32, 'same') + conv(x_ofdm(3, :), h33, 'same');
        
        r_1 = conv(x_ofdm(1, :), h11) + conv(x_ofdm(2, :), h12) + conv(x_ofdm(3, :), h13);
        r_2 = conv(x_ofdm(1, :), h21) + conv(x_ofdm(2, :), h22) + conv(x_ofdm(3, :), h23);
        r_3 = conv(x_ofdm(1, :), h31) + conv(x_ofdm(2, :), h32) + conv(x_ofdm(3, :), h33);

        r_time = [r_1; r_2; r_3];      

        r_channel = r_time + randn(N_r, length(r_1)) * sqrt(0.5 * sigma2) * (1+1i); 
        r_block = zeros(N_r, N_fft+length_CP, block);
        for i = 1:block
            r_block(:, :, i) = r_channel(:, (i-1)*(N_fft+length_CP) + 1:i*(N_fft+length_CP));
        end
        r_ofdm = r_block(:, length_CP+1:length_CP+N_fft, :);
        
        %% MIMO-OFDM�������       
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % ��������Լ��Ĵ��룬�Խ��յ���OFDM-MIMO���Ž��������MMSE������о������㲢�ָ��������
        % ���룺3·����OFDM����r_ofdm
        % ����������ı�������msg_r
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        r_fft = fft(r_ofdm, [], 2);

        r_com = zeros(N_r, N_fft, block);
        r_eq = zeros(N_r, N_fft, block);
        % combine
        for i = 1:N_fft
            r_com(:, i, :) = pagemtimes(U(:, :, i)', r_fft(:, i, :));        
            % equalization
            if(equal_method == 0)
                r_eq(:, i, :) = pagemtimes(inv(S(:, :, i)), r_com(:, i, :));
            elseif(equal_method == 1)
                r_eq(:, i, :) = pagemtimes(inv((S(:, :, i)' * S(:, :, i) + sigma2 * eye(N_t))), pagemtimes(S(:, :, i)', r_com(:, i)));
            end
        end
        
        % delete zeros
        r_data = sqrt(N_t) * reshape(pagetranspose(r_eq(:, 1:N_sc, :)), [], 1).';
        r_data = r_data(1:N_ldpc);
        
        % demod
        % if(Q == 4)
        %     msg_demod = pskdemod(r_data, Q);
        % elseif(Q == 16)
        %     msg_demod = qamdemod(sqrt(10) * r_data, Q);
        % end
        
        % calculate llr of BPSK
        % llr_data = -(1/sigma2)*(abs(r_data - 1).^2 - abs(r_data + 1).^2);
        llr_data = 2*real(r_data)/sigma2;
        % puncture
        llr_data = [zeros(1, 512), llr_data(513:end)];

        % LDPC decoding
        [msg_r, success, iter] = LDPC_decoder(llr_data.', H_sparse, iter_num);
        
        % msg_r = reshape(int2bit(msg_demod, B), 1, []);
        
        %% �������ͳ��
        BER_count(snr_count, count) = BER_count(snr_count, count) + sum(abs(msg_r.' - msg));
        
    end
end

% �����ʼ���
BER = BER_count / (N_info);
BER = sum(BER, 2) / block_num;

%% BER���߻���
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% �������Ĵ��룬���ð�������꣬ʹ��semilogy��������BER-SNR����
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
semilogy(SNR_list, BER, 'LineWidth', 1.5);
xlabel('SNR/dB');
ylabel('BER'); 
grid on;

