function [dec_bits, success, iter] = LDPC_decoder(llr_channel, H, maxiter)
% llr_channel:column vector
% H:sparse matrix
% maxiter:scalar

bitlen = length(llr_channel);
[row, column] = size(H);

% i/j:row/culumn index of check node to calculate for each variable node
% max(i_c) <= row
[i_v, j_v, ~] = find(H);

% i/j:column/row index of variable node to calculate for each check node
% max(i_v) <= column
[j_c, i_c, ~] = find(H');

% initial variable and check node matrix
llr_mat = zeros(row, column);

% record if the input code is successfully decoded
success = 0;

% updated llr in each iter to decode
L_total = zeros(bitlen, 1);

% decode process
for iter = 1:maxiter

    % variable node process
    % for idx_v = 1:length(i_v)
    %     if((idx_v == 1) ||  (j_v(idx_v-1) ~= j_v(idx_v)))
    %         L_v = sum(llr_mat(:, j_v(idx_v)).*H(:, j_v(idx_v)), 1);
    %         L_total(j_v(idx_v)) = L_v + llr_channel(j_v(idx_v));
    %     end
    %     % calculate new llr of the variable node
    %     llr_mat(i_v(idx_v), j_v(idx_v)) = L_total(j_v(idx_v)) - llr_mat(i_v(idx_v), j_v(idx_v));
    % end

    % variable node process with matrix process
    L_total = sum(llr_mat, 1).' + llr_channel;
    llr_mat_sum = H*spdiags(L_total, 0, column, column);
    llr_mat = llr_mat_sum - llr_mat;
    
    sign_all = 0;
    min_0 = 0;
    min_1 = 0;
    idx_min = 0;

    % check node process
    for idx_c = 1:length(i_c)
        if((idx_c == 1) || (i_c(idx_c-1) ~= i_c(idx_c)))
            % sign of all variable node linked to the check node i_c(idx_c)
            sign_all = -(2*mod(sum(llr_mat(i_c(idx_c), :) < 0), 2)-1);
            % find min and second min llr of a check node
            [val, idx] = sort(abs(llr_mat(i_c(idx_c), :))); 
            magmins_idx = find(val, 2);
            mins = val(magmins_idx);
            minidxs = idx(magmins_idx);
            if(length(mins) > 1)
                idx_min = minidxs(1);
                min_0 = mins(1);
                min_1 = mins(2);
            elseif(length(mins) == 1)
                idx_min = minidxs(1);
                min_0 = mins(1);
                min_1 = 0;
            end         
        end
        % calculate the llr of the check node
        if(j_c(idx_c) ~= idx_min)
            if(llr_mat(i_c(idx_c), j_c(idx_c)) ~= 0)
                llr_mat(i_c(idx_c), j_c(idx_c)) = (2*(sign_all*llr_mat(i_c(idx_c), j_c(idx_c)) > 0)-1)*min_0;
            elseif(llr_mat(i_c(idx_c), j_c(idx_c)) == 0)
                llr_mat(i_c(idx_c), j_c(idx_c)) = sign*min_0;
            end
        else
            if(llr_mat(i_c(idx_c), j_c(idx_c)) ~= 0)
                llr_mat(i_c(idx_c), j_c(idx_c)) = (2*(sign_all*llr_mat(i_c(idx_c), j_c(idx_c)) > 0)-1)*min_1;
            elseif(llr_mat(i_c(idx_c), j_c(idx_c)) == 0)
                llr_mat(i_c(idx_c), j_c(idx_c)) = sign*min_1;
            end
        end

    end
    
    % process for each check node(row in H)
    % for idx_c = 1:row
    %     node_c = llr_mat(idx_c, :);
    %     sign_all = -(2*mod(sum(node_c < 0), 2)-1);
    %     [val, idx] = sort(abs(node_c)); 
    %     nonzero_idx = find(val);
    %     if(length(nonzero_idx) == 1)
    %         min_1 = 0;
    %         min_0 = val(idx(nonzero_idx(1)));
    %         idx_min = idx(nonzero_idx(1));
    %     else
    %         min_1 = val(idx(nonzero_idx(2)));
    %         min_0 = val(idx(nonzero_idx(1)));
    %         idx_min = idx(nonzero_idx(1));
    %     end
    %     sign_bit = sign(node_c);
    %     check_one = find(H(idx_c, :));
    %     one_in_one = find(sign_bit(check_one));
    %     zero_in_one = check_one(~ismember(check_one, one_in_one));
    %     if(~isempty(zero_in_one))
    %         sign_bit(zero_in_one) = 1;
    %     end
    %     sign_bit_c = sign_bit*sign_all;
    %     sign_bit_c(idx_min) = sign_bit_c(idx_min)*min_1;
    %     idx_exmin = check_one(~ismember(check_one, idx_min));
    %     sign_bit_c(idx_exmin) = sign_bit_c(idx_exmin)*min_0;
    %     llr_mat(idx_c, :) = sign_bit_c;
    % end

    decoded_bits = (L_total <= 0);
    check = sum(mod(H * decoded_bits, 2));
    if(check == 0)
        success = 1;
        dec_bits = decoded_bits(1:column-row);
        break;
    else
        success = 0;
        dec_bits = decoded_bits(1:column-row);
    end

end

end