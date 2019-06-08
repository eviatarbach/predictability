function [] = granger_multiple(offset, length)

nobs = 14235;
max_delay = 360;
max_order = 45;
delays = 0:max_delay;

cd mvgc_v1.0

startup

cd ..

sst = importdata('data/sst01_365.mat');
vort = importdata('data/vort01_365.mat');
div = importdata('data/div01_365.mat');
sp = importdata('data/sp01_365.mat');
temp = importdata('data/temp01_365.mat');
q = importdata('data/q01_365.mat');

F_atmos_to_sst = NaN(length, max_delay + 1);
F_sst_to_atmos = NaN(length, max_delay + 1);

mspe_full_atmos_to_sst = NaN(length, max_delay + 1);
mspe_reduced_atmos_to_sst = NaN(length, max_delay + 1);
mspe_full_sst_to_atmos = NaN(length, max_delay + 1);
mspe_reduced_sst_to_atmos = NaN(length, max_delay + 1);

sig_atmos_to_sst = NaN(length, max_delay + 1);
sig_sst_to_atmos = NaN(length, max_delay + 1);

times = NaN(length);

'offset', offset

j = 1;

for delay = delays
    F_atmos_to_sst_j = F_atmos_to_sst(:, j);
    F_sst_to_atmos_j = F_sst_to_atmos(:, j);
    mspe_full_atmos_to_sst_j = mspe_full_atmos_to_sst(:, j);
    mspe_reduced_atmos_to_sst_j = mspe_reduced_atmos_to_sst(:, j);
    mspe_full_sst_to_atmos_j = mspe_full_sst_to_atmos(:, j);
    mspe_reduced_sst_to_atmos_j = mspe_reduced_sst_to_atmos(:, j);
    sig_atmos_to_sst_j = sig_atmos_to_sst(:, j);
    sig_sst_to_atmos_j = sig_sst_to_atmos(:, j);

    for cell = 1:length
        if offset+cell-1 > 88838
            continue
        end
        sst_ts = sst(:, offset+cell-1)';
        vort_ts = vort(:, offset+cell-1)';
        div_ts = div(:, offset+cell-1)';
        sp_ts = sp(:, offset+cell-1)';
        temp_ts = temp(:, offset+cell-1)';
        q_ts = q(:, offset+cell-1)';

        if (sst_ts(1) == 9999) || (vort_ts(1) == 9999) || (div_ts(1) == 9999) || (sp_ts(1) == 9999) || (temp_ts(1) == 9999) || (q_ts(1) == 9999)
            continue
        end

        if delay == 0
            [AIC, ~, ~, ~] = tsdata_to_infocrit([sst_ts; vort_ts; div_ts; sp_ts; temp_ts; q_ts]);
            [~, idx] = sort(AIC(2:end));
            moAIC = idx(1) + 1;
            times(cell) = moAIC;
        else
            moAIC = times(cell);
        end

        sst_ts_effect = sst_ts(1+delay:end);
        sst_ts_cause = sst_ts(1:end-delay);

        vort_ts_effect = vort_ts(1+delay:end);
        vort_ts_cause = vort_ts(1:end-delay);

        div_ts_effect = div_ts(1+delay:end);
        div_ts_cause = div_ts(1:end-delay);

        sp_ts_effect = sp_ts(1+delay:end);
        sp_ts_cause = sp_ts(1:end-delay);

        temp_ts_effect = temp_ts(1+delay:end);
        temp_ts_cause = temp_ts(1:end-delay);

        q_ts_effect = q_ts(1+delay:end);
        q_ts_cause = q_ts(1:end-delay);

        X_sst_cause = [sst_ts_cause; vort_ts_effect; div_ts_effect; sp_ts_effect; temp_ts_effect; q_ts_effect];
        X_atmos_cause = [sst_ts_effect; vort_ts_cause; div_ts_cause; sp_ts_cause; temp_ts_cause; q_ts_cause];

        [A_sst_cause, SIG_sst_cause] = tsdata_to_var(X_sst_cause, moAIC);
        assert(~isbad(A_sst_cause), 'VAR estimation failed');
        [G_sst_cause, info_sst_cause] = var_to_autocov(A_sst_cause, SIG_sst_cause);
        acerr_sst_cause = var_info(info_sst_cause, false);

        [A_atmos_cause, SIG_atmos_cause] = tsdata_to_var(X_atmos_cause, moAIC);
        assert(~isbad(A_atmos_cause), 'VAR estimation failed');
        [G_atmos_cause, info_atmos_cause] = var_to_autocov(A_atmos_cause, SIG_atmos_cause);
        acerr_atmos_cause = var_info(info_atmos_cause, false);

        if acerr_sst_cause || acerr_atmos_cause
            continue
        end

        [F_atmos_cause, mspe_full_atmos_cause, mspe_reduced_atmos_cause] = autocov_to_mvgc(G_atmos_cause, [1], [2 3 4 5 6]);  % atmos -> sst

        [F_sst_cause, mspe_full_sst_cause, mspe_reduced_sst_cause] = autocov_to_mvgc(G_sst_cause, [2 3 4 5 6], [1]);  % sst -> atmos

        F_atmos_to_sst_ij(cell) = F_atmos_cause;
        F_sst_to_atmos_ij(cell) = F_sst_cause;

        if isreal(F_atmos_cause)
            mspe_full_atmos_to_sst_j(cell) = mspe_full_atmos_cause;
            mspe_reduced_atmos_to_sst_j(cell) = mspe_reduced_atmos_cause;
            pval_atmos_cause = mvgc_pval(F_atmos_cause, moAIC, nobs - delay, 1, 1, 5, 0);
            sig_atmos_to_sst_j(cell) = pval_atmos_cause;
        end
        if isreal(F_sst_cause)
            mspe_full_sst_to_atmos_j(cell) = mspe_full_sst_cause;
            mspe_reduced_sst_to_atmos_j(cell) = mspe_reduced_sst_cause;
            pval_sst_cause = mvgc_pval(F_sst_cause, moAIC, nobs - delay, 1, 5, 1, 0);
            sig_sst_to_atmos_j(cell) = pval_sst_cause;
        end
    end
    F_atmos_to_sst(:, j) = F_atmos_to_sst_j;
    F_sst_to_atmos(:, j) = F_sst_to_atmos_j;
    mspe_full_atmos_to_sst(:, j) = mspe_full_atmos_to_sst_j;
    mspe_reduced_atmos_to_sst(:, j) = mspe_reduced_atmos_to_sst_j;
    mspe_full_sst_to_atmos(:, j) = mspe_full_sst_to_atmos_j;
    mspe_reduced_sst_to_atmos(:, j) = mspe_reduced_sst_to_atmos_j;
    sig_atmos_to_sst(:, j) = sig_atmos_to_sst_j;
    sig_sst_to_atmos(:, j) = sig_sst_to_atmos_j;

    j = j + 1;
end

save(['data_atmos/F_atmos_to_sst_' num2str(offset) '.mat'], 'F_atmos_to_sst');
save(['data_atmos/F_sst_to_atmos_' num2str(offset) '.mat'], 'F_sst_to_atmos');
save(['data_atmos/mspe_full_atmos_to_sst_' num2str(offset) '.mat'], 'mspe_full_atmos_to_sst');
save(['data_atmos/mspe_reduced_atmos_to_sst_' num2str(offset) '.mat'], 'mspe_reduced_atmos_to_sst');
save(['data_atmos/mspe_full_sst_to_atmos_' num2str(offset) '.mat'], 'mspe_full_sst_to_atmos');
save(['data_atmos/mspe_reduced_sst_to_atmos_' num2str(offset) '.mat'], 'mspe_reduced_sst_to_atmos');
save(['data_atmos/sig_atmos_to_sst_' num2str(offset) '.mat'], 'sig_atmos_to_sst');
save(['data_atmos/sig_sst_to_atmos_' num2str(offset) '.mat'], 'sig_sst_to_atmos');
save(['data_atmos/times_atmos_' num2str(offset) '.mat'], 'times');
