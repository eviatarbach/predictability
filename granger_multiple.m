function [] = granger_multiple(offset, length)
cd /lustre/ebach/causality

nobs = {14235, 2847, 949};
max_delay = {360, 180, 120};
max_order = {45, 9, 3};
delays = {0:max_delay{1}, 0:max_delay{2}, 0:max_delay{3}};

cd mvgc_v1.0

startup

cd ..

sst_05 = importdata('data/sst05_365.mat');
vort_05 = importdata('data/vort05_365.mat');
div_05 = importdata('data/div05_365.mat');
sp_05 = importdata('data/sp05_365.mat');
temp_05 = importdata('data/temp05_365.mat');
q_05 = importdata('data/q05_365.mat');

sst_15 = importdata('data/sst15_365.mat');
vort_15 = importdata('data/vort15_365.mat');
div_15 = importdata('data/div15_365.mat');
sp_15 = importdata('data/sp15_365.mat');
temp_15 = importdata('data/temp15_365.mat');
q_15 = importdata('data/q15_365.mat');

sst = {sst_05, sst_05, sst_15};
vort = {vort_05, vort_05, vort_15};
div = {div_05, div_05, div_15};
sp = {sp_05, sp_05, sp_15};
temp = {temp_05, temp_05, temp_15};
q = {q_05, q_05, q_15};

F_atmos_to_sst = {NaN(length, max_delay{1} + 1), NaN(length, max_delay{2} + 1), NaN(length, max_delay{3} + 1)};
F_sst_to_atmos = {NaN(length, max_delay{1} + 1), NaN(length, max_delay{2} + 1), NaN(length, max_delay{3} + 1)};

mspe_full_atmos_to_sst = {NaN(length, max_delay{1} + 1), NaN(length, max_delay{2} + 1), NaN(length, max_delay{3} + 1)};
mspe_reduced_atmos_to_sst = {NaN(length, max_delay{1} + 1), NaN(length, max_delay{2} + 1), NaN(length, max_delay{3} + 1)};
mspe_full_sst_to_atmos = {NaN(length, max_delay{1} + 1), NaN(length, max_delay{2} + 1), NaN(length, max_delay{3} + 1)};
mspe_reduced_sst_to_atmos = {NaN(length, max_delay{1} + 1), NaN(length, max_delay{2} + 1), NaN(length, max_delay{3} + 1)};

sig_atmos_to_sst = {NaN(length, max_delay{1} + 1), NaN(length, max_delay{2} + 1), NaN(length, max_delay{3} + 1)};
sig_sst_to_atmos = {NaN(length, max_delay{1} + 1), NaN(length, max_delay{2} + 1), NaN(length, max_delay{3} + 1)};

times = {NaN(length), NaN(length), NaN(length)};
%times = importdata('data/times.mat');

'offset', offset

for i = 2:3
    nobs_i = nobs{i};
    sst_i = sst{i};
    vort_i = vort{i};
    div_i = div{i};
    sp_i = sp{i};
    temp_i = temp{i};
    q_i = q{i};
    delays_i = delays{i};

    times_i = times{i};
    F_atmos_to_sst_i = F_atmos_to_sst{i};
    F_sst_to_atmos_i = F_sst_to_atmos{i};
    mspe_full_atmos_to_sst_i = mspe_full_atmos_to_sst{i};
    mspe_reduced_atmos_to_sst_i = mspe_reduced_atmos_to_sst{i};
    mspe_full_sst_to_atmos_i = mspe_full_sst_to_atmos{i};
    mspe_reduced_sst_to_atmos_i = mspe_reduced_sst_to_atmos{i};
    sig_sst_to_atmos_i = sig_sst_to_atmos{i};
    sig_atmos_to_sst_i = sig_atmos_to_sst{i};

    j = 1;

    for delay = delays_i
        F_atmos_to_sst_ij = F_atmos_to_sst_i(:, j);
        F_sst_to_atmos_ij = F_sst_to_atmos_i(:, j);
        mspe_full_atmos_to_sst_ij = mspe_full_atmos_to_sst_i(:, j);
        mspe_reduced_atmos_to_sst_ij = mspe_reduced_atmos_to_sst_i(:, j);
        mspe_full_sst_to_atmos_ij = mspe_full_sst_to_atmos_i(:, j);
        mspe_reduced_sst_to_atmos_ij = mspe_reduced_sst_to_atmos_i(:, j);
        sig_atmos_to_sst_ij = sig_atmos_to_sst_i(:, j);
        sig_sst_to_atmos_ij = sig_sst_to_atmos_i(:, j);

        for cell = 1:length
            if offset+cell-1 > 88838
                continue
            end
            sst_ts = sst_i(:, offset+cell-1)';
            vort_ts = vort_i(:, offset+cell-1)';
            div_ts = div_i(:, offset+cell-1)';
            sp_ts = sp_i(:, offset+cell-1)';
            temp_ts = temp_i(:, offset+cell-1)';
            q_ts = q_i(:, offset+cell-1)';

            if (sst_ts(1) == 9999) || (vort_ts(1) == 9999) || (div_ts(1) == 9999) || (sp_ts(1) == 9999) || (temp_ts(1) == 9999) || (q_ts(1) == 9999)
                continue
            end

            % sst_ts = detrend(sst_ts);  % remove global warming signal

            if delay == 0
                [AIC, ~, ~, ~] = tsdata_to_infocrit([sst_ts; vort_ts; div_ts; sp_ts; temp_ts; q_ts], max_order{i});
                [~, idx] = sort(AIC(2:end));
                moAIC = idx(1) + 1;
                times_i(cell) = moAIC;
            else
                moAIC = times_i(cell);
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
                mspe_full_atmos_to_sst_ij(cell) = mspe_full_atmos_cause;
                mspe_reduced_atmos_to_sst_ij(cell) = mspe_reduced_atmos_cause;
                pval_atmos_cause = mvgc_pval(F_atmos_cause, moAIC, nobs_i - delay, 1, 1, 5, 0);
                sig_atmos_to_sst_ij(cell) = pval_atmos_cause;
            end
            if isreal(F_sst_cause)
                mspe_full_sst_to_atmos_ij(cell) = mspe_full_sst_cause;
                mspe_reduced_sst_to_atmos_ij(cell) = mspe_reduced_sst_cause;
                pval_sst_cause = mvgc_pval(F_sst_cause, moAIC, nobs_i - delay, 1, 5, 1, 0);
                sig_sst_to_atmos_ij(cell) = pval_sst_cause;
            end
        end
        F_atmos_to_sst_i(:, j) = F_atmos_to_sst_ij;
        F_sst_to_atmos_i(:, j) = F_sst_to_atmos_ij;
        mspe_full_atmos_to_sst_i(:, j) = mspe_full_atmos_to_sst_ij;
        mspe_reduced_atmos_to_sst_i(:, j) = mspe_reduced_atmos_to_sst_ij;
        mspe_full_sst_to_atmos_i(:, j) = mspe_full_sst_to_atmos_ij;
        mspe_reduced_sst_to_atmos_i(:, j) = mspe_reduced_sst_to_atmos_ij;
        sig_atmos_to_sst_i(:, j) = sig_atmos_to_sst_ij;
        sig_sst_to_atmos_i(:, j) = sig_sst_to_atmos_ij;

        j = j + 1;
    end
    F_atmos_to_sst{i} = F_atmos_to_sst_i;
    F_sst_to_atmos{i} = F_sst_to_atmos_i;
    mspe_full_atmos_to_sst{i} = mspe_full_atmos_to_sst_i;
    mspe_reduced_atmos_to_sst{i} = mspe_reduced_atmos_to_sst_i;
    mspe_full_sst_to_atmos{i} = mspe_full_sst_to_atmos_i;
    mspe_reduced_sst_to_atmos{i} = mspe_reduced_sst_to_atmos_i;
    sig_sst_to_atmos{i} = sig_sst_to_atmos_i;
    sig_atmos_to_sst{i} = sig_atmos_to_sst_i;
    times{i} = times_i;
end

save(['data_atmos3/F_atmos_to_sst_' num2str(offset) '.mat'], 'F_atmos_to_sst');
save(['data_atmos3/F_sst_to_atmos_' num2str(offset) '.mat'], 'F_sst_to_atmos');
save(['data_atmos3/mspe_full_atmos_to_sst_' num2str(offset) '.mat'], 'mspe_full_atmos_to_sst');
save(['data_atmos3/mspe_reduced_atmos_to_sst_' num2str(offset) '.mat'], 'mspe_reduced_atmos_to_sst');
save(['data_atmos3/mspe_full_sst_to_atmos_' num2str(offset) '.mat'], 'mspe_full_sst_to_atmos');
save(['data_atmos3/mspe_reduced_sst_to_atmos_' num2str(offset) '.mat'], 'mspe_reduced_sst_to_atmos');
save(['data_atmos3/sig_atmos_to_sst_' num2str(offset) '.mat'], 'sig_atmos_to_sst');
save(['data_atmos3/sig_sst_to_atmos_' num2str(offset) '.mat'], 'sig_sst_to_atmos');
save(['data_atmos3/times_atmos_' num2str(offset) '.mat'], 'times');
