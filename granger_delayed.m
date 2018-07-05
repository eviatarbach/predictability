function granger_delayed(lon, lat)
lat = 145;
lon = 245;

nobs = {13870, 2774, 924};
max_delay = {30, 6, 2};
delays = {0:max_delay{1}, 0:max_delay{2}, 0:max_delay{3}};
mhtc = 'FDR';
alpha = 0.10;

cd mvgc_v1.0

startup

cd ..

sst_01 = importdata('data/sst_01.mat');
vort_01 = importdata('data/vort_01.mat');
sst_05 = importdata('data/sst_05.mat');
vort_05 = importdata('data/vort_05.mat');
sst_15 = importdata('data/sst_15.mat');
vort_15 = importdata('data/vort_15.mat');

'loaded'

sst = {sst_01, sst_05, sst_15};
vort = {vort_01, vort_05, vort_15};

F_vort_to_sst = {NaN(1, max_delay{1} + 1), NaN(1, max_delay{2} + 1), NaN(1, max_delay{3} + 1)};
F_sst_to_vort = {NaN(1, max_delay{1} + 1), NaN(1, max_delay{2} + 1), NaN(1, max_delay{3} + 1)};
sig_vort_to_sst = {NaN(1, max_delay{1} + 1), NaN(1, max_delay{2} + 1), NaN(1, max_delay{3} + 1)};
sig_sst_to_vort = {NaN(1, max_delay{1} + 1), NaN(1, max_delay{2} + 1), NaN(1, max_delay{3} + 1)};

times = importdata('data/times.mat');

for i = 1:3
    nobs_i = nobs{i};
    sst_i = sst{i};
    vort_i = vort{i};
    delays_i = delays{i};

    times_i = times{i};
    F_vort_to_sst_i = F_vort_to_sst{i};
    F_sst_to_vort_i = F_sst_to_vort{i};

    sig_sst_to_vort_i = sig_sst_to_vort{i};
    sig_vort_to_sst_i = sig_vort_to_sst{i};

    pval_sst_cause = NaN(1, max_delay{i} + 1);
    pval_vort_cause = NaN(1, max_delay{i} + 1);

    parfor lon=1:480
        for lat=1:241
            for delay = delays_i
                sst_ts = reshape(sst_i(lon, lat, :), [1, nobs_i]);
                vort_ts = reshape(vort_i(lon, lat, :), [1, nobs_i]);

                if (sst_ts(1) == -9999) || (vort_ts(1) == -9999)
                    continue
                end

                sst_ts = detrend(sst_ts);  % remove global warming signal

                moAIC = times_i(lon, lat);

                sst_ts_effect = sst_ts(1+delay:end);
                sst_ts_cause = sst_ts(1:end-delay);

                vort_ts_effect = vort_ts(1+delay:end);
                vort_ts_cause = vort_ts(1:end-delay);

                X_sst_cause = [sst_ts_cause; vort_ts_effect];
                X_vort_cause = [sst_ts_effect; vort_ts_cause];

                [A_sst_cause, SIG_sst_cause] = tsdata_to_var(X_sst_cause, moAIC);
                assert(~isbad(A_sst_cause), 'VAR estimation failed');
                [G_sst_cause, info_sst_cause] = var_to_autocov(A_sst_cause, SIG_sst_cause);
                acerr_sst_cause = var_info(info_sst_cause, false);

                [A_vort_cause, SIG_vort_cause] = tsdata_to_var(X_vort_cause, moAIC);
                assert(~isbad(A_vort_cause), 'VAR estimation failed');
                [G_vort_cause, info_vort_cause] = var_to_autocov(A_vort_cause, SIG_vort_cause);
                acerr_vort_cause = var_info(info_vort_cause, false);

                if acerr_sst_cause || acerr_vort_cause
                    continue
                end

                F_1 = autocov_to_mvgc(G_vort_cause, [1], [2]);  % vort -> sst
                F_2 = autocov_to_mvgc(G_sst_cause, [2], [1]);  % sst -> vort

                F_vort_to_sst_i(delay + 1) = F_1;
                F_sst_to_vort_i(delay + 1) = F_2;

                pval_vort_cause(delay + 1) = mvgc_pval(F_1, moAIC, nobs_i - delay, 1, 1, 1, 0);
                pval_sst_cause(delay + 1) = mvgc_pval(F_2, moAIC, nobs_i - delay, 1, 1, 1, 0);
            end

        %    sig = significance([pval_vort_cause pval_sst_cause], alpha, mhtc);
        %    sig_vort_cause = sig(1:max_delay{i} + 1)
        %    sig_sst_cause = sig(max_delay{i} + 2:end)
            sig_vort_cause = pval_vort_cause;
            sig_sst_cause = pval_sst_cause;
            sig_vort_to_sst_i = sig_vort_cause;
            sig_sst_to_vort_i = sig_sst_cause;

            F_vort_to_sst{i} = F_vort_to_sst_i;
            F_sst_to_vort{i} = F_sst_to_vort_i;
            sig_sst_to_vort{i} = sig_sst_to_vort_i;
            sig_vort_to_sst{i} = sig_vort_to_sst_i;
        end
    end
end

save(strcat('data/F_vort_to_sst_delay_', int2str(lat), '_', int2str(lon), '_.mat'), 'F_vort_to_sst');
save(strcat('data/F_sst_to_vort_delay_', int2str(lat), '_', int2str(lon), '_.mat'), 'F_sst_to_vort');
save(strcat('data/sig_vort_to_sst_delay_', int2str(lat), '_', int2str(lon), '_.mat'), 'sig_vort_to_sst');
save(strcat('data/sig_sst_to_vort_delay_', int2str(lat), '_', int2str(lon), '_.mat'), 'sig_sst_to_vort');
