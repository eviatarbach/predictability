nobs = {13870, 2774, 924};
mhtc = 'FDR';

cd mvgc_v1.0

startup

cd ..

sst_01 = importdata('data/sst_01.mat');
vort_01 = importdata('data/vort_01.mat');
sst_05 = importdata('data/sst_05.mat');
vort_05 = importdata('data/vort_05.mat');
sst_15 = importdata('data/sst_15.mat');
vort_15 = importdata('data/vort_15.mat');

sst = {sst_01, sst_05, sst_15};
vort = {vort_01, vort_05, vort_15};

F_vort_to_sst = {NaN(480, 241), NaN(480, 241), NaN(480, 241)};
F_sst_to_vort = {NaN(480, 241), NaN(480, 241), NaN(480, 241)};
sig_90_vort_to_sst = {NaN(480, 241), NaN(480, 241), NaN(480, 241)};
sig_90_sst_to_vort = {NaN(480, 241), NaN(480, 241), NaN(480, 241)};
sig_95_vort_to_sst = {NaN(480, 241), NaN(480, 241), NaN(480, 241)};
sig_95_sst_to_vort = {NaN(480, 241), NaN(480, 241), NaN(480, 241)};
times = {NaN(480, 241), NaN(480, 241), NaN(480, 241)};

parpool(20)

parfor lat = 1:241
    for lon = 1:480
        for i = 1:3
            lat, lon
            sst = reshape(sst{i}(lon, lat, :), [1, nobs{i}]);
            vort = reshape(vort{i}(lon, lat, :), [1, nobs{i}]);

            if (sst(1) == -9999) | (vort(1) == -9999)
                continue
            end

            sst = detrend(sst);  % remove global warming signal

            X = [sst; vort];

            [~, ~, moAIC, ~] = tsdata_to_infocrit(X, 50);
	    times{i}(lon, lat) = moAIC;

            [A, SIG] = tsdata_to_var(X, moAIC);
            assert(~isbad(A), 'VAR estimation failed');
            [G, info] = var_to_autocov(A, SIG);
            var_info(info, true);

            F = autocov_to_pwcgc(G);

            F_1 = F(1, 2);  % vort -> sst
            F_2 = F(2, 1);  % sst -> vort

	    F_vort_to_sst{i}(lon, lat) = F_1;
	    F_sst_to_vort{i}(lon, lat) = F_2;

            pval = mvgc_pval(F, moAIC, nobs{i}, 1, 1, 1, 0);

	    for alpha = [0.1, 0.05]
                sig = significance(pval, alpha, mhtc);
		sig_1 = sig(1, 2);
                sig_2 = sig(2, 1);
		if alpha == 0.1
                    sig_90_vort_to_sst{i}(lon, lat) = sig_1;
                    sig_90_sst_to_vort{i}(lon, lat) = sig_2;
                else
                    sig_95_vort_to_sst{i}(lon, lat) = sig_1;
                    sig_95_sst_to_vort{i}(lon, lat) = sig_2;
		end
            end
        end
    end
end

save('data/F_vort_to_sst.mat', 'F_vort_to_sst');
save('data/F_sst_to_vort.mat', 'F_sst_to_vort');
save('data/sig_90_vort_to_sst.mat', 'sig_90_vort_to_sst');
save('data/sig_90_sst_to_vort.mat', 'sig_90_sst_to_vort');
save('data/sig_95_vort_to_sst.mat', 'sig_95_vort_to_sst');
save('data/sig_95_sst_to_vort.mat', 'sig_95_sst_to_vort');
save('data/times.mat', 'times');
