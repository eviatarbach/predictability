nobs = {13870, 2774, 924};
morder = 10;
alpha = 0.1;
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

F_vort_to_sst = {zeros(480, 241), zeros(480, 241), zeros(480, 241)};
F_sst_to_vort = {zeros(480, 241), zeros(480, 241), zeros(480, 241)};
sig_90_vort_to_sst = {zeros(480, 241), zeros(480, 241), zeros(480, 241)};
sig_90_sst_to_vort = {zeros(480, 241), zeros(480, 241), zeros(480, 241)};
sig_95_vort_to_sst = {zeros(480, 241), zeros(480, 241), zeros(480, 241)};
sig_95_sst_to_vort = {zeros(480, 241), zeros(480, 241), zeros(480, 241)};

parpool(20)

parfor lat = 1:241
    for lon = 1:480
        for i = 1:3
            lat, lon
            sst = reshape(sst{i}(lon, lat, :), [1, nobs{i}]);
            vort = reshape(vort{i}(lon, lat, :), [1, nobs{i}]);
            if sst(1) == -9999
                ratios(lon, lat) = NaN;
                continue
            end

            sst = detrend(sst);  % remove global warming signal

            X = [sst; vort];
            [A, SIG] = tsdata_to_var(X, morder, 'OLS');
            assert(~isbad(A),'VAR estimation failed');
            [G,info] = var_to_autocov(A,SIG);
            var_info(info,true);

            F = autocov_to_pwcgc(G);
            pval = mvgc_pval(F,morder,nobs{i},1,1,1,0);
            sig = significance(pval,alpha,mhtc);

            F_1 = F(1, 2);  % vort -> sst
            F_2 = F(2, 1);  % sst -> vort

            sig_1 = sig(1, 2);
            sig_2 = sig(2, 1);
        end
    end
end
