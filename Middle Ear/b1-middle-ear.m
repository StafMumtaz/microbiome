%Mustafa Mumtaz
%Middle Ear Analysis

% Load (keep headers exactly as in the CSV)
T = readtable('a3-middle-ear.csv','VariableNamingRule','preserve');
idName = T.Properties.VariableNames{1};

% Use only numeric abundance columns (skip ID col)
numMask = varfun(@(x) isnumeric(x), T, 'OutputFormat','uniform');
numMask(1) = false;
genusLabels = T.Properties.VariableNames(numMask);              % for legend (original text)
genusVars   = matlab.lang.makeValidName(genusLabels);           % safe names for table vars

% Relative abundance to percent (robust to counts or proportions)
A  = T{:,numMask}; A(~isfinite(A)) = 0;
rs = sum(A,2,'omitnan'); z = rs>0;
if any(rs > 1.0001), A(z,:) = 100*A(z,:)./rs(z);
elseif max(A,[],'all') <= 1, A = 100*A;
end

% Split MEN vs others
ids  = string(T.(idName));
idxN = startsWith(ids,'MEN');

tblN = [table(ids(idxN),'VariableNames',{idName}), array2table(A(idxN,:),'VariableNames',genusVars)];
tblD = [table(ids(~idxN),'VariableNames',{idName}), array2table(A(~idxN,:),'VariableNames',genusVars)];

% Spacer row (same schema as first available row)
if ~isempty(tblN), template = tblN(1,:); elseif ~isempty(tblD), template = tblD(1,:); else, error('No rows to plot.'); end
blank = template; blank.(idName) = ""; blank{:,2:end} = NaN;

plotTbl = [tblN; blank; tblD];

% Plot (classic "lines" palette)
figure('Color','w');
b = bar(plotTbl{:,2:end}, 'stacked', 'EdgeColor',[1 1 1], 'LineWidth',0.25, 'BarWidth',0.98);
cmap = lines(width(plotTbl)-1); for k = 1:width(plotTbl)-1, b(k).FaceColor = cmap(k,:); end
xticks(1:height(plotTbl)); xticklabels(string(plotTbl{:,1})); xtickangle(90);
ylabel('% Relative Abundance'); ylim([0 100]); box off; set(gca,'FontSize',10,'TickDir','out');
legend(genusLabels, 'Location','eastoutside', 'Box','off', 'Interpreter','none');
title({'Relative Abundance Analysis — Middle Ear', 'MEN* (left)  |  spacer  |  Others (right)'});

%%
%Richness and diversity comparisons

A = readtable('a3b-me-no-agg.csv','VariableNamingRule','preserve');
B = readtable('b2-demographics.csv','VariableNamingRule','preserve');
leftKey = A.Properties.VariableNames{1};
k  = B.Properties.VariableNames{2};
d  = B.Properties.VariableNames(7:12);
[~, ia] = unique(string(B{:,k}),'stable');  % de-dup on the key by name (no dot-access)
B = B(ia,:);
B = B(:,[k d]);
J = outerjoin(A,B,'LeftKeys',leftKey,'RightKeys',k,'Type','left','MergeKeys',false,'RightVariables',d);
rest = setdiff(A.Properties.VariableNames,leftKey,'stable');
T2   = [J(:,leftKey) J(:,d) J(:,rest)];
T2 = T2(~isnan(T2{:,2}), :);

% ===== From T2: diversity & richness per binary factor =====
A  = table2array(T2(:,8:end)); A(~isfinite(A)) = 0;
rs = sum(A,2); nz = rs>0;
P  = zeros(size(A)); P(nz,:) = A(nz,:)./rs(nz);
rich = sum(A>0,2);
T   = P.*log(P); T(P==0) = 0;
H   = -sum(T,2); H(~nz) = NaN;

gcols  = [2 3 4 5 7];                                % use cols 2,3,4,5,7 (skip 6)
labels = {'CWD','Antibiotics','Otitis Media','Mastoidectomy','Cholesteatoma'};

for k = 1:numel(gcols)
    g = T2{:,gcols(k)}; m = ismember(g,[0 1]);
    r = rich(m); h = H(m); g = g(m);
    r0 = r(g==0); r1 = r(g==1);
    h0 = h(g==0); h1 = h(g==1);
    pR = ranksum(r0,r1);
    pH = ranksum(h0(~isnan(h0)), h1(~isnan(h1)));

    figure('Color','w'); tiledlayout(1,2,'Padding','compact','TileSpacing','compact');

    % Richness
    nexttile;
    boxplot([r0; r1],[zeros(numel(r0),1); ones(numel(r1),1)]);
    set(gca,'XTickLabel',{'No','Yes'}); ylabel('OTU Richness'); title('Richness');
    yl = ylim; y = yl(2)*0.95; hold on; plot([1 1 2 2],[y*0.98 y y y*0.98],'k');
    text(1.5,y,sprintf('p = %.2f',pR),'horiz','center','vert','bottom','FontWeight','bold'); hold off;

    % Shannon
    nexttile;
    boxplot([h0; h1],[zeros(numel(h0),1); ones(numel(h1),1)]);
    set(gca,'XTickLabel',{'No','Yes'}); ylabel('Shannon Index'); title('Shannon Diversity');
    yl = ylim; y = yl(2)*0.95; hold on; plot([1 1 2 2],[y*0.98 y y y*0.98],'k');
    text(1.5,y,sprintf('p = %.2f',pH),'horiz','center','vert','bottom','FontWeight','bold'); hold off;

    sgtitle(strrep(labels{k},'_','\_'));
end

%% ===== One-sided combine + severity trend (use cols 2,3,4,5; skip 6; include 7 later) =====
labels4 = {'CWD','Antibiotics','Otitis Media','Mastoidectomy'};
g4      = [2 3 4 5];

p1   = nan(1,numel(g4)); zval = nan(1,numel(g4)); dmed = nan(1,numel(g4));
for i = 1:numel(g4)
    g = T2{:,g4(i)}; m = ismember(g,[0 1]) & isfinite(H);
    h0 = H(m & g==0); h1 = H(m & g==1);
    if isempty(h0) || isempty(h1), continue; end
    dmed(i) = median(h1) - median(h0);                                  % Yes - No
    [p1(i),~,st] = ranksum(h1,h0,'tail','left','method','approximate'); % test Yes < No
    zval(i) = st.zval;
end

good = isfinite(p1) & p1>0;
X = -2*sum(log(p1(good)));
p_fisher = 1 - chi2cdf(X, 2*nnz(good));

rng('default'); B = 10000; Zobs = -nansum(zval);
Zsim = nan(B,1); n = numel(H);
for b = 1:B
    Hp = H(randperm(n)); zb = 0;
    for i = 1:numel(g4)
        g = T2{:,g4(i)}; m = ismember(g,[0 1]) & isfinite(Hp);
        x = Hp(m & g==1); y = Hp(m & g==0);
        if isempty(x) || isempty(y), continue; end
        [~,~,st] = ranksum(x,y,'tail','left','method','approximate');
        zb = zb - st.zval;
    end
    Zsim(b) = zb;
end
p_perm = mean(Zsim >= Zobs);

% ----- Severity from cols [2 3 4 5 7] (NaN -> 0) -----
sevcols = [2 3 4 5 7];
G = T2{:,sevcols}; M = (G==0)|(G==1)|isnan(G); keep = all(M,2) & isfinite(H);
G = G(keep,:); Hk = H(keep);
G(~isfinite(G)) = 0; sev = sum(G,2);

% Spearman trend
[rho_s, p_spear] = corr(sev, Hk, 'type','Spearman','rows','complete');

% ===== VISUALS =====
figure('Color','w');

% A) Per-factor effect summary
subplot(1,2,1);
K = numel(g4);
plot(dmed, 1:K, 'o','MarkerFaceColor',[0.2 0.2 0.8],'MarkerEdgeColor','k'); hold on;
xline(0,'k:'); set(gca,'YTick',1:K,'YTickLabel',labels4);
xlabel('Median difference: Yes - No (Shannon)'); title('Per-factor effect (negative = lower)');
box off; grid on; grid minor;
xlims = xlim;
text(xlims(1)+0.02*range(xlims), K+0.6, { ...
    sprintf('Fisher (1-sided) p = %.3g', p_fisher), ...
    sprintf('Permutation combined p = %.3g', p_perm)}, 'VerticalAlign','bottom');

% B) Severity regression with 95%% CI (jittered points)
subplot(1,2,2); hold on;
j = 0.15; sev_j = sev + (rand(size(sev))-0.5)*2*j;
scatter(sev_j, Hk, 18, 'filled', 'MarkerFaceAlpha', 0.6);
mdl = fitlm(sev, Hk);
xg  = linspace(min(sev), max(sev), 200)'; [yg, yci] = predict(mdl, xg);
fill([xg; flipud(xg)], [yci(:,1); flipud(yci(:,2))], [0.85 0.9 1], 'EdgeColor','none','FaceAlpha',0.5);
plot(xg, yg, 'k-', 'LineWidth', 1.5);
xlabel('Severity (sum of cols 2,3,4,5,7)'); ylabel('Shannon Index');
title(sprintf('Trend: Spearman \\rho = %.2f, p = %.3g', rho_s, p_spear));
box off; grid on; grid minor; hold off;

figure('Color','w');
K   = numel(labels4);
pad = 0.8;                                 % increase to pull closer
y   = linspace(1+pad, K-pad, K);           % squeezed y-positions
plot(dmed, y, 'o','MarkerFaceColor',[0.2 0.2 0.8],'MarkerEdgeColor','k'); hold on;
xline(0,'k:');
set(gca,'YTick',y,'YTickLabel',labels4);
ylim([1 K]);                               % keep margins above/below
xlabel('Median difference: Yes - No (Shannon)');
title('Per-factor effect (negative = lower)'); box off; grid on; grid minor;

%%
%Spearman for richness
% ----- Severity (cols 2,3,4,5,7) vs Richness -----
A = table2array(T2(:,8:end)); A(~isfinite(A)) = 0;
rich = sum(A>0,2);                       % use >0; change to >1e-6 if you prefer

sevcols = [2 3 4 5 7];
G = T2{:,sevcols};
valid = (G==0)|(G==1)|isnan(G);
keep = all(valid,2) & isfinite(rich);
G = G(keep,:); r = rich(keep);
G(~isfinite(G)) = 0;                     % NaN -> 0
sev = sum(G,2);

% Spearman
[rho_s, p_spear] = corr(sev, r, 'type','Spearman', 'rows','complete');

% Plot with jitter + linear fit and 95% CI
figure('Color','w'); hold on;
j = 0.15; sev_j = sev + (rand(size(sev))-0.5)*2*j;
scatter(sev_j, r, 18, 'filled', 'MarkerFaceAlpha', 0.6);
mdl = fitlm(sev, r);
xg = linspace(min(sev), max(sev), 200)'; [yg, yci] = predict(mdl, xg);
fill([xg; flipud(xg)], [yci(:,1); flipud(yci(:,2))], [0.85 0.9 1], ...
     'EdgeColor','none','FaceAlpha',0.5);
plot(xg, yg, 'k-', 'LineWidth', 1.5);
xlabel('Severity (sum of cols 2,3,4,5,7)'); ylabel('OTU Richness');
title(sprintf('Trend: Spearman \\rho = %.2f, p = %.3g', rho_s, p_spear));
box off; grid on; grid minor; hold off;

%%

% --- Load & normalize ---
T = readtable('a3-middle-ear.csv','VariableNamingRule','preserve');
M = T{:,2:end};                        % numeric abundances only
M(~isfinite(M)) = 0;
rs = sum(M,2,'omitnan'); nz = rs>0;
if any(rs > 1.0001)                    % counts -> proportions
    M(nz,:) = M(nz,:)./rs(nz);
elseif max(M,[],'all') > 1.0001        % percents -> proportions
    M = M./100;
end

% --- A priori contribution: sum of table columns 2,5,9 ---
sel_tbl = [2 5 9];                     % table column numbers (1st col = IDs)
y = sum(M(:, sel_tbl-1), 2, 'omitnan');  % -1 because M starts at table col 2
ids = upper(string(T{:,1}));

% --- Split MEN (Normal) vs others (Diseased) ---
isMEN = startsWith(ids,'MEN');
xN = y(isMEN);   xN = xN(isfinite(xN)); % MEN (Normal)
xD = y(~isMEN);  xD = xD(isfinite(xD)); % MED (Diseased)

% --- Summary stats & tests ---
nN = numel(xN);  nD = numel(xD);
mN = mean(xN);   mD = mean(xD);
sN = std(xN);    sD = std(xD);
pMW = ranksum(xN, xD);                 % Mann-Whitney (two-sided)
[~, pWT] = ttest2(xN, xD, 'Vartype','unequal');  % Welch t-test

% --- Plot ---
figure('Color','w'); hold on;
boxplot([xN; xD], [zeros(nN,1); ones(nD,1)], ...
        'Labels', {'MEN (Normal)','MED (Diseased)'});

% overlay jittered points
j = 0.12;
scatter(1 + (rand(nN,1)-0.5)*2*j, xN, 18, 'filled', 'MarkerFaceAlpha', 0.7);
scatter(2 + (rand(nD,1)-0.5)*2*j, xD, 18, 'filled', 'MarkerFaceAlpha', 0.7);

% mark means
plot(1, mN, 'kd','MarkerFaceColor','k');
plot(2, mD, 'kd','MarkerFaceColor','k');

ylabel('Contribution from A. Priori Bacterium');
ylim([0 1]); xlim([0.5 2.5]); box off; set(gca,'TickDir','out','FontSize',11);

title({ ...
  'Differences in Relative Contribution of Corynebacterium, Cutibacterium, and Pseudomonas to Bacterial Makeup in Diseased vs Normal States', ...
  sprintf('MEN (Normal): n=%d, mean=%.3f, SD=%.3f   |   MED (Diseased): n=%d, mean=%.3f, SD=%.3f', nN, mN, sN, nD, mD, sD), ...
  sprintf('Mann-Whitney p=%.4g   |   Welch t-test p=%.4g', pMW, pWT) ...
});
