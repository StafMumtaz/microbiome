% Relative abundance by SN | SDP | SD (not SDP)
sinus_data = readtable('a3-sinus.csv');

idName   = sinus_data.Properties.VariableNames{1};
genusVar = sinus_data.Properties.VariableNames(2:end);
ids      = upper(string(sinus_data.(idName)));

% Groups
idxSN  = startsWith(ids, "SN");
idxSDP = startsWith(ids, "SDP");
idxSD  = startsWith(ids, "SD") & ~idxSDP;

% Convert to percentages row-wise
M = sinus_data{:,2:end};
rowSums = sum(M,2,'omitnan');
if any(rowSums > 1.0001), M = 100*M./rowSums;
elseif max(M,[],'all') <= 1, M = 100*M;
end
sinus_data{:,2:end} = M;

% Build table with blank spacers between groups
blank = array2table(nan(1,width(sinus_data)));
blank.Properties.VariableNames = sinus_data.Properties.VariableNames;
blank{1,1} = "";

chunks = {};
if any(idxSN),  chunks{end+1} = sinus_data(idxSN ,:); end
if any(idxSN) && (any(idxSDP) || any(idxSD)), chunks{end+1} = blank; end
if any(idxSDP), chunks{end+1} = sinus_data(idxSDP,:); end
if any(idxSDP) && any(idxSD), chunks{end+1} = blank; end
if any(idxSD),  chunks{end+1} = sinus_data(idxSD ,:); end
plotTbl = vertcat(chunks{:});

% Plot stacked relative abundance
figure('Color','w');
b = bar(plotTbl{:,2:end}, 'stacked', 'EdgeColor',[1 1 1], 'LineWidth',0.25, 'BarWidth',0.98);

% Consistent colors
nGenus = numel(genusVar);
cmap = lines(nGenus);
for k = 1:nGenus, b(k).FaceColor = cmap(k,:); end

% Axes/labels
xticks(1:height(plotTbl));
xticklabels(string(plotTbl{:,1}));
xtickangle(90);
ylabel('% Relative Abundance'); ylim([0 100]);
box off; set(gca,'FontSize',10,'TickDir','out');

legend(genusVar,'Location','eastoutside','Box','off','Interpreter','none');
title('Relative Abundance by Genus — SN | SDP | SD');

%%
%Relative Abundance Aggregated

% Compute group means (percent scale already in sinus_data{:,2:end})
mean_SN  = mean(sinus_data{idxSN , 2:end}, 1, 'omitnan');
mean_SDP = mean(sinus_data{idxSDP, 2:end}, 1, 'omitnan');
mean_SD  = mean(sinus_data{idxSD , 2:end}, 1, 'omitnan');

groupNames = {'SN','SDP','SD'};
meansMat   = [mean_SN; mean_SDP; mean_SD];

% Keep only groups that exist
valid = [any(idxSN), any(idxSDP), any(idxSD)];
meansMat   = meansMat(valid, :);
groupNames = groupNames(valid);

% Table of means (optional display)
aggTbl = array2table(meansMat, 'VariableNames', genusVar, 'RowNames', groupNames);
disp('Mean relative abundance (%) by group:');
disp(aggTbl);

% Stacked bar of group means (uses same 'cmap' defined earlier)
figure('Color','w');
b2 = bar(meansMat, 'stacked', 'EdgeColor', [1 1 1], 'LineWidth', 0.25);
for k = 1:numel(genusVar)
    b2(k).FaceColor = cmap(k,:);
end
xticks(1:numel(groupNames));
xticklabels(groupNames);
ylabel('Mean % Relative Abundance');
ylim([0 100]);
box off; set(gca,'FontSize',10,'TickDir','out');
legend(genusVar, 'Location', 'eastoutside', 'Box', 'off', 'Interpreter', 'none');
title('Mean Relative Abundance by Group (SN | SDP | SD)');

%%
% Richness & Shannon: SN vs SD
T = readtable('a3-sinus.csv');

ids = upper(string(T{:,1}));
X   = T{:,2:end};                   % abundance matrix
X(~isfinite(X)) = 0;

% Groups
idxSN = startsWith(ids,'SN');
idxSD = startsWith(ids,'SD');       % includes SDP
% If you want SD but NOT SDP, use:
% idxSD = startsWith(ids,'SD') & ~startsWith(ids,'SDP');

% ----- Richness -----
rich = sum(X>0, 2);

% ----- Shannon diversity (use proportions) -----
rowsum = sum(X,2);
nz = rowsum > 0;
P = zeros(size(X));
P(nz,:) = X(nz,:) ./ rowsum(nz);
Tlog = P .* log(P); Tlog(P==0) = 0;
H = -sum(Tlog, 2); H(~nz) = NaN;

% Split groups
r1 = rich(idxSN); r2 = rich(idxSD);
h1 = H(idxSN);    h2 = H(idxSD);

% Wilcoxon rank-sum p-values
pR = NaN; if ~isempty(r1) && ~isempty(r2), pR = ranksum(r1, r2); end
pH = NaN; 
if ~isempty(h1) && ~isempty(h2)
    pH = ranksum(h1(~isnan(h1)), h2(~isnan(h2)));
end

% ----- Plot -----
figure('Color','w'); tiledlayout(1,2,'Padding','compact','TileSpacing','compact');

% Richness
nexttile;
boxplot([r1; r2], [zeros(numel(r1),1); ones(numel(r2),1)]);
set(gca,'XTickLabel',{'SN','SD'}); ylabel('OTU Richness'); title('Richness');
yl = ylim; y = yl(2)*0.95; hold on; plot([1 1 2 2],[y*0.98 y y y*0.98],'k');
txt = 'p = NA'; if ~isnan(pR), txt = sprintf('p = %.2f', pR); end
text(1.5, y, txt, 'horiz','center','vert','bottom','FontWeight','bold'); hold off;

% Shannon
nexttile;
boxplot([h1; h2], [zeros(numel(h1),1); ones(numel(h2),1)]);
set(gca,'XTickLabel',{'SN','SD'}); ylabel('Shannon Index'); title('Shannon Diversity');
yl = ylim; y = yl(2)*0.95; hold on; plot([1 1 2 2],[y*0.98 y y y*0.98],'k');
txt = 'p = NA'; if ~isnan(pH), txt = sprintf('p = %.2f', pH); end
text(1.5, y, txt, 'horiz','center','vert','bottom','FontWeight','bold'); hold off;

%%
%Repeat Diversity and Richness for different table

sinus_data_no_agg = readtable('a3b-sinus-no-aggregate.csv');

ids = upper(string(sinus_data_no_agg{:,1}));
X   = sinus_data_no_agg{:,2:end};
X(~isfinite(X)) = 0;

% Groups
idxSN = startsWith(ids,'SN');
idxSD = startsWith(ids,'SD');   % includes SDP
% If you want SD but NOT SDP, use:
% idxSD = startsWith(ids,'SD') & ~startsWith(ids,'SDP');

% ----- Richness -----
rich = sum(X > 0, 2);  % change to (X > 1e-6) if you want a small threshold

% ----- Shannon diversity (use proportions) -----
rowsum = sum(X,2);
nz = rowsum > 0;
P = zeros(size(X));
P(nz,:) = X(nz,:) ./ rowsum(nz);
Tlog = P .* log(P); Tlog(P==0) = 0;
H = -sum(Tlog, 2); H(~nz) = NaN;

% Split groups
r1 = rich(idxSN); r2 = rich(idxSD);
h1 = H(idxSN);    h2 = H(idxSD);

% Wilcoxon rank-sum p-values
pR = NaN; if ~isempty(r1) && ~isempty(r2), pR = ranksum(r1, r2); end
pH = NaN;
if ~isempty(h1) && ~isempty(h2)
    pH = ranksum(h1(~isnan(h1)), h2(~isnan(h2)));
end

% ----- Plot -----
figure('Color','w'); tiledlayout(1,2,'Padding','compact','TileSpacing','compact');

% Richness
nexttile;
boxplot([r1; r2], [zeros(numel(r1),1); ones(numel(r2),1)]);
set(gca,'XTickLabel',{'SN','SD'}); ylabel('OTU Richness'); title('Richness (no aggregation)');
yl = ylim; y = yl(2)*0.95; hold on; plot([1 1 2 2],[y*0.98 y y y*0.98],'k');
txt = 'p = NA'; if ~isnan(pR), txt = sprintf('p = %.2f', pR); end
text(1.5, y, txt, 'horiz','center','vert','bottom','FontWeight','bold'); hold off;

% Shannon
nexttile;
boxplot([h1; h2], [zeros(numel(h1),1); ones(numel(h2),1)]);
set(gca,'XTickLabel',{'SN','SD'}); ylabel('Shannon Index'); title('Shannon Diversity (no aggregation)');
yl = ylim; y = yl(2)*0.95; hold on; plot([1 1 2 2],[y*0.98 y y y*0.98],'k');
txt = 'p = NA'; if ~isnan(pH), txt = sprintf('p = %.2f', pH); end
text(1.5, y, txt, 'horiz','center','vert','bottom','FontWeight','bold'); hold off;



%% 
% PCoA (Bray–Curtis) — Sinus Samples (SN | SDP | SD-not-SDP)
T = sinus_data;                           % or: sinus_data_no_agg
ids = upper(string(T{:,1}));
X  = T{:,2:end};
X(~isfinite(X)) = 0;

% Normalize rows to proportions for Bray–Curtis
rowsum = sum(X,2);
nz = rowsum > 0;
P = zeros(size(X));
P(nz,:) = X(nz,:) ./ rowsum(nz);

% Bray–Curtis distance (0.5 * L1 on proportions) and PCoA via cmdscale
D = pdist(P,'cityblock')/2;
[Y, evals] = cmdscale(D, 2);             % Nx2 coordinates

% Percent variance explained by the first two positive eigenvalues
expl = [NaN NaN];
pos  = evals > 0;
if any(pos)
    tot = sum(evals(pos));
    if tot > 0
        if numel(evals) >= 1 && evals(1) > 0, expl(1) = 100*evals(1)/tot; end
        if numel(evals) >= 2 && evals(2) > 0, expl(2) = 100*evals(2)/tot; end
    end
end

% Groups: SN, SDP, SD (not SDP)
idxSDP = startsWith(ids,'SDP');
idxSN  = startsWith(ids,'SN');
idxSD  = startsWith(ids,'SD') & ~idxSDP;

% --- Plot (paper style) ---
figure('Color','w'); hold on; grid on;
set(gca,'FontSize',11,'LineWidth',0.75,'TickDir','out','Box','off');

C   = lines(3);                           % consistent palette
Mks = {'o','^','s'};                      % SN, SDP, SD-not-SDP
Lbl = {'SN','SDP','SD (not SDP)'};
grpIdx = {idxSN, idxSDP, idxSD};
h = gobjects(1,3);

for g = 1:3
    idx = grpIdx{g};
    if ~any(idx), continue; end
    h(g) = scatter(Y(idx,1), Y(idx,2), 42, ...
        'Marker', Mks{g}, ...
        'MarkerFaceColor', C(g,:), ...
        'MarkerEdgeColor', [0 0 0], ...
        'MarkerFaceAlpha', 0.75, ...
        'MarkerEdgeAlpha', 0.35, ...
        'DisplayName', sprintf('%s (n=%d)', Lbl{g}, nnz(idx)));

    % 95% confidence ellipse if enough points
    G = Y(idx,:); G = G(all(isfinite(G),2),:);
    if size(G,1) >= 3
        mu = mean(G,1);
        S  = cov(G);
        if all(isfinite(S(:))) && rank(S) == 2
            chi2_95 = 5.9915; k = sqrt(chi2_95);
            [V,L] = eig(S);
            t = linspace(0,2*pi,200);
            circ = [cos(t); sin(t)];
            E = (V*sqrt(L))*circ*k;    % 2x200
            E = (E' + mu);             % 200x2
            patch(E(:,1), E(:,2), C(g,:), ...
                  'FaceAlpha', 0.12, 'EdgeColor', C(g,:), 'LineWidth', 1.2);
            plot(mu(1), mu(2), 'x', 'Color', C(g,:), 'LineWidth', 1.25, 'MarkerSize', 10);
        end
    end
end

% Axis labels, legend, title
xlab = 'PCoA 1'; if ~isnan(expl(1)), xlab = sprintf('PCoA 1 (%.1f%%)', expl(1)); end
ylab = 'PCoA 2'; if ~isnan(expl(2)), ylab = sprintf('PCoA 2 (%.1f%%)', expl(2)); end
xlabel(xlab);
ylabel(ylab);

legend(h(isgraphics(h)), 'Location','eastoutside','Box','off');
title('Principal Coordinate Analysis — Sinus Samples (SN | SDP | SD-not-SDP)');

axis equal tight;
ax = axis; pad = 0.05 * max([ax(2)-ax(1), ax(4)-ax(3)]);
axis([ax(1)-pad, ax(2)+pad, ax(3)-pad, ax(4)+pad]);
hold off;
