%Mustafa Mumtaz

%Normal Principle Coordinate Analysis
T = readtable('a1-main-data.csv','VariableNamingRule','preserve');
ids = upper(string(T{:,1}));
X   = T{:,2:end};                     % features
X(~isfinite(X)) = 0;

% Keep only requested prefixes
PFX = {'MEN','OEN','TN','LN','SN'};
keep = false(size(ids));
for k = 1:numel(PFX), keep = keep | startsWith(ids,PFX{k}); end
T   = T(keep,:); ids = ids(keep);
X   = X(keep,:);

% Normalize rows to proportions (Bray–Curtis uses composition)
rowsum = sum(X,2); nz = rowsum > 0;
P = zeros(size(X)); P(nz,:) = X(nz,:) ./ rowsum(nz);

% Bray–Curtis distance and PCoA
D = pdist(P,'cityblock')/2;
[Y, evals] = cmdscale(D, 2);

% Percent variance explained by first two positive eigenvalues
expl = [NaN NaN];
pos  = evals > 0;
if any(pos)
    tot = sum(evals(pos));
    if tot > 0
        if numel(evals) >= 1 && evals(1) > 0, expl(1) = 100*evals(1)/tot; end
        if numel(evals) >= 2 && evals(2) > 0, expl(2) = 100*evals(2)/tot; end
    end
end

% Group indices (in filtered table)
grpIdx = cell(1,numel(PFX));
for k = 1:numel(PFX), grpIdx{k} = startsWith(ids,PFX{k}); end

% --- Plot (paper style, matching your example) ---
figure('Color','w'); hold on; grid on;
set(gca,'FontSize',11,'LineWidth',0.75,'TickDir','out','Box','off');

C   = lines(numel(PFX));
Mks = {'o','^','s','d','>'};
h   = gobjects(1,numel(PFX));

for g = 1:numel(PFX)
    idx = grpIdx{g};
    if ~any(idx), continue; end

    h(g) = scatter(Y(idx,1), Y(idx,2), 42, ...
        'Marker', Mks{min(g,numel(Mks))}, ...
        'MarkerFaceColor', C(g,:), ...
        'MarkerEdgeColor', [0 0 0], ...
        'MarkerFaceAlpha', 0.75, ...
        'MarkerEdgeAlpha', 0.35, ...
        'DisplayName', sprintf('%s (n=%d)', PFX{g}, nnz(idx)));

    % 95% confidence ellipse (if enough points / full-rank covariance)
    G = Y(idx,:); G = G(all(isfinite(G),2),:);
    if size(G,1) >= 3
        mu = mean(G,1); S = cov(G);
        if all(isfinite(S(:))) && rank(S) == 2
            chi2_95 = 5.9915; k = sqrt(chi2_95);
            [V,L] = eig(S);
            t = linspace(0,2*pi,200); E = (V*sqrt(L))*[cos(t); sin(t)]*k;
            E = (E' + mu);
            patch(E(:,1), E(:,2), C(g,:), 'FaceAlpha', 0.12, ...
                  'EdgeColor', C(g,:), 'LineWidth', 1.2);
            plot(mu(1), mu(2), 'x', 'Color', C(g,:), 'LineWidth', 1.25, 'MarkerSize', 10);
        end
    end
end

% Labels/legend/title
xlab = 'PCoA 1'; if ~isnan(expl(1)), xlab = sprintf('PCoA 1 (%.1f%%)', expl(1)); end
ylab = 'PCoA 2'; if ~isnan(expl(2)), ylab = sprintf('PCoA 2 (%.1f%%)', expl(2)); end
xlabel(xlab); ylabel(ylab);
legend(h(isgraphics(h)), 'Location','eastoutside','Box','off');
title('Principal Coordinate Analysis — Selected Groups (MEN | OEN | TN | LN | SN)');

axis equal tight;
ax = axis; pad = 0.05 * max([ax(2)-ax(1), ax(4)-ax(3)]);
axis([ax(1)-pad, ax(2)+pad, ax(3)-pad, ax(4)+pad]);
hold off;

%% 
% Relative abundance — OE samples (OEN vs OE-other), using the unfiltered table
T_all  = readtable('a1-main-data.csv','VariableNamingRule','preserve');  % fresh, unfiltered
idName = T_all.Properties.VariableNames{1};

ids0 = upper(string(T_all.(idName)));
ids0 = strtrim(erase(ids0, char(160)));   % trim + remove NBSPs

M = T_all{:,2:end}; M(~isfinite(M)) = 0;

% Row-wise to %
rs = sum(M,2,'omitnan'); z = rs>0;
if any(rs > 1.0001), M(z,:) = 100*M(z,:)./rs(z);
elseif max(M,[],'all') <= 1, M = 100*M;
end

% OE universe, split OEN vs other OE*
idxOE   = startsWith(ids0,'OE');
idxOEN  = startsWith(ids0,'OEN');
idxOE_O = idxOE & ~idxOEN;

if ~any(idxOE), error('No OE* samples found.'); end

% Optional: collapse to top-K taxa for readability
featVars = T_all.Properties.VariableNames(2:end);
TOPK = 25;
mOE  = mean(M(idxOE,:),1,'omitnan');
[~,ix] = sort(mOE,'descend');
keepC  = ix(1:min(TOPK-1,numel(ix)));
other  = sum(M(:,setdiff(1:numel(featVars), keepC)),2,'omitnan');

% Avoid duplicate legend name if a column is already named "Other"
otherName = cellstr(matlab.lang.makeUniqueStrings("Other", featVars));
Mkeep = [M(:,keepC) other];
leg   = [featVars(keepC), otherName];

% Build order with a blank spacer between groups
n1   = nnz(idxOEN);
Mplt = [Mkeep(idxOEN,:); nan(1,size(Mkeep,2)); Mkeep(idxOE_O,:)];
labs = [ids0(idxOEN); " "; ids0(idxOE_O)];

figure('Color','w');
bar(Mplt,'stacked','EdgeColor','none'); colormap(lines(size(Mplt,2)));
ylabel('Relative abundance (%)'); ylim([0 100]);
set(gca,'FontSize',11,'LineWidth',0.75,'TickDir','out','Box','off', ...
        'XTick',1:numel(labs),'XTickLabel',labs,'TickLabelInterpreter','none');
xtickangle(90);
xline(n1+0.5,'k:','LineWidth',0.8);
title('Relative Abundance — OE Samples (OEN vs OE-other)');
legend(leg,'Location','eastoutside','Box','off');

%%
%Outer Ear Average Relative Abundance
