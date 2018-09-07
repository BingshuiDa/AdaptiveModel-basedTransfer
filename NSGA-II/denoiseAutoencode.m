function [offspring, n_offspring] = denoiseAutoencode(pc, generation)
%Generate offspring solutions from source problems
load allmodels_AE.mat
n_source = length(allmodels_AE);
pop = size(pc, 1);
pc = pc';
pcc = [pc; ones(1, size(pc, 2))];
offspring = [];
for i = 1:n_source
    source = allmodels_AE{i};
    pp = source.history(:, :, generation);
    pp = pp';
    ppp = [pp; ones(1, size(pp, 2))];
%     d_pp = size(pp, 2);
    bs = source.best; bs = bs';

    M = (pcc*(ppp'))/(ppp*(ppp'));
    
    M = M(1:end-1, 1:end-1);
    ls = M*bs;
    offspring = [offspring; ls'];
    
%     if d_pc <= d_pp
%         pc = [pc, zeros(size(pc, 1), d_pp - d_pc)];
%         M = (pc*(pp'))*pinv(pp*(pp'));
%         ls = M*bs;
%         offspring = [offspring; ls(:, 1:d_pc)];
%     else
%         pp = [pp, zeros(size(pp, 1), d_pc - d_pp)];
%         M = (pc*(pp'))*pinv(pp*(pp'));
%         bs = [bs, zeros(size(bs, 1), d_pc - d_pp)];
%         ls = M*bs;
%         offspring = [offspring; ls];
%     end
end

n_offspring = size(offspring, 1);
if n_offspring > pop
    k = randperm(n_offspring, pop);
    offspring = offspring(k, :);
    n_offspring = pop;
end
offspring(offspring > 1) = 1;
offspring(offspring < 0) = 0;