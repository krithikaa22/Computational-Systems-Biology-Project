fileName = 'Recon2M.2_BiGG_Entrez_Gene.xml';

%import model
model = readCbModel(fileName);

%initial solution before integrating expression data
init_soln = optimizeCbModel(model)

%obtaining gene-reaction association matrix : r x m where r is number of
%reactions and m is number of genes
gene_reaction_matrix = zeros(numel(model.rules), numel(model.genes));
rule_condition = zeros(numel(model.rules),1);
for i = 1:numel(model.rules)
    rule = model.rules{i};
    %using regex, parse the index of the gene associated with the reaction
    genes_in_rule = regexp(rule, '\d+', 'match');
    for j = 1:numel(genes_in_rule)
        % set the association as 1 
        gene_reaction_matrix(i, j) = 1;
    end
    %parse type of interaction between multiple genes from the rule
    if strfind(rule, '|') > 0
        rule_condition(i) = 0;
    elseif strfind(rule, '&') > 0
            rule_condition(i) = 1;
    else
            rule_condition(i) = -1;
    end
end

%HEALTHY TRANSCRIPTOME ANALYSIS

%list of genes and their constraints have been imported after analysis
%using python
id = upperboundshel.(1);
upper = upperboundshel.(2);

for i=1:1150
        % for each gene, find the reactions associated with it from the
        % gene-reaction matrix
        gene = id(i);
        col = find(strcmp(model.genes, num2str(gene)));
        rxns = find(gene_reaction_matrix(:, col) == 1);
        % set the upper bound for each reaction
        for r=1:numel(rxns)
            if model.ub(r) ~= 1000
                if rule_condition(r) == 0
                    model.ub(r) = min(upper(i), model.ub(r));
                else
                    model.ub(r) = model.ub(r) + upper(i);
                end
            end
            model.ub(r) = upper(i)
        end
end
healthy_soln = optimizeCbModel(model)

%AML TRANSCRIPTOME

%list of genes and their constraints have been imported after analysis
%using python
id = upperboundsdis.(1);
upper = upperboundsdis.(2);

for i=1:1150
        % for each gene, find the reactions associated with it from the
        % gene-reaction matrix
        gene = id(i);
        col = find(strcmp(model.genes, num2str(gene)));
        rxns = find(gene_reaction_matrix(:, col) == 1);
        % set the upper bound for each reaction
        for r=1:numel(rxns)
            if model.ub(r) ~= 1000
                if rule_condition(r) == 0
                    model.ub(r) = min(upper(i), model.ub(r));
                else
                    model.ub(r) = model.ub(r) + upper(i);
                end
            end
            model.ub(r) = upper(i)
        end
end
dis_soln = optimizeCbModel(model)

%filter upregulated and downregulated reactions using threshold
upregulated = (dis_soln.v - healthy_soln.v) > 0.0001;
downregulated = (dis_soln.v - healthy_soln.v) < -0.0001;

up_genes = [];
down_genes = [];

%obtain gene ids for the genes associated with upregulated and
%downregulated reactions
for i=1:5842
    if upregulated(i) == 1
        associated_genes = find(gene_reaction_matrix(i, :) == 1);
        for idx=1:numel(associated_genes)
            up_genes = [up_genes model.genes(idx)];
        end
    end
    if downregulated(i) == 1
        associated_genes = find(gene_reaction_matrix(i, :) == 1);
        for idx=1:numel(associated_genes)
            down_genes = [down_genes model.genes(idx)];
        end
    end
end
