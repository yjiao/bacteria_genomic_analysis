%% Comments

% This script was used to analyse sequencing data for Jiao, Y and Baym M.
% et al, "Population diversity reverses the outcome of antibiotic cycling".
% Manuscrpt in progress.

% fastq files from Illumina were first processed using breseq
% (https://github.com/barricklab/breseq).

% This script provides further filtering and allows visualization of breseq
% parameters for quality control.

%% Set paths
cd ~/Dropbox/Kishony_Lab/2015seqdata/breseq/
folder = 'fullRun20151221/postProcess/04merged';
fcfolder = 'fullRun20151221/postProcess/06forceCalled';
unfolder = 'fullRun20151221/postProcess/07unknowns';

%% import 0/1 matrix of whether a strain has a certain mutation
mutmat = readtable([folder '/matrix.txt'], 'delimiter', '\t');
strains = mutmat.Properties.VariableNames;
fcmat = readtable([fcfolder '/matrix.txt'], 'delimiter', '\t', 'ReadVariableNames', 0);
unmat = readtable([unfolder '/matrix.txt'], 'delimiter', '\t', 'ReadVariableNames', 0);

%% import mutation information
mutations = readtable([folder '/mutations.txt'], 'delimiter', '\t');

%% Import Parameters

% Import variant frequency information (new in cycling)
freq = LoadBreseqParam([folder '/minmax_params/variant_frequency']);
freqFC = LoadBreseqParam([fcfolder '/minmax_params/variant_frequency']);

%%%%% consensus information
% breseq determines the base with the highest value of L, and records read alignment evidence if this base is different
% from the reference base. This evidence is assigned log10 L minus the log10 of the cumulative length of all reference
% sequences as a consensus E-value score for this consensus base prediction.
consensus = LoadBreseqParam([folder '/minmax_params/' '/consensus_score']);
consensusFC = LoadBreseqParam([fcfolder '/minmax_params/' '/consensus_score']);

%%%%% polymorphism information
% tests the hypothesis that reads aligned to each reference position (and base insertions relative to the reference)
% support a model that is a mixture of a major and minor mutational variant as opposed to a model that all disagreements
% % with the reference sequence (or consensus change predicted as above) are due to sequencing errors.

% breseq then tests the statistical support for the model having only one reference base in the sequenced sample versus
% the model with one additional free parameter consisting of mixture of two alleles using a likelihood-ratio test. That
% is, twice the natural logarithm of the probability of the mixture model divided by the probability of the one-base
% model is compared to a chi-squared distribution with 1 degree of freedom.
% 
% As for the case of consensus mutation prediction from read alignment evidence, the p-value significance of the
% likelihood-ratio test is finally converted to a polymorphism E-value score by multiplying by the total number of
% reference positions.
polymorphism = LoadBreseqParam([folder '/minmax_params/' '/polymorphism_score']);
polymorphismFC = LoadBreseqParam([fcfolder '/minmax_params/' '/polymorphism_score']);

%%%%% strand bias
% Fisher Exact test
strandBias = LoadBreseqParam([folder '/minmax_params/' '/fisher_strand_p_value']);
strandBiasFC = LoadBreseqParam([fcfolder '/minmax_params/' '/fisher_strand_p_value']);

%%%%% quality bias
% KS test on quality values for ref vs alt
ksTest = LoadBreseqParam([folder '/minmax_params/' '/ks_quality_p_value']);
ksTestFC = LoadBreseqParam([fcfolder '/minmax_params/' '/ks_quality_p_value']);

%%%%% bias_e_value
% not sure what this is...
biasE = LoadBreseqParam([folder '/minmax_params/' '/bias_e_value']);
biasEFC = LoadBreseqParam([fcfolder '/minmax_params/' '/bias_e_value']);

%%%%% bias_p_value
% not sure what this is...
biasP = LoadBreseqParam([folder '/minmax_params/' '/bias_p_value']);
biasPFC = LoadBreseqParam([fcfolder '/minmax_params/' '/bias_p_value']);

%%%%% major allele frequency
majorAFreq = LoadBreseqParam([folder '/minmax_params/' '/major_frequency']);
majorAFreqFC = LoadBreseqParam([fcfolder '/minmax_params/' '/major_frequency']);

%%%%% total coverage
covTotal = LoadBreseqParam([folder '/minmax_params/' '/total_cov_total']);
covTotalFC = LoadBreseqParam([fcfolder '/minmax_params/' '/total_cov_total']);

%%%%% top and bottom coverage
covTop = LoadBreseqParam([folder '/minmax_params/' '/total_cov_top']);
covBot = LoadBreseqParam([folder '/minmax_params/' '/total_cov_bottom']);
covTopFC = LoadBreseqParam([fcfolder '/minmax_params/' '/total_cov_top']);
covBotFC = LoadBreseqParam([fcfolder '/minmax_params/' '/total_cov_bottom']);

%%%%% html links
html = readtable([folder '/evidencehtml.txt'], 'delimiter', '\t', 'ReadVariableNames', 1, 'TreatAsEmpty', 'none');

%%%%% Total # of reads
% need a lot of manuvering in order to get this to columns = strains in the
% same order as all other matrices. Initially this is read in where columns
% are different data types.
summary = readtable('fullRun20151221/postProcess/summary.txt', 'delimiter', '\t', 'ReadVariableNames', 0);
summary.Properties.VariableNames = {'Strain','Reads','Percent_Alignment'};

summaryStrains = summary.Strain;
reorderidx = zeros(length(summaryStrains),1);

% add an 'x' to the start of strain IDs that start with numbers, since
% matlab did this automatically earlier.
for i = 1:length(strains)
    str=summaryStrains{i};
    if isstrprop(str(1),'digit')
        summaryStrains{i} = ['x' str];
    end
end

for i = 1:length(strains)
    str = strains{i};
    match = find(strcmp(summaryStrains, str));
    reorderidx(i) =  match;
end

summaryStrains = summaryStrains(reorderidx);
summaryReads = summary.Reads(reorderidx);
summaryAlign = summary.Percent_Alignment(reorderidx);

% this section just double checks that the reindexing was done correctly
for i = 1:length(strains)
    if ~strcmp(summaryStrains{i}, strains{i})
        {summaryStrains{i}, strains{i}}
    end
end

summary = array2table(summaryReads');
summary = [summary; array2table(summaryAlign')];
summary.Properties.RowNames = {'Reads','Percent_Alignment'};
summary.Properties.VariableNames = summaryStrains;

%% Filter out strains with too few reads or too low alignment, wt calls
min_reads = 100000;
filter_reads = summaryReads < min_reads;

min_align = 90;
filter_align = summaryAlign < min_align;

filter_strains = filter_reads | filter_align;

% filter out wildtype
wt_flag = mutmat.WT==1;
% mutmat{wt_flag, :} = 0;
'strain filter done'
%% Use optimal row/col filter based on previous section
clc

orthresh = 128;
% fcthresh = 16; % best 6
% unthresh = 19; % best 19
% matchthresh = 13; % best 13

ormat = fcmat{:,:} > 0 | unmat{:,:} > 0 | mutmat{:,:} > 0;
filter_probablyAncestral =  sum(ormat, 2) >= orthresh;

filter_row = ~wt_flag & ~filter_probablyAncestral;
filter_col = ~filter_strains;
'row filter done'
%% Export mutation filter
mutation_filter = mutations(filter_row,:);
save('mutation_filter_20160526.mat', 'mutation_filter');


%% apply row/col filters to everything
mutmat = mutmat(filter_row, filter_col);
mutmat = mutmat{:,:};
strains = strains(filter_col);
fcmat = fcmat(filter_row, filter_col);
fcmat = fcmat{:,:};
unmat = unmat(filter_row, filter_col);
unmat = unmat{:,:};
mutations = mutations(filter_row,:);

html = html(filter_row, filter_col);

param_freq = freq.maxMat(filter_row, filter_col);
param_freq(isnan(param_freq))=0;
param_freqFC = freqFC.maxMat(filter_row, filter_col);
param_freqFC(isnan(param_freqFC))=0;

param_consensus = consensus.maxMat(filter_row, filter_col);
param_consensusFC = consensusFC.maxMat(filter_row, filter_col);
param_polymorphism = polymorphism.maxMat(filter_row, filter_col);
param_polymorphismFC = polymorphismFC.maxMat(filter_row, filter_col);
param_strandBias = strandBias.maxMat(filter_row, filter_col);
param_strandBiasFC = strandBiasFC.maxMat(filter_row, filter_col);
param_ksTest = ksTest.maxMat(filter_row, filter_col);
param_ksTestFC = ksTestFC.maxMat(filter_row, filter_col);
param_biasE = biasE.maxMat(filter_row, filter_col);
param_biasEFC = biasEFC.maxMat(filter_row, filter_col);
param_biasP = biasP.maxMat(filter_row, filter_col);
param_biasPFC = biasPFC.maxMat(filter_row, filter_col);
param_majorAFreq = majorAFreq.maxMat(filter_row, filter_col);
param_majorAFreqFC = majorAFreqFC.maxMat(filter_row, filter_col);
param_covTotal = covTotal.maxMat(filter_row, filter_col);
param_covTotalFC = covTotalFC.maxMat(filter_row, filter_col);
param_covTop = covTop.maxMat(filter_row, filter_col);
param_covBot = covBot.maxMat(filter_row, filter_col);
param_covTopFC = covTopFC.maxMat(filter_row, filter_col);
param_covBotFC = covBotFC.maxMat(filter_row, filter_col);

summary = summary(:, filter_col);
summaryStrains = summaryStrains(filter_col);
summaryReads = summaryReads(filter_col);
summaryAlign = summaryAlign(filter_col);
'apply filter done'
%% create strain filter for cycling
filter_notcycling = cellfun(@isempty, strfind(strains, 'XA'));
filter_notcycling = filter_notcycling & cellfun(@isempty, strfind(strains, 'XB'));
filter_notcycling = filter_notcycling & cellfun(@isempty, strfind(strains, 'XC'));
filter_notcycling = filter_notcycling & cellfun(@isempty, strfind(strains, 'CXB'));
filter_notcycling = filter_notcycling & cellfun(@isempty, strfind(strains, 'CXC'));
filter_notcycling = filter_notcycling & cellfun(@isempty, strfind(strains, 'CXA'));

filter_notcycling = filter_notcycling & cellfun(@isempty, strfind(strains, 'NA'));
filter_notcycling = filter_notcycling & cellfun(@isempty, strfind(strains, 'NB'));
filter_notcycling = filter_notcycling & cellfun(@isempty, strfind(strains, 'NC'));
filter_notcycling = filter_notcycling & cellfun(@isempty, strfind(strains, 'CNB'));
filter_notcycling = filter_notcycling & cellfun(@isempty, strfind(strains, 'CNC'));
filter_notcycling = filter_notcycling & cellfun(@isempty, strfind(strains, 'CNA'));
filter_notcycling = filter_notcycling & cellfun(@isempty, strfind(strains, 'CN_A'));
filter_notcycling = filter_notcycling & cellfun(@isempty, strfind(strains, 'XN_A'));
filter_notcycling = filter_notcycling & cellfun(@isempty, strfind(strains, 'XNB'));
filter_notcycling = filter_notcycling & cellfun(@isempty, strfind(strains, 'XNC'));

filter_notcycling = filter_notcycling & cellfun(@isempty, strfind(strains, 'HL'));

filter_step1cycling = [...
    find(strncmpi(strains, 'XA1',3))...
    find(strncmpi(strains, 'XB1',3))...
    find(strncmpi(strains, 'XC1',3))...
    find(strncmpi(strains, 'NA1',3))...
    find(strncmpi(strains, 'NB1',3))...
    find(strncmpi(strains, 'NC1',3))...
    ];

% these are old:
filter_replacedStrains = [...
    find(strncmpi(strains, 'XA_1',5))...
    find(strncmpi(strains, 'XB_1',5))...
    find(strncmpi(strains, 'XC_1',5))...
    find(strncmpi(strains, 'NA_1',5))...
    find(strncmpi(strains, 'NB_1',5))...
    find(strncmpi(strains, 'NC_1',5))...
    ];


step2 = cellfun(@isempty,regexp(strains, 'XE[69]_NOV*'))';
% step2: 1 = NOT STEP 2
% filter_notcycling: 1 = not cycling sample
% filter_step1cycling: indices where we have first step of cycling
% filter_replacedStrains: indices where strains have been resequenced

% 1 step only
onlystep1 = step2 & filter_notcycling'; 
onlystep1(filter_step1cycling) = 1;

% 2 or more steps
onlystep2 = ~step2 | ~filter_notcycling'; 
onlystep2(filter_replacedStrains) = 0;

%% Create matrix filters
clc
strandthresh = 7;
filter_strandCoverage = param_covTop < strandthresh | param_covBot < strandthresh;
mutmat(filter_strandCoverage) = 0;
'matrix filter done'
%% manual filters: misalignment hotspots
orig = mutmat;

targets = [292584 294027 2723571 289820 289773 289809 289784 293026 2295869 550312];
buffer = 100;
pos = mutations.position;

for i = 1:length(targets)
    target = targets(i);
    flag = (pos >= target-buffer) & (pos <= target+buffer);
    nflag = sum(orig(flag,:),1) > 0;
    
    mutmat(flag, :) = zeros(sum(flag), size(mutmat, 2));
    
    disp('-----')
    sum(flag)
    ncalls1  = sum(orig, 1);
    before = ncalls1(nflag)
    ncalls2  = sum(mutmat, 1);
    after = ncalls2(nflag)

%     origcalls = sum(orig, 1);
%     origcalls(flag)
%     strains(flag)
end

% filter out random polymorphic calls??
mutmat(211, 12) = 0;

% test = mutmat{:,:};
% test(filter_strandCoverage) = 0;
% sum(sum(test(:,onlystep2),1)==0)
% sum(sum(test(:,onlystep2),1)==1)
% sum(sum(test(:,onlystep2),1)==2)
% sum(sum(test(:,onlystep2),1)>=3)
'manual filter done'
%% Automatic filter of misalignment hotspots
% for each column in mut mat, look at diff of positions

for i = 1:size(mutmat,2)
    mutCol = mutmat(:,i);
    pos = mutations.position;
    buffer = 100;
    flag = diff(pos(mutCol==1)) < buffer;
    if sum(flag) > 0
        i
        strains{i}
        sum(mutCol)
        idx = find(flag);
        idx = unique([idx; idx + 1]);
        hasmutpos = find(mutCol);
        filterpos = hasmutpos(idx);
        mutmat(filterpos,:) = 0;
        sum(mutmat(:,i))
    end

end
'hotspot filter done'

%% Get rid of mutations that are now empty
filter_emptyRows = ~sum(mutmat,2)==0;
filter_emptyCols = ~sum(mutmat,1)==0;

mutmat = mutmat(filter_emptyRows, filter_emptyCols);
strains = strains(filter_emptyCols);
fcmat = fcmat(filter_emptyRows, filter_emptyCols);
unmat = unmat(filter_emptyRows, filter_emptyCols);
mutations = mutations(filter_emptyRows,:);

html = html(filter_emptyRows, filter_emptyCols);

param_consensus = param_consensus(filter_emptyRows, filter_emptyCols);
param_consensusFC = param_consensusFC(filter_emptyRows, filter_emptyCols);
param_polymorphism = param_polymorphism(filter_emptyRows, filter_emptyCols);
param_polymorphismFC = param_polymorphismFC(filter_emptyRows, filter_emptyCols);
param_strandBias = param_strandBias(filter_emptyRows, filter_emptyCols);
param_strandBiasFC = param_strandBiasFC(filter_emptyRows, filter_emptyCols);
param_ksTest = param_ksTest(filter_emptyRows, filter_emptyCols);
param_ksTestFC = param_ksTestFC(filter_emptyRows, filter_emptyCols);
param_biasE = param_biasE(filter_emptyRows, filter_emptyCols);
param_biasEFC = param_biasEFC(filter_emptyRows, filter_emptyCols);
param_biasP = param_biasP(filter_emptyRows, filter_emptyCols);
param_biasPFC = param_biasPFC(filter_emptyRows, filter_emptyCols);
param_majorAFreq = param_majorAFreq(filter_emptyRows, filter_emptyCols);
param_majorAFreqFC = param_majorAFreqFC(filter_emptyRows, filter_emptyCols);
param_covTotal = param_covTotal(filter_emptyRows, filter_emptyCols);
param_covTotalFC = param_covTotalFC(filter_emptyRows, filter_emptyCols);
param_covTop = param_covTop(filter_emptyRows, filter_emptyCols);
param_covBot = param_covBot(filter_emptyRows, filter_emptyCols);
param_covTopFC = param_covTopFC(filter_emptyRows, filter_emptyCols);
param_covBotFC = param_covBotFC(filter_emptyRows, filter_emptyCols);

summary = summary(:, filter_emptyCols);
summaryStrains = summaryStrains(filter_emptyCols);
summaryReads = summaryReads(filter_emptyCols);
summaryAlign = summaryAlign(filter_emptyCols);

%% Filter out polymorphic sites in HTML files (matrix filter)
% do this after "get rid of mutations" since otherwise it will take forever
% url = 'file:///Users/Joy/Dropbox/Kishony_Lab/2015seqdata/breseq/fullRun20151221/CNC5/output/evidence/SNP_3.html';
polymflag = zeros(size(html));
for i = 1:size(html, 1)
    for j = 1:size(html,2)
        
        url = html{i,j};
        url = ['file:///Users/Joy/Dropbox/Kishony_Lab/2015seqdata/breseq/' url{1}];
        if ~isempty(url) && mutmat(i,j) == 1
            htmltxt = urlread(url);
            pattern = '<tr class="polymorphism_table_row">';
            polymflag(i,j) = ~isempty(strfind(htmltxt, pattern));
        end
    end
end

mutmat(polymflag==1) = 0;
'HTML filter done'

%% Get rid of mutations that are now empty--again
filter_emptyRows = ~sum(mutmat,2)==0;
filter_emptyRows = filter_emptyRows;

filter_emptyCols = ~sum(mutmat,1)==0;

mutmat = mutmat(filter_emptyRows, filter_emptyCols);
strains = strains(filter_emptyCols);
fcmat = fcmat(filter_emptyRows, filter_emptyCols);
unmat = unmat(filter_emptyRows, filter_emptyCols);
mutations = mutations(filter_emptyRows,:);

html = html(filter_emptyRows, filter_emptyCols);

param_freq = param_freq(filter_emptyRows, filter_emptyCols);
param_freqFC = param_freqFC(filter_emptyRows, filter_emptyCols);
param_consensus = param_consensus(filter_emptyRows, filter_emptyCols);
param_consensusFC = param_consensusFC(filter_emptyRows, filter_emptyCols);
param_polymorphism = param_polymorphism(filter_emptyRows, filter_emptyCols);
param_polymorphismFC = param_polymorphismFC(filter_emptyRows, filter_emptyCols);
param_strandBias = param_strandBias(filter_emptyRows, filter_emptyCols);
param_strandBiasFC = param_strandBiasFC(filter_emptyRows, filter_emptyCols);
param_ksTest = param_ksTest(filter_emptyRows, filter_emptyCols);
param_ksTestFC = param_ksTestFC(filter_emptyRows, filter_emptyCols);
param_biasE = param_biasE(filter_emptyRows, filter_emptyCols);
param_biasEFC = param_biasEFC(filter_emptyRows, filter_emptyCols);
param_biasP = param_biasP(filter_emptyRows, filter_emptyCols);
param_biasPFC = param_biasPFC(filter_emptyRows, filter_emptyCols);
param_majorAFreq = param_majorAFreq(filter_emptyRows, filter_emptyCols);
param_majorAFreqFC = param_majorAFreqFC(filter_emptyRows, filter_emptyCols);
param_covTotal = param_covTotal(filter_emptyRows, filter_emptyCols);
param_covTotalFC = param_covTotalFC(filter_emptyRows, filter_emptyCols);
param_covTop = param_covTop(filter_emptyRows, filter_emptyCols);
param_covBot = param_covBot(filter_emptyRows, filter_emptyCols);
param_covTopFC = param_covTopFC(filter_emptyRows, filter_emptyCols);
param_covBotFC = param_covBotFC(filter_emptyRows, filter_emptyCols);

summary = summary(:, filter_emptyCols);
summaryStrains = summaryStrains(filter_emptyCols);
summaryReads = summaryReads(filter_emptyCols);
summaryAlign = summaryAlign(filter_emptyCols);
'get rid of empties done'
%% import KEGG information
kegg = readtable('kegg/sao_parsed.txt', 'delimiter', '\t', 'ReadVariableNames', 0);
kegg.Properties.VariableNames = {'LocusTag','Gene','Product','Description','KO','C','B','A'};

% create dictionary of locusTag to gene name
dict_Locus2Gene = containers.Map(kegg.LocusTag, kegg.Gene);
dict_Locus2Product = containers.Map(kegg.LocusTag, kegg.Product);
dict_Locus2Description = containers.Map(kegg.LocusTag, kegg.Description);
dict_Locus2KO = containers.Map(kegg.LocusTag, kegg.KO);
dict_Locus2A = containers.Map(kegg.LocusTag, kegg.A);
dict_Locus2B = containers.Map(kegg.LocusTag, kegg.B);
dict_Locus2C = containers.Map(kegg.LocusTag, kegg.C);

keggStruct = struct;
keggStruct.dict_Locus2Gene = dict_Locus2Gene;
keggStruct.dict_Locus2Product = dict_Locus2Product;
keggStruct.dict_Locus2Description = dict_Locus2Description;
keggStruct.dict_Locus2KO = dict_Locus2KO;
keggStruct.dict_Locus2A = dict_Locus2A;
keggStruct.dict_Locus2B = dict_Locus2B;
keggStruct.dict_Locus2C = dict_Locus2C;
'kegg done'
%% Create Mutation objects
MUTATIONS=Mutation.empty(size(mutations,1), 0);
for i = 1:size(mutations,1)
    MUTATIONS(i) = Mutation(mutations(i,:), keggStruct);
end
'mutation obj done'
%% Create Strain objects
clc
STRAINS = Strain.empty(length(strains),0); % this creates an empty array of class Strain

% preload databases to save some time
database_oxa = loadLibrary('OXA', 3);
database_nov = loadLibrary('NOV', 3);

% preload cycling mics
filenames = dir('cyclingData/*.mat');
cyclingMICs = struct();
for i = 1:length(filenames)
    file = filenames(i).name;
    oxaphen = load(['cyclingData/', file]);
    field = fieldnames(oxaphen);
    if strcmp(field{1}, 'data') == 1
        culture = strrep(file, '.mat', '');
        culture = strrep(culture, 'Cycling', '');
    else
        culture = upper(strrep(field{1}, 'data_', ''));
    end
    cyclingMICs.(culture) = oxaphen.(field{1});
end

for i = 1:length(strains)
    str = strains{i};
    % fix matlab auto formatting: this gets rid of the number in front of a
    % bunch of "XPos" strains. Ex. "0XPos" -> 'XPos'
    if str(1) == 'x'
        str = str(strfind(str,'XPos'):end);
    end
    % fix weird mistake writing 'XN' instead of 'CN' for control cycling
    % populations
    if strcmp(str(1:2), 'XN') && ~strcmp(str(1:4), 'Xneg') && ~strcmp(str(1:4), 'XNeg')
        str(1:2) = 'CX';
    end
    
    if ~isempty(strfind(str, 'XE6')) || ~isempty(strfind(str, 'XE9')) || ~isempty(strfind(str, 'C9_'))
    %%%%%%%% use E2 parser
        if  ~isempty(strfind(str, 'C9_'))
            pattern='^(?<drug>NOV|OXA)(?<conc>\d)(?<culture>C9)_(?<phenotype>\S+)_(?<set>e2)_(?<plate>UL\d*|UR\d*|LL\d*|LR\d*)_(?<well>\w\d{1,2})';
        elseif ~isempty(strfind(str, 'TSB'))
            pattern='^(?<culture>XE6|XE9)_(?<drug>TSB)_(?<phenotype>\S+)_(?<set>xn1)_(?<plate>UL\d*|UR\d*|LL\d*|LR\d*)_(?<well>\w\d{1,2})';
        else
            pattern='^(?<culture>XE6|XE9)_(?<drug>NOV|OXA)(?<conc>\d)_(?<phenotype>\S+)_(?<set>xn2)_(?<plate>UL\d*|UR\d*|LL\d*|LR\d*)_(?<well>\w\d{1,2})';
        end
        tokens=regexp(str,pattern,'names');
        if isempty(tokens)
            disp('not found in E2:')
            disp(str)
        end
        mic = subsetE2(tokens); %[oxa, nov]
        
        % fill in the strainStruct structure for passing to the Strain
        % Class constructor
        strainStruct = struct();
        strainStruct.name=str;
        if ~isempty(strfind(str, 'XE6')) || ~isempty(strfind(str, 'XE9'))
            strainStruct.library = 'OXA';
        else
            strainStruct.library = tokens.drug;
        end
        
        strainStruct.quadrant = tokens.plate;
        strainStruct.wellname96 = tokens.well;
        strainStruct.phenotype = upper(tokens.phenotype); % convert to all upper case
        
        if strcmp(tokens.drug, 'TSB') || strcmp(tokens.culture, 'C9')
           strainStruct.conc = -1;
           strainStruct.step = 1;
        else
            strainStruct.conc = str2double(tokens.conc);
            strainStruct.step = 2;
        end
        
        
        if strcmp(tokens.drug, 'TSB')
            tokens.drug='OXA';
        end
        strainStruct.experiment = 'E2';
        strainStruct.culture = tokens.culture;
        strainStruct.platename = tokens.set;
        strainStruct.MIC = mic;
        strainStruct.mic_OXA = mic(1);
        strainStruct.mic_NOV = mic(2);
        
        STRAINS(i) = Strain(strainStruct);
    
    elseif ~isempty(strfind(str, 'X0')) || ~isempty(strfind(str, 'XNeg')) || ~isempty(strfind(str, 'XPos')) ||...
            ~isempty(strfind(str, 'Xneg')) || ~isempty(strfind(str, 'Xpos')) ||...
            ~isempty(strfind(str, 'N0')) || ~isempty(strfind(str, 'NNeg')) || ~isempty(strfind(str, 'NPos')) ||...
            ~isempty(strfind(str, 'Nneg')) || ~isempty(strfind(str, 'Npos'))
    %%%%%%%% use E1 parser
    
        if ~isempty(strfind(str, 'X0')) || ~isempty(strfind(str, 'N0'))
            pattern='^(?<phenotype>X0|N0)_(?<plate>UL\d*|UR\d*|LL\d*|LR\d*)_(?<well>\w_\d{1,2})';
            tokens=regexp(str,pattern,'names');
            tokens.MIC = 0;
        elseif ~isempty(strfind(str, 'pos')) || ~isempty(strfind(str, 'Pos'))
            pattern='^(?<phenotype>\w+)_(?<plate>UL\d*|UR\d*|LL\d*|LR\d*)_(?<well>\w_\d{1,2})';
            tokens=regexp(str,pattern,'names');
            tokens.MIC = 1;
        else
            pattern='^(?<phenotype>\w+)_(?<MIC>\d)_(?<plate>UL\d*|UR\d*|LL\d*|LR\d*)_(?<well>\w_\d{1,2})';
            tokens=regexp(str,pattern,'names');
        end
        
        tokens.well = strrep(tokens.well, '_', '');
        
        if isempty(tokens)
            disp('not found in E1:')
            disp(str)
        end
        
        if ~isempty(strfind(tokens.phenotype, 'X'))
            drug = 'OXA';
            database = database_oxa;
        elseif ~isempty(strfind(tokens.phenotype, 'N'))
            drug = 'NOV';
            database = database_nov;
        else
            disp(['drug not found' drug])
        end
        
        dataEntry = subsetDatabase(database, tokens);
        
        % fill in the strainStruct structure for passing to the Strain
        % Class constructor
        strainStruct = struct();
        strainStruct.name=str;
        strainStruct.library = dataEntry.drugname;
        strainStruct.quadrant = dataEntry.quadrant;
        strainStruct.wellname96 = dataEntry.wellName96;
        strainStruct.phenotype = upper(tokens.phenotype); % convert to all upper case
        strainStruct.conc = str2double(dataEntry.drugConc);
        strainStruct.experiment = 'E1';
        strainStruct.step = 1;
        strainStruct.culture = dataEntry.culture;
        strainStruct.platename = dataEntry.plateName;
        strainStruct.MIC = [dataEntry.micOXA_wt, dataEntry.micNOV_wt];
        strainStruct.mic_OXA = dataEntry.micOXA_wt;
        strainStruct.mic_NOV = dataEntry.micNOV_wt;
        strainStruct.mic_CPR = dataEntry.micCPR_wt;
        strainStruct.mic_DOX = dataEntry.micDOX_wt;
        strainStruct.mic_AMK = dataEntry.micAMK_wt;
        strainStruct.mic_GEN = dataEntry.micGEN_wt;
        
        STRAINS(i) = Strain(strainStruct);
        
    else
        %%%% use cycling parser
        if isempty(strfind(str, 'HL'))
            pattern = '^(?<culture>\w+[^_])_*(?<step>\d)';
            tokens=regexp(str,pattern,'names');
            step = str2double(tokens.step);
            
            
            if isempty(tokens)
                disp('not found for cycling:')
                disp(str)
            end
            
            row = cyclingMICs.(strrep(tokens.culture,'_', ''));
            row = row(step+1,:); %there's 00 at the beginning for plotting, that's why we need the +1
            
            strainStruct = struct();
            strainStruct.step = step;
            strainStruct.name=str;
            if strfind(str, 'X')
                strainStruct.library = 'OXA';
            elseif strfind(str, 'N')
                strainStruct.library = 'NOV';
            end
            strainStruct.culture = tokens.culture;
            strainStruct.experiment = 'Cycling';
            strainStruct.MIC = row(1:2);
            strainStruct.MICstd = row(3:4);
            strainStruct.mic_OXA = row(1);
            strainStruct.mic_NOV = row(2);
            
            if row(1) > 0 && row(2)>0
                strainStruct.phenotype='POS';
            elseif row(1) >0 && row(2)<0
                strainStruct.phenotype='XNEG';
            elseif row(1) >0 && row(2)==0
                strainStruct.phenotype='X0';
            elseif row(2) >0 && row(1)<0
                strainStruct.phenotype='NNEG';
            elseif row(2) >0 && row(1)==0
                strainStruct.phenotype='N0';
            else
                'here'
                strainStruct
            end
            
            
            STRAINS(i) = Strain(strainStruct);
            
        end
        
    end
    
    try
        strain_summary = summary.(strains{i});
        
        STRAINS(i).reads = strain_summary(1);
        STRAINS(i).align = strain_summary(2);
    catch
        disp('error in adding read summaries')
        i, strains{i}, str
    end
end
'strain obj done'
%% Create Event objects
allparams = struct();
allparams.strains = strains;
allparams.summaryReads = summaryReads;
allparams.summaryAlign = summaryAlign;
allparams.html = table2cell(html);
allparams.mutmat = mutmat;
allparams.param_consensus = param_consensus;
allparams.param_polymorphism = param_polymorphism;
allparams.param_strandBias = param_strandBias;
allparams.param_ksTest = param_ksTest;
allparams.param_biasE = param_biasE;
allparams.param_biasP = param_biasP;
allparams.param_majorAFreq = param_majorAFreq;
allparams.param_covTotal = param_covTotal;
allparams.param_covTop = param_covTop;
allparams.param_covBot = param_covBot;
allparams.param_freq = param_freq;

allparams.param_consensusFC = param_consensusFC;
allparams.param_polymorphismFC = param_polymorphismFC;
allparams.param_strandBiasFC = param_strandBiasFC;
allparams.param_ksTestFC = param_ksTestFC;
allparams.param_biasEFC = param_biasEFC;
allparams.param_biasPFC = param_biasPFC;
allparams.param_majorAFreqFC = param_majorAFreqFC;
allparams.param_covTotalFC = param_covTotalFC;
allparams.param_covTopFC = param_covTopFC;
allparams.param_covBotFC = param_covBotFC;
allparams.param_freqFC = param_freqFC;


EVENTMATRIX = Event.empty(size(mutmat,1), size(mutmat,2),0);
empty = Event();
for i = 1:size(mutmat,1)
    for j = 1:size(mutmat,2)
        if mutmat(i,j) == 1
            EVENTMATRIX(i,j,1) = Event(allparams, i,j, MUTATIONS(i), STRAINS(j));
            % also add mutation obj to strain obj and vice versa
            addMutations(STRAINS(j), MUTATIONS(i))
            addStrains(MUTATIONS(i), STRAINS(j))
            % note strains in mutations won't exactly match STRAINS--some of them
            % might have been filtered out due to low read counts earlier
        else
            EVENTMATRIX(i,j,1) = empty;
        end
    end
end
'event obj done'
%% Get rid of likely contaminations
EVENTMATRIX = EVENTMATRIX(:, ~strcmp({STRAINS.name},'XNeg_3_LR1_D_4'));
STRAINS = STRAINS(~strcmp({STRAINS.name},'XNeg_3_LR1_D_4'));
'elim contam done'

%% Add Event info to MUTATIONS and STRAINS
for i = 1:size(EVENTMATRIX,1)
    for j = 1:size(EVENTMATRIX,2)
        if EVENTMATRIX(i,j).call
            EVENTMATRIX(i,j).STRAIN.EVENTS = [EVENTMATRIX(i,j).STRAIN.EVENTS EVENTMATRIX(i,j)];
            EVENTMATRIX(i,j).MUTATION.EVENTS = [EVENTMATRIX(i,j).MUTATION.EVENTS EVENTMATRIX(i,j)];
        end
    end
end
'assign events done'

%% Export passed mutations to csv file for mpileup forcecalling
writetable(table(test.mutations.position, test.mutations.ref, test.mutations.alt),...
    'mutationlist_20160528.txt',...
    'Delimiter', ',',...
    'WriteVariableNames', 0,...
    'WriteRowNames',0);

%% View all breseq parameters in a Matlab GUI
viewer(allparams)



