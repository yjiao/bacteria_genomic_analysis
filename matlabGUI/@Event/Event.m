% http://www.mathworks.com/help/matlab/matlab_oop/user-defined-classes.html
% use a handle class here because we never want more than one copy of each
% mutation
classdef Event < handle
    properties
        call = 0;
        html = '';
        param_consensus
        param_polymorphism
        param_strandBias
        param_ksTest
        param_biasE
        param_biasP
        param_majorAFreq
        param_covTotal
        param_covTop
        param_covBot
        
        param_consensusFC
        param_polymorphismFC
        param_strandBiasFC
        param_ksTestFC
        param_biasEFC
        param_biasPFC
        param_majorAFreqFC
        param_covTotalFC
        param_covTopFC
        param_covBotFC
        
        param_freq
        param_freqFC
        
        STRAIN
        MUTATION
    end
    
    methods
        function obj = Event(varargin)
            if nargin > 0
                allparams = varargin{1};
                i = varargin{2};
                j = varargin{3};
                obj.MUTATION = varargin{4};
                obj.STRAIN = varargin{5};
                
                obj.call = allparams.mutmat(i,j);
                obj.html = allparams.html(i,j);
                obj.param_consensus = allparams.param_consensus(i,j);
                obj.param_polymorphism = allparams.param_polymorphism(i,j);
                obj.param_strandBias = allparams.param_strandBias(i,j);
                obj.param_ksTest = allparams.param_ksTest(i,j);
                obj.param_biasE = allparams.param_biasE(i,j);
                obj.param_biasP = allparams.param_biasP(i,j);
                obj.param_majorAFreq = allparams.param_majorAFreq(i,j);
                obj.param_covTotal = allparams.param_covTotal(i,j);
                obj.param_covTop = allparams.param_covTop(i,j);
                obj.param_covBot = allparams.param_covBot(i,j);
                obj.param_consensusFC = allparams.param_consensusFC(i,j);
                obj.param_polymorphismFC = allparams.param_polymorphismFC(i,j);
                obj.param_strandBiasFC = allparams.param_strandBiasFC(i,j);
                obj.param_ksTestFC = allparams.param_ksTestFC(i,j);
                obj.param_biasEFC = allparams.param_biasEFC(i,j);
                obj.param_biasPFC = allparams.param_biasPFC(i,j);
                obj.param_majorAFreqFC = allparams.param_majorAFreqFC(i,j);
                obj.param_covTotalFC = allparams.param_covTotalFC(i,j);
                obj.param_covTopFC = allparams.param_covTopFC(i,j);
                obj.param_covBotFC = allparams.param_covBotFC(i,j);
                
                try
                    obj.param_freq = allparams.param_freq(i,j);
                    obj.param_freqFC = allparams.param_freqFC(i,j);
                end
                
            end
                
        end
    end
end