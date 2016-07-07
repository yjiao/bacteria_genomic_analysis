classdef Strain < handle
    properties %(Access = private)
        name='';
        library='' % which drug it was selected in, corresponds to "drugname" in mutdatabase
        quadrant = '' % UL1/ LR2 etc
        wellname96 = '' % 96 well into which it was picked
        phenotype = ''
        conc = -1
        experiment = ''
        step = -1
        
        
        culture = ''
        platename = ''
        
        MIC = nan
        MICstd = [-1 -1];
        MUTATIONS=[];
        EVENTS=[];
        reads = -1
        align = -1
        
        mic_CPR = nan
        mic_DOX = nan
        mic_OXA = nan
        mic_NOV = nan
        mic_AMK = nan
        mic_GEN = nan
    end
    
    methods
        function obj = Strain(strainStruct)
            if nargin == 1
                
                % mandatory fields
                obj.name = strainStruct.name;
                obj.mic_OXA = strainStruct.mic_OXA;
                obj.mic_NOV = strainStruct.mic_NOV;
                obj.MIC = strainStruct.MIC;
                
                % Optional Fields
                
                if isfield(strainStruct, 'wellname96')
                    obj.wellname96 = strainStruct.wellname96;
                end
                
                if isfield(strainStruct, 'conc')
                    obj.conc = strainStruct.conc;
                end
                
                
                if isfield(strainStruct, 'phenotype')
                    obj.phenotype = strainStruct.phenotype;
                end
                
                if isfield(strainStruct, 'quadrant')
                    obj.quadrant = strainStruct.quadrant;
                end
                
                if isfield(strainStruct, 'experiment')
                    obj.experiment = strainStruct.experiment;
                end
                
                if isfield(strainStruct, 'library')
                    obj.library = strainStruct.library;
                end
                
                if isfield(strainStruct, 'culture')
                    obj.culture = strainStruct.culture;
                end
                
                if isfield(strainStruct, 'platename')
                    obj.platename = strainStruct.platename;
                end
                
                if isfield(strainStruct, 'MICstd')
                    obj.MICstd = strainStruct.MICstd;
                end
                
                if isfield(strainStruct, 'MUTATIONS')
                    obj.MUTATIONS = strainStruct.MUTATIONS;
                end
                
                if isfield(strainStruct, 'mic_CPR')
                    obj.mic_CPR = strainStruct.mic_CPR;
                end
                
                if isfield(strainStruct, 'mic_DOX')
                    obj.mic_DOX = strainStruct.mic_DOX;
                end
                
                if isfield(strainStruct, 'mic_AMK')
                    obj.mic_AMK = strainStruct.mic_AMK;
                end
                
                if isfield(strainStruct, 'mic_GEN')
                    obj.mic_GEN = strainStruct.mic_GEN;
                end
                
                if isfield(strainStruct, 'step')
                    obj.step = strainStruct.step;
                end
                
            end
        end
        
        % Input = array of Mutation objects
        % Add Mutation obj to strain, also update parent property of each Mutation object
        function addMutations(this,mutationsArray)
            this.MUTATIONS=[this.MUTATIONS mutationsArray];
        end
%         
%         % display MICs of a strain collection as a nx2 matrix, col1=oxa, col2=nov
%         function mic=prettyMIC(sarray)
%             mic=cell2mat({sarray.MIC}');
%         end
%         
%         % Input: collection of Strain objects, all mutations
%         % Output: table summary of strain objects
%         tab=summarize(sarray,SNP,INDEL,varargin)
%         plotSNPs(sarray,varargin)
    end
end
