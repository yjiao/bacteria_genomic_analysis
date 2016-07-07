% http://www.mathworks.com/help/matlab/matlab_oop/user-defined-classes.html
% use a handle class here because we never want more than one copy of each
% mutation
classdef Mutation < handle
    properties
        snp_type
        aa_ref_seq
        alt
        size
        locus_tag
        repeat_length
        gene_strand
        codon_position
        type
        repeat_new_copies
        aa_position
        aa_new_seq
        codon_ref_seq
        repeat_ref_copies
        ID
        gene_position
        ref
        codon_new_seq
        strains % strain names
        STRAINS % Strain objects
        gene_product
        position
        gene_list
        gene_name
        repeat_seq
        
        count
        
        description
        KeggOntology
        KeggA
        KeggB
        KeggC
        
        EVENTS=[];
        
    end
    
    methods
        function obj = Mutation(row, keggStruct)
            obj.snp_type = row.snp_type;
            obj.aa_ref_seq = row.aa_ref_seq;
            obj.alt = row.alt;
            obj.size = row.size;
            obj.locus_tag = row.locus_tag;
            obj.repeat_length = row.repeat_length;
            obj.gene_strand = row.gene_strand;
            obj.codon_position = row.codon_position;
            obj.type = row.type;
            obj.repeat_new_copies = row.repeat_new_copies;
            obj.aa_position = row.aa_position;
            obj.aa_new_seq = row.aa_new_seq;
            obj.codon_ref_seq = row.codon_ref_seq;
            obj.repeat_ref_copies = row.repeat_ref_copies;
            obj.ID = row.ID;
            obj.gene_position = row.gene_position;
            obj.ref = row.ref;
            obj.codon_new_seq = row.codon_new_seq;
            obj.strains = strsplit(row.strains{:},',');
            obj.count = length(obj.strains);
            obj.gene_product = row.gene_product;
            obj.position = row.position;
            obj.gene_list = row.gene_list;
            obj.gene_name = row.gene_name;
            obj.repeat_seq = row.repeat_seq;
            
            if ~isempty(strfind(obj.gene_name{1}, 'SAOUHSC')) && isKey(keggStruct.dict_Locus2Gene,obj.locus_tag{1})
                obj.gene_name = keggStruct.dict_Locus2Gene(obj.locus_tag{1});
            end
            
            if isKey(keggStruct.dict_Locus2Gene,obj.locus_tag{1})
                obj.description = keggStruct.dict_Locus2Description(obj.locus_tag{1});
                obj.KeggOntology = keggStruct.dict_Locus2KO(obj.locus_tag{1});
                obj.KeggA = keggStruct.dict_Locus2A(obj.locus_tag{1});
                obj.KeggB = keggStruct.dict_Locus2B(obj.locus_tag{1});
                obj.KeggC = keggStruct.dict_Locus2C(obj.locus_tag{1});
            end
        end
        
        function addStrains(this,StrainObj)
            this.STRAINS=[this.STRAINS StrainObj];
        end
    end
end