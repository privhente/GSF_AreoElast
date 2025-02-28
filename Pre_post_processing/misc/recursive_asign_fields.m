function [model] = recursive_asign_fields(field,pa_field,nn,model)
    eval([pa_field '= field;']);
    inf_field = whos('pa_field');
    eval(['str_field   = ' inf_field.name ';']);
    fields = fieldnames(field);
    for i = 1:size(fields,1)
        eval(['str_field_t = ' '''' pa_field '.' fields{i} '''' ';']);
        eval(['temp_field=' pa_field '.' fields{i} ';' ]);
        inf_field_i = whos('temp_field');
        if strcmp(inf_field_i.class,'struct')
            str_field = [str_field_t];
            [model] = recursive_asign_fields(temp_field,str_field,nn,model);
        else
            eval(['model.' str_field_t(nn+2:end) '=' str_field_t ';']);
        end
    end
end