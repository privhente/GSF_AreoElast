function [ofast_res,strInpVar] = fun_LoadOpenFastResFile(strfilename)
    ofast_res = [];
    q   = [];
    
    fid = fopen(strfilename);
    if fid~= -1
        ofast_res = cell2mat(textscan(fid,'','headerlines',8));
        fclose(fid);
    end
    
    strInpVar = {};
    fid = fopen(strfilename);
    if fid~= -1
        for i = 1:6; strval = fgetl(fid); end
        strval = fgetl(fid);
        strInpVar = textscan(strval,'%s'); 
        fclose(fid);
    end