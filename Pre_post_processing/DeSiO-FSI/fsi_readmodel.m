function model = fsi_readmodel()
% =================================================================================================================
model = [];
fid1 = fopen(['simulationinput_fsi.txt'],'r');
fid2 = fopen(['simulationinput.txt'],'r');
fid = -1;
    if fid1 ~= -1
        fid = fid1;
    elseif fid2 ~= -1
        fid = fid2;
    end
    if fid >= 0 
        for i = 1:3
            tline = fgets(fid);
        end
        tline = fgetl(fid);
        [str_n,count] = sscanf(tline,'%s%c');
        str = str_n;
        model.strSimName.fsi       = '';
        model.strSimName.aero      = '';
        model.strSimName.structure = '';
        if count>1
            str_n = textscan(str_n,'%s');
            str = cellstr(str_n{1}(1)); str = str{1};
            model.strSimName = str;
%             str = cellstr(str_n{1}(2)); str = str{1};
%             model.strSimName.aero      = str;
%             str = cellstr(str_n{1}(3)); str = str{1};
%             model.strSimName.structure = str;
        else
            model.strSimName = str;
        end
        fclose(fid);
    end
fid1 = fopen([model.strSimName '_uvlm_models.dres'],'r');
fid2 = fopen([model.strSimName '_uvlm_models.txt'],'r');
fid = -1;
    if fid1 ~= -1
        fid = fid1;
    elseif fid2 ~= -1
        fid = fid2;
    end
if fid >= 0 
    tline = fgetl(fid);
    model.nsurfaces = sscanf(tline,'%i');
    i3 = 0;
    for i = 1:model.nsurfaces
        tline = fgetl(fid);
        arrline = sscanf(tline,'%i',5);
        model.surfaces(i).nnodes     = arrline(1);
        model.surfaces(i).nelements = arrline(4);
        for j = 1:model.surfaces(i).nelements
            tline   = fgetl(fid);
            model.surfaces(i).connectivity(j,:) = sscanf(tline,'%i',4);
        end
        for j = 1:model.surfaces(i).nnodes
            model.surfaces(i).nodes(j).indices_q = [3*(j-1)+1:3*(j-1)+3] + i3;
        end
        i3 = i3 + 3*model.surfaces(i).nnodes;
    end
fclose(fid);
end
% =================================================================================================================
return