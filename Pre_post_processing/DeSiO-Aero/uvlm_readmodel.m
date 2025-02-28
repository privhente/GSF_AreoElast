function model = uvlm_readmodel()
% =================================================================================================================
model = [];
fid1 = -1;
currDirr = pwd;
if exist([currDirr '\simulationinput_aero.txt'], 'file')
    fid1 = fopen(['simulationinput_aero.txt'],'r');
end
fid2 = -1;
if exist([currDirr '\simulationinput_fsi.txt'], 'file')
    fid2 = fopen(['simulationinput_fsi.txt'],'r');
end
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
        if count>1
            str_n = textscan(str_n,'%s');
            str = cellstr(str_n{1}(1));
            str = str{1};
        end
        model.strSimName = str;
        for i = 1:3
            tline = fgets(fid);
        end
        tline = fgetl(fid); tline = strrep(tline,'d','e');
        model.simulationsettings = sscanf(tline,'%f',3);
        for i = 1:3
            tline = fgets(fid);
        end
        tline = fgetl(fid); tline = strrep(tline(10:end),'d','e');
        model.windData = sscanf(tline,'%f',6);
        model.density  = model.windData(1);
        model.winddir  = model.windData(4:6);
        model.vinf     = model.windData(2)*model.winddir;
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

% model_inp = model;
% model_inp.surfaces = uvlm_readsurfaceinput();
% inf_model_inp = whos('model_inp');
% [model] = recursive_asign_fields(model_inp,inf_model_inp.name,length(inf_model_inp.name),model);

% =================================================================================================================
return