function [M] = uvlm_m_surface(model,t,qs,dp,x0)
% =================================================================================================================
    M = [];
    for k = 1:length(t)
        a = 0;
        mcp = [0,0,0];
        for i = 1:model.nsurfaces
            nodes = model.surfaces(i).nodes;
            for j = 1:model.surfaces(i).nelements
                c   = model.surfaces(i).connectivity(j,:);
                n1  = c(1); n2 = c(2); n3 = c(3); n4 = c(4);
                % Area of surface
                coord_n1 = qs(k,nodes(n1).indices_q);
                coord_n2 = qs(k,nodes(n2).indices_q);
                coord_n3 = qs(k,nodes(n3).indices_q);
                coord_n4 = qs(k,nodes(n4).indices_q);

                d21 = coord_n2-coord_n1;
                d23 = coord_n2-coord_n3;
                d41 = coord_n4-coord_n1;
                d43 = coord_n4-coord_n3;
                d24 = coord_n2-coord_n4;
                d13 = coord_n1-coord_n3;
                
                L1  = cross(d21,d41);
                L2  = cross(d43,d23);
                Aj = 0.5 * (norm(L1) + norm(L2));
                
                % calculating unit vectors of (trapezodial) ring surface
                h1 = d24-d13;
                h2 = -(d13+d24);
                e1 = h1/norm(h1);
                e2 = h2/norm(h2);
    
                % calculating plane normal of ring surface
                nj = cross(e1,e2);

                % coordinates of control point
                xcpj = 0.25*( coord_n1 + coord_n2 + coord_n3 + coord_n4 )';
                xaj  = xcpj - x0;
                % resultant force/dynamic pressure due to pressure
                Fj  = dp(k,j+a)*Aj*nj';
                mj  = cross(xaj,Fj)';
                mcp = mcp + mj;
            end
            % global lift = projecting Cp_vector in lift direction
            a = a + model.surfaces(i).nelements;
        end
        M(k,1:3) = mcp;
    end
% =================================================================================================================
return