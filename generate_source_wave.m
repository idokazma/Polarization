if (params.tm)
    
    
    %E_inc_z = exp(-1i*k0*params.X);
    params.E_inc_z = scalar_green([params.source_loc_x;params.source_loc_y],[params.X(:),params.Y(:)]',params,0);
    params.E_inc_z = reshape(params.E_inc_z,length(params.X),[]);
    params.E_inc_x = 0*params.E_inc_z;
    params.E_inc_y = 0*params.E_inc_z;
    
    %     H_inc_y = 1i/(omega*(mu0*mr_out))*(-1i*k0*exp(-1i*k0*params.X));
    %     H_inc_x = 0*(H_inc_y);
    %     H_inc_z = 0*(H_inc_y);
    [params.H_inc_x,params.H_inc_y] = dyiadic_green([params.source_loc_x;params.source_loc_y],[params.X(:),params.Y(:)]',params,0);
    params.H_inc_x = reshape(1*params.H_inc_x,length(params.X),[]);
    params.H_inc_y = reshape(1*params.H_inc_y,length(params.X),[]);
    params.H_inc_z = 0*params.H_inc_y;
    
    if (false)
    %symbolic
    syms xsym ysym zsym k0sym
    params.Xsym = [xsym ysym zsym];
    %     Esym = [0, 0,exp(-1i*k0sym*xsym)]; % plane wave
    %     Xsym = [xsym ysym zsym];
    %     Hsym = -1i/(omega*(mu0*mr_out))*curl(Esym,Xsym);
    % end
    ez = 1i/4 * besselh (0 , 1 , k0sym *sqrt((xsym - params.source_loc_x).^2+(ysym- params.source_loc_y).^2));
    grot =  exp(1i*params.k0*params.OMEGA/params.c*(params.source_loc_x.*ysym-params.source_loc_y.*xsym));
    params.Esym = [0, 0, ez*grot];
    %     params.Hsym = 1/(1i*params.omega*(params.mu0*params.mr_out))*curl(params.Esym,params.Xsym) - transpose(params.Xsym.*params.Esym(3)*params.OMEGA/(params.c^2*(params.mu0*params.mr_out))) ;
    
    params.Hsym_curl =curl(params.Esym,params.Xsym);
    %     params.transpose = -transpose(params.Xsym.*params.Esym(3)*params.OMEGA/(params.c^2*(params.mu0*params.mr_out)));
    params.transpose = (-1*(1i*params.Xsym.*params.omega.*params.OMEGA/(params.c^2))).';
    
    params.Hsym = 1/(1i*params.omega*(params.mu0*params.mr_out)) * (params.Hsym_curl + params.Esym(3).*params.transpose);
    end
    if (false)
        tic
        E_inc = zeros(length(params.X(:)),3);
        H_inc = E_inc;
        for i = 1:length(params.X(:))
            E_inc(i,:) = double(subs(params.Esym,[xsym,ysym,k0sym] ,[params.X(i),params.Y(i),params.k0]));
            %     H_inc(i,:) = double(subs(params.Hsym,[xsym,ysym,zsym,k0sym] ,[params.X(i),params.Y(i),0,params.k0]));
            H_inc(i,:) = double(subs(params.Hsym,[xsym,ysym,zsym,k0sym] ,[params.X(i),params.Y(i),0,params.k0]));
        end
        toc
        check_symbolic_vs_green;
    end
end

if (params.te)
    
    %E_inc_z = exp(-1i*k0*params.X);
    params.H_inc_z = scalar_green([params.source_loc_x;params.source_loc_y],[params.X(:),params.Y(:)]',params.n_out,params.k0,params.c,params.OMEGA);
    params.H_inc_z = reshape(params.H_inc_z,length(params.X),[]);
    params.H_inc_x = 0*params.H_inc_z;
    params.H_inc_y = 0*params.H_inc_z;
    
    %     H_inc_y = 1i/(omega*(mu0*mr_out))*(-1i*k0*exp(-1i*k0*params.X));
    %     H_inc_x = 0*(H_inc_y);
    %     H_inc_z = 0*(H_inc_y);
    [params.E_inc_x,params.E_inc_y] = dyiadic_green([params.source_loc_x;params.source_loc_y],[params.X(:),params.Y(:)]',params.n_out,params.k0, c,params.OMEGA);
    params.E_inc_x = reshape(1i/(params.omega*(params.e0*params.er_out))*params.E_inc_x,length(params.X),[]);
    params.E_inc_y = reshape(1i/(params.omega*(params.e0*params.er_out))*params.E_inc_y,length(params.X),[]);
    params.E_inc_z = 0*params.E_inc_y;
    
    if (false)
        %symbolic
        syms xsym ysym zsym k0sym
        params.Xsym = [xsym ysym zsym];
        %     Esym = [0, 0,exp(-1i*k0sym*xsym)]; % plane wave
        %     Xsym = [xsym ysym zsym];
        %     Hsym = -1i/(omega*(mu0*mr_out))*curl(Esym,Xsym);
        % end
        hz = 1i/4 * besselh (0 , 1 , k0sym *sqrt((xsym - params.source_loc_x).^2+(ysym- params.source_loc_y).^2));
        grot =  exp(1i*params.k0*params.OMEGA/params.c*(params.source_loc_x.*ysym-params.source_loc_y.*xsym));
        params.Hsym = [0, 0, hz*grot];
        
        params.Esym_curl =curl(params.Hsym,params.Xsym);
        params.transpose = (1*(1i*params.Xsym.*params.omega.*params.OMEGA/(params.c^2))).';
        
        params.Esym = 1/(1i*params.omega*(params.e0*params.er_out)) * (-1*params.Esym_curl + params.Hsym(3).*params.transpose);
        if (false)
            tic
            E_inc = zeros(length(params.X(:)),3);
            H_inc = E_inc;
            for i = 1:length(params.X(:))
                E_inc(i,:) = double(subs(params.Esym,[xsym,ysym,zsym,k0sym] ,[params.X(i),params.Y(i),0,params.k0]));
                %     H_inc(i,:) = double(subs(params.Hsym,[xsym,ysym,zsym,k0sym] ,[params.X(i),params.Y(i),0,params.k0]));
                H_inc(i,:) = double(subs(params.Hsym,[xsym,ysym,k0sym] ,[params.X(i),params.Y(i),params.k0]));
            end
            toc
            
            check_symbolic_vs_green_TE;
        end
    end
    
end