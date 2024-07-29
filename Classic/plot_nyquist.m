function nyquist_data = plot_nyquist(sys, options)
    arguments
        sys
        options.clear_plot = true;
        options.plot = true
        options.plot_seprated = false
        options.plot_path = false
        options.plot_var = false
        
        options.omega_inf = [-200; 200]

        options.nOmega_infPath = 1e5 + 1
        options.nOmega_avoidPath = 1e5 + 1
        options.nOmega_imagPath = 1e5 + 1

        options.tol_abs_vanish = 1e-9
        options.tol_re_vanish = 1e-10

        options.turnsign { mustBeMember(options.turnsign , [1, -1])} = -1

        options.eps1 {mustBePositive} = 5e-4% Scale of Avoidance Curve
        options.eps2 {mustBeNonnegative} = 1e-4 % Scale of Asymmetry around root on Imag
        options.decatrate = 1.02

        options.avoidcurve {...
                            mustBeMember(options.avoidcurve, ...
                                {'circle', 'polar1', 'polar2'})} = 'circle'

        options.plot_turns = true
        options.phase_turncount = 3*pi/4
    end
    
    %% General Parameters
    omega_inf = options.omega_inf;    
    N_infPath = options.nOmega_infPath;
    N_avoidPath = options.nOmega_avoidPath;
    N_imagPath = options.nOmega_imagPath;

    tol_abs_vanish = options.tol_abs_vanish;
    tol_re_vanish = options.tol_re_vanish;

    eps1 = options.eps1;
    eps2 = options.eps2;

    decayrate = options.decatrate;

    phase_turncount = options.phase_turncount;
    turnsign = options.turnsign;
    avoidcurve = options.avoidcurve;

    %% System
    sys_zpk = zpk(sys);
    pole_arr = sys_zpk.P{1};
    zero_arr = sys_zpk.Z{1};
    num_deg = length(zero_arr);
    den_deg = length(pole_arr);    
    
    sys_tf = tf(sys);
    num_poly = sys_tf.Numerator{1};
    den_poly = sys_tf.Denominator{1};
    
    %% I) Extract Poles and Zeros in Origin, Called Origin Excluded as Pure
    zero_inOrig_idx = abs(zero_arr) < tol_abs_vanish ;
    pole_inOrig_idx = abs(pole_arr) < tol_abs_vanish ;
    zero_inOrig_arr = zero_arr(zero_inOrig_idx);
    pole_inOrig_arr = pole_arr(pole_inOrig_idx);

    d_num = length(zero_inOrig_arr);
    d_den = length(pole_inOrig_arr);
    if d_num ~= 0 && d_den ~= 0
        error('Zero Pole Cancellation at Orign')
    end
    d = d_den - d_num;
    
    num_deg_pure = num_deg - d_num;
    den_deg_pure = den_deg - d_den;
    num_poly_pure = num_poly(end - num_deg:end - d_num);
    den_poly_pure = den_poly(end - den_deg:end - d_den);
    
    num_pure_fun = @(s) sum(num_poly_pure.*(s.^(num_deg_pure:-1:0)), 2);
    den_pure_fun = @(s) sum(den_poly_pure.*(s.^(den_deg_pure:-1:0)), 2);
    
    G = @(s) (1./s.^d).*num_pure_fun(s)./den_pure_fun(s);

    zero_array_pure = zero_arr(~zero_inOrig_idx);
    pole_array_pure = pole_arr(~pole_inOrig_idx);

    for iZero_test = 1:length(zero_array_pure)
        zero_test = zero_array_pure(iZero_test);
        for jPole_test = 1:length(pole_array_pure)
            pole_test = pole_array_pure(jPole_test);
            if abs(zero_test - pole_test) < tol_abs_vanish
                error('Zero Pole Cancellation')
            end
        end
    end

    %% II) Cross Points On Real and Imaginary Axis
    syms omega_sym real
    assume(omega_sym >=0)
    
    num_pure_comp_sym = num_pure_fun(omega_sym*1j);
    num_pure_re_sym = real(num_pure_comp_sym);
    num_pure_im_sym = imag(num_pure_comp_sym);
    
    den_pure_comp_sym = den_pure_fun(omega_sym*1j);
    den_pure_re_sym = real(den_pure_comp_sym);
    den_pure_im_sym = imag(den_pure_comp_sym);

    alpha_sym = num_pure_re_sym*den_pure_re_sym ...
              + num_pure_im_sym*den_pure_im_sym;
    beta_sym = den_pure_re_sym*num_pure_im_sym ...
             - num_pure_re_sym*den_pure_im_sym;    

    G_prim_comp_sym = omega_sym^max(0, -d)...
                    *(alpha_sym + beta_sym*1j)/(1j)^d;
    G_prim_re_sym = real(G_prim_comp_sym);
    G_prim_im_sym = imag(G_prim_comp_sym);

    den_pure_abs_sym = den_pure_re_sym^2 + den_pure_im_sym^2;
    delta_sym = abs(omega_sym)^max(0, d)*den_pure_abs_sym;
    
    omega_imcross_arr = vpa(solve(G_prim_re_sym, omega_sym));
    G_im_imcross_arr = omega_imcross_arr*0;
    omega_crossim_idx = 0;
    while omega_crossim_idx < length(omega_imcross_arr)
        omega_crossim_idx = omega_crossim_idx + 1;
        omega_imcross = omega_imcross_arr(omega_crossim_idx);
        delta = vpa(subs(delta_sym, omega_sym, omega_imcross));
        if delta <= tol_abs_vanish
            omega_imcross_arr(omega_crossim_idx) = [];
            G_im_imcross_arr(omega_crossim_idx) = [];
            omega_crossim_idx = omega_crossim_idx - 1;
        else
            G_prim_im = subs(G_prim_im_sym, omega_sym, omega_imcross);
            G_im_imcross_arr(omega_crossim_idx) = vpa(G_prim_im/delta);
        end
    end

    omega_recross_arr = vpa(solve(G_prim_im_sym, omega_sym));
    G_re_recross_arr = omega_recross_arr*0;
    omega_crossre_idx = 0;
    while omega_crossre_idx < length(omega_recross_arr)
        omega_crossre_idx = omega_crossre_idx + 1;
        omega_recross = omega_recross_arr(omega_crossre_idx);
        delta = vpa(subs(delta_sym, omega_sym, omega_recross));
        if delta <= tol_abs_vanish
            omega_recross_arr(omega_crossre_idx) = [];
            G_re_recross_arr(omega_crossre_idx) = [];
            omega_crossre_idx = omega_crossre_idx - 1;
        else
            G_prim_re = subs(G_prim_re_sym, omega_sym, omega_recross);
            G_re_recross_arr(omega_crossre_idx) = vpa(G_prim_re/delta);
        end
    end

    omega_crtitical = vpasolve(G(omega_sym*1j) + 1 == 0, omega_sym);
    omega_crtitical(abs(imag(omega_crtitical)) > tol_re_vanish) = [];
    if ~isempty(omega_crtitical)
        crit_stable_pole = true;
    else
        crit_stable_pole = false;
    end

    G_comp_sym_data.alpha1_sym = num_pure_re_sym;
    G_comp_sym_data.beta1_sym = num_pure_im_sym;
    G_comp_sym_data.alpha2_sym = den_pure_re_sym;
    G_comp_sym_data.beta2_sym = den_pure_im_sym;
    G_comp_sym_data.G_prim_re_sym = G_prim_re_sym;
    G_comp_sym_data.G_prim_im_sym = G_prim_im_sym;
    G_comp_sym_data.delta_sym = delta_sym;
    G_comp_sym_data.omega_sym = omega_sym;

    %% III) Nyquist Path: 1. Modify Turns
    % zero_im_onPath_idx = abs(real(zero_arr)) < tol_re_vanish;
    pole_im_onPath_idx = abs(real(pole_arr)) < tol_re_vanish;
    % zero_im_onPath_arr = imag(zero_arr(zero_im_onPath_idx));
    % zero_im_onPath_arr = [];
    pole_im_onPath_arr = imag(pole_arr(pole_im_onPath_idx));
    
    pt_im_onPath_arr = sort(unique(pole_im_onPath_arr));
    nAvoids = length(pt_im_onPath_arr);
    nImSects = 1 + 2*nAvoids;
    omegaspan_arr = cell(nImSects + 1, 1);
    G_comp_arr = cell(nImSects + 1, 1);
    
    for i = 1:nAvoids
        p_im_onPath = pt_im_onPath_arr(i);
        [re, im] = get_avoidpath(p_im_onPath, eps1, eps2, ...
                                       N_avoidPath, avoidcurve, turnsign);
        eps1 = decayrate*eps1;
        eps2 = decayrate*eps2;
        omegaspan_sect = re + im*1i;
        omegaspan_arr{2*i} = omegaspan_sect;
        G_comp_arr{2*i} = G(omegaspan_sect);
    end

    P_nyq = sum(real(pole_array_pure)  >= tol_re_vanish)  ...
        + length(pole_im_onPath_arr)*max(0, options.turnsign);

    %% III) Nyquist Path: 2. Imaginary Axis
    omega_init = omega_inf(1);
    for i = 1:nAvoids
        omega_nextsection = omegaspan_arr{2*i};
        omega_fin = imag(omega_nextsection(1));
        omegaspan_sect = linspace(omega_init, omega_fin, N_imagPath)'*1i;
        omegaspan_arr{2*i - 1} = omegaspan_sect;
        G_comp_arr{2*i - 1} = G(omegaspan_sect);
        omega_init = imag(omega_nextsection(end));
    end
    omegaspan_sect = linspace(omega_init, omega_inf(2), N_imagPath)'*1i;
    omegaspan_arr{nImSects} = omegaspan_sect;
    G_comp_arr{nImSects} = G(omegaspan_sect);
                        
    %% III) Nyquist Path: 3. Infinite Path
    [re, im] = get_infpath(omega_inf, N_infPath);
    omegaspan_sect = re + im*1i;
    omegaspan_arr{nImSects + 1} = omegaspan_sect;
    G_comp_arr{nImSects + 1} = G(omegaspan_sect);

    %% III)  Nyquist Path: Gather Path Sections
    omegaspan = cell2mat(omegaspan_arr);
    G_comp = G(omegaspan);

    %% IV) Contact Points
    G_zegond = 1 + G_comp;
    G_zegond_phase = angle(G_zegond);
    contact_idx = find([0; diff(sign(G_zegond_phase - phase_turncount))]);
    contact_data = repmat(contact_idx, [1, 3]);
    if ~crit_stable_pole
        i = 0;
        while i < length(contact_idx)
            i = i + 1;
            idx = contact_idx(i);
            if abs(G_zegond_phase(idx) - phase_turncount) < pi/6 && ...
                idx ~= 0 && idx ~= length(omegaspan)
                d1 = G_zegond_phase(idx + 1) - G_zegond_phase(idx - 1);
                d2 = G_comp(idx + 4) - G_comp(idx - 4);
                if d1 > 0
                    contact_data(i, :) = [real(d2), imag(d2), -1];
                else
                    contact_data(i, :) = [real(d2), imag(d2), 1];
                end                
            else
                contact_idx(i) = [];
                contact_data(i, :) = [];
                i = i - 1;
            end
        end
        nContacts = i;
        N_nyq = sum(contact_data(:, 3));
        Z_nyq = N_nyq + P_nyq;        
    else
        nContacts = 0;
        contact_data = [];
        N_nyq = NaN;
        Z_nyq = NaN;  
    end
    

    %% Plots
    if options.clear_plot
        close all
    end

    if options.plot_path
        figure("Name", "Nyquist Path", "NumberTitle", "off")
        hold on, grid on, box on
        title(['Modified Nyquiest Path, P = ', num2str(P_nyq)])
        for i = 1:nImSects
            omegaspan_sect = omegaspan_arr{i};
            omegaspan_sect_re = real(omegaspan_sect);
            omegaspan_sect_im = imag(omegaspan_sect);
            if mod(i, 2) == 1
                plot(omegaspan_sect_re, omegaspan_sect_im, ...
                    'LineWidth', 1.5, 'Color', 'c', 'LineStyle', ':')
            else
                plot(omegaspan_sect_re, omegaspan_sect_im, ...
                    'LineWidth', 2, 'Color', 'b')
                txt_str = [num2str(pt_im_onPath_arr(i/2)), '^-'];
                text(omegaspan_sect_re(1), omegaspan_sect_im(1), ...
                    txt_str, 'FontSize',14)
                txt_str = [num2str(pt_im_onPath_arr(i/2)), '^+'];
                text(omegaspan_sect_re(end), omegaspan_sect_im(end), ...
                    txt_str, 'FontSize',14)
            end
            if i == 1
                text(omegaspan_sect_re(1), omegaspan_sect_im(1), ...
                    ['-\infty'])
            end
            if i == nImSects
                text(omegaspan_sect_re(end), omegaspan_sect_im(end), ...
                    ['\infty'])
            end    
        end
        omegaspan_sect = omegaspan_arr{end};
        omegaspan_sect_re = real(omegaspan_sect);
        omegaspan_sect_im = imag(omegaspan_sect);
        plot(omegaspan_sect_re, omegaspan_sect_im, ...
            'LineWidth', 1.5, 'Color', 'r', 'LineStyle', '-.')
        xlabel('Re'), ylabel('Im')
    end

    if options.plot  
        figure("Name", "Nyquist Diagram", "NumberTitle", "off")
        hold on, grid on, box on
        title(['Nyquist Diagram, ' ...
                'N=', num2str(N_nyq), [', ' ...
                'P='], num2str(P_nyq), [' ==> ' ...
                'Z=N+P='], num2str(Z_nyq)])
        for i = 1:nImSects
            G_comp_sect = G_comp_arr{i};
            G_re_sect = real(G_comp_sect);
            G_im_sect = imag(G_comp_sect);
            if mod(i, 2) == 1 % On Imaginary Axis
                plot(G_re_sect, G_im_sect, 'LineWidth', 1.5, ...
                    'Color', 'c', 'LineStyle', ':')
            else
                plot(G_re_sect, G_im_sect, 'LineWidth', 2, 'Color', 'b')
                plot(G_re_sect, G_im_sect, 'LineWidth', 2, 'Color', 'b')
                txt_str = [num2str(pt_im_onPath_arr(i/2)), '^-'];
                text(G_re_sect(1), G_im_sect(1), txt_str, 'FontSize', 14)
                txt_str = [num2str(pt_im_onPath_arr(i/2)), '^+'];
                text(G_re_sect(end), G_im_sect(end), txt_str, ...
                                                          'FontSize', 14)        
            end
            if i == 1
                text(G_re_sect(1), G_im_sect(1), ['-\infty'])
            end
            if i == nImSects
                 text(G_re_sect(end), G_im_sect(end), ['\infty'])
            end 
        end
        omegaspan_sect = omegaspan_arr{end};
        G_comp_sect = G(omegaspan_sect);
        G_re_sect = real(G_comp_sect);
        G_im_sect = imag(G_comp_sect);
        plot(G_re_sect, G_im_sect, ...
            'LineWidth', 1.5, 'Color', 'r', 'LineStyle', '-.')
        xlabel('Re'), ylabel('Im')
        if options.plot_turns
            for i = 1:nContacts
                G_contact = G(omegaspan(contact_idx(i)));
                re = real(G_contact);
                im = imag(G_contact);
                plot(re, im, 'LineStyle', 'none', ...
                    'Marker', '*', 'MarkerSize', 10, 'LineWidth', 2)
                quiver(re, im, contact_data(i, 1), contact_data(i, 2), ...
                                                        'LineWidth', 3)
            end
        end
    end
    

    if options.plot_seprated
        prc_1 = ceil(-log10(eps1));
        prc_2 = prc_1 + 1;
        nExtraPlots = (nImSects - 1)/2;
        omega_init = '-\infty';
        txt_base = 'G(j \omega) \omega \in (';
        if nExtraPlots > 0
            figure("Name", "Nyquist Diagram Seperated (Negative Imaginary)", ...
                   "NumberTitle", "off")
        end
        for i = 1:nExtraPlots
            subplot(1, nExtraPlots, i)
            hold on, grid on, box on
            omegaspan_sect = omegaspan_arr{i};
            G_comp_sect = G(omegaspan_sect);
            G_re_sect = real(G_comp_sect);
            G_im_sect = imag(G_comp_sect);
            if mod(i, 2) == 1  % On Imaginary Axis
                col = 'c';
                omega_fin = [num2str(pt_im_onPath_arr((i+1)/2)), '^-'];
            else
                col = 'b';
                omega_fin = [num2str(pt_im_onPath_arr(i/2)), '^+'];
            end
            txt = [txt_base, omega_init, ',', omega_fin, ')'];
            title(txt)
            omega_init = omega_fin;
                             
            plot(G_re_sect, G_im_sect, ...
                   'LineWidth', 1.5, 'Color', col, 'LineStyle', '-.')
            txt_str = ['\omega_{start}=', ...
                num2str(round(omegaspan_sect(1), prc_1), prc_2)];
            text(G_re_sect(1), G_im_sect(1),txt_str)
            txt_str = ['\omega_{end}=', ...
                num2str(round(omegaspan_sect(end), prc_1), prc_2)];
            text(G_re_sect(end), G_im_sect(end), txt_str)
            xlabel('Re'), ylabel('Im')
            axis square
        end
        
        figure("Name", "Nyquist Diagram Seperated (Middle)", ...
            "NumberTitle", "off")
        hold on, grid on, box on
        i = i + 1;
        if nExtraPlots == 0
            omega_fin = '\infty';
        else
            if mod(i, 2) == 1  % On Imaginary Axis
                col = 'c';
                omega_fin = [num2str(pt_im_onPath_arr((i+1)/2)), '^-'];
            else
                col = 'b';
                omega_fin = [num2str(pt_im_onPath_arr(i/2)), '^+'];
            end
        end
        
        txt = [txt_base, omega_init, ',', omega_fin, ')'];
        title(txt)
        omega_init = omega_fin;
        omegaspan_sect = omegaspan_arr{nExtraPlots + 1};
        G_comp_sect = G(omegaspan_sect);
        G_re_sect = real(G_comp_sect);
        G_im_sect = imag(G_comp_sect);
        plot(G_re_sect, G_im_sect, ...
               'LineWidth', 1.5, 'Color', 'c', 'LineStyle', '-.')
        txt_str = ['\omega_{start}=', ...
            num2str(round(omegaspan_sect(1), prc_1), prc_2)];
        text(G_re_sect(1), G_im_sect(1),txt_str)
        txt_str = ['\omega_{end}=', ...
            num2str(round(omegaspan_sect(end), prc_1), prc_2)];
        text(G_re_sect(end), G_im_sect(end), txt_str)
        xlabel('Re'), ylabel('Im')
        axis square

        if nExtraPlots > 0
            figure("Name", "Nyquist Diagram Seperated (Positive Imaginary)", ...
                   "NumberTitle", "off")
        end
        for i = nExtraPlots + 2:2*nExtraPlots + 1
            subplot(1, nExtraPlots, i - nExtraPlots - 1)
            hold on, grid on, box on
            omegaspan_sect = omegaspan_arr{i};
            G_comp_sect = G(omegaspan_sect);
            G_re_sect = real(G_comp_sect);
            G_im_sect = imag(G_comp_sect);
            if i == 2*nExtraPlots + 1
                omega_fin = '+\infty';
            else
                if mod(i, 2) == 1
                    col = 'c';
                    omega_fin = [num2str(pt_im_onPath_arr((i+1)/2)), '^-'];
                else
                    col = 'b';
                    omega_fin = [num2str(pt_im_onPath_arr(i/2)), '^+'];
                end            
            end

            txt = [txt_base, omega_init, ',', omega_fin, ')'];
            title(txt)
            omega_init = omega_fin;
            plot(G_re_sect, G_im_sect, ...
                   'LineWidth', 1.5, 'Color',  col, 'LineStyle', '-.')
            txt_str = ['\omega_{start}=', ...
                num2str(round(omegaspan_sect(1), prc_1), prc_2)];
            text(G_re_sect(1), G_im_sect(1),txt_str)
            txt_str = ['\omega_{end}=', ...
                num2str(round(omegaspan_sect(end), prc_1), prc_2)];
            text(G_re_sect(end), G_im_sect(end), txt_str)
            xlabel('Re'), ylabel('Im')
            axis square
        end
    end

    if options.plot_var
        figure("Name", "Nyquist Diagram (system)", "NumberTitle", "off")
        nyquist(sys)        
    end

    %% Gather Output
    nyquist_data.pole_im_onPath_arr = pole_im_onPath_arr;
    
    nyquist_data.omegaspan_arr = omegaspan_arr;
    nyquist_data.G_comp_arr = G_comp_arr;

    nyquist_data.omegaspan = omegaspan;
    nyquist_data.G_comp = G_comp;
    
    nyquist_data.omega_imcross_arr = omega_imcross_arr;
    nyquist_data.G_im_imcross_arr = G_im_imcross_arr;
    nyquist_data.omega_recross_arr = omega_recross_arr;
    nyquist_data.G_re_recross_arr = G_re_recross_arr;

    nyquist_data.P_nyq = P_nyq;
    nyquist_data.N_nyq = N_nyq;
    nyquist_data.Z_nyq = Z_nyq;
    nyquist_data.crit_stable_pole = crit_stable_pole;
    nyquist_data.contact_data = contact_data;

    nyquist_data.G_comp_sym_data = G_comp_sym_data;

end

function [re, im] = get_avoidpath(p_im_onPath, epsilon1, epsilon2, N, ...
                                  type, turnsign)
    if epsilon2 >= epsilon1/2
        error('Wrong epsilon values')
    end
    thetaspan = linspace(-pi/2, pi/2, N)';
    switch type
        case 'circle'
            re_curve = -cos(thetaspan);
            im_curve = sin(thetaspan);            
        case 'polar1'
            
            r_fun = @(theta) 1 + cos(theta);
            re_curve = -r_fun(thetaspan).*cos(thetaspan);
            im_curve = r_fun(thetaspan).*sin(thetaspan);            
        case 'polar2'
            thetaspan = linspace(-pi/2, pi/2, N)';
            r_fun = @(theta) 1 + sin(theta) + cos(theta);
            re_curve = -r_fun(thetaspan).*cos(thetaspan);
            im_curve = r_fun(thetaspan).*sin(thetaspan) - 1;            
        otherwise
            error('Wrong Type')
    end
    re = epsilon1*re_curve*turnsign;
    im = epsilon1*im_curve + epsilon2 + p_im_onPath;
end

function [re, im] = get_infpath(omega_inf, N)
    thetaspan= linspace(pi/2, -pi/2, N)';
    center = (omega_inf(2) + omega_inf(1))/2;
    radius = (omega_inf(2) - omega_inf(1))/2;
    re = radius*(cos(thetaspan));
    im = center + radius*(sin(thetaspan));
end
