function [routhMat, routhTable] =  routh(polynomial, opt)
    arguments
        polynomial (1, :)
        opt.Tol = 1e-10;
        opt.Verbos = false;
        opt.simplify = true;
        opt.limit = false;
    end
    syms s    
    deg = size(polynomial, 2) - 1;
    e = zeros(deg + 1, 1, 'sym');
    for i = 0:deg
        e(i + 1) = sym(['e', num2str(i)]);
        assume([0 < e(i + 1); e(i + 1) < opt.Tol])
    end    

    %% Manage Dimensions of Routh Table
    if rem(deg, 2) == 1
        tableWidth = (deg + 1)/2;
        secRowFill = [];         
    else
        tableWidth = deg/2 + 1;
        secRowFill = 0;     
    end
    routhMat = zeros(deg + 1, tableWidth, 'sym');
    descList = cell(deg + 1, 1);
    powerList = cell(deg + 1, 1);

    %% Fill Initial Values
    for i = 0:deg
        descList{i + 1, 1} = "";
        powerList{i + 1, 1} = char(s^(deg - i));
    end
    routhMat(1, 1:tableWidth) = polynomial(1:2:end);
    routhMat(2, 1:tableWidth) = [polynomial(2:2:end), secRowFill];     
    epsSubed = false;

    %% Handle Second Row Special Cases   
    zeroRow = false;
    zeroElm = false;
    if isempty(symvar(routhMat(2, :)))
        checkfun = @(item) all(abs(item) < opt.Tol);
    else
        checkfun = @(item) all(item == 0);
    end
    if checkfun(routhMat(2, 1))
        if checkfun(routhMat(2, 2:tableWidth))
            zeroElm = false;
            zeroRow = true;
        else
            zeroElm = true;
            epsSubed = true;
        end
    end
    if zeroRow
        charEq = zeros(1, 'sym');
        for j = deg:-2:0
            charEq = charEq + routhMat(1, (deg - j)/2 + 1)*s^(j);
        end
        if opt.Verbos
            fprintf('Even Polynomial Found as:\n')
            disp(charEq)
            fprintf('Roots:\n')
            disp(solve(charEq, s))                
        end
        diffAuxSym = diff(charEq, s); 
        descList{2, 1} = "Derivative Substituted";
        diffAuxpoly = coeffs(diffAuxSym, s, 'all');
        routhMat(2, 1:deg/2) = diffAuxpoly(1:2:end - 1);
    end
    if zeroElm
        routhMat(2, 1) = e(2);
        epsSubed = true;
        descList{2, 1} = 'Epsilon Substituted';
    end
    
    %% Fill All Table
    for i = 3:deg + 1
        for j = 1:tableWidth - 1
            routhMat(i,j) = (routhMat(i - 1, 1)*routhMat(i - 2, j + 1) - ...
                             routhMat(i - 2, 1)*routhMat(i - 1, j + 1))/...
                             routhMat(i - 1, 1);
        end
        zeroRow = false;
        zeroElm = false;

        if isempty(symvar(routhMat(i, :)))
            checkfun = @(item) all(abs(item) < opt.Tol);
        else
            checkfun = @(item) all(item == 0);
        end
        if checkfun(routhMat(i, 1))
            if checkfun(routhMat(i, 2:tableWidth)) && i <= deg
                zeroElm = false;
                zeroRow = true;
            else
                zeroElm = true;
                epsSubed = true;
            end
        end

        % If entire row = 0 ceq = characteristic equation
        if zeroRow
            auxDeg = deg - i + 2;
            charEq = zeros(1, 'sym');
            for j = auxDeg:-2:0
                charEq = charEq + routhMat(i - 1, (auxDeg - j)/2 + 1)*s^(j);
            end
            if opt.Verbos
                if rem(i, 2) == 0    
                    fprintf('Even Polynomial Found as:\n')
                else
                    fprintf('Odd Polynomial Found as:\n')
                end
                disp(charEq)
                fprintf('Roots:\n')
                disp(solve(charEq, s))
            end            
            diffAuxSym = diff(charEq, s); 
            descList{i, 1} = 'Derivative Substituted';
            diffAuxpoly = coeffs(diffAuxSym, s, 'all');
            routhMat(i, 1:auxDeg/2) = diffAuxpoly(1:2:end - 1);
        end

        %if only first number of row = 0
        if zeroElm
            routhMat(i, 1) = e(i);
            descList{i, 1} = 'Epsilon Substituted';
        end       
    end
    
    if opt.limit && epsSubed
        for i = 1:deg
            routhMat = limit(routhMat, e(i), 0, 'Right');
        end
    end
    if opt.simplify
        routhMat = simplify(routhMat);
    end
    columnNames = cell(1, tableWidth + 1);
    columnsData = cell(1, tableWidth + 1);
    for i = 1:tableWidth
        columnsData{i} = string(routhMat(:, i));
        columnNames{i} = ['Column ', num2str(i)];
    end
    columnNames{1, end} = 'Description';
    columnsData{1, end} = string(descList);
    
    routhTable = table(columnsData{:}, ...
                     'RowNames', powerList, 'VariableNames', columnNames);    
end