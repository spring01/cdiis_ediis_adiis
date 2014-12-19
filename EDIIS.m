classdef EDIIS < handle
    
    properties (Access = private)
        
        oeiVector;
        
        fockVectors;
        densVectors;
        
    end
    
    methods
        
        function obj = EDIIS(oeiVector, numVectors)
            oeiVector = reshape(oeiVector, [], 1);
            if(nargin < 2)
                numVectors = 5;
            end
            obj.fockVectors = zeros(length(oeiVector), numVectors);
            obj.densVectors = zeros(length(oeiVector), numVectors);
            
            obj.oeiVector = oeiVector;
        end
        
        function Push(obj, newFockVector, newDensVector)
            % push in matrices or rows or columns
            newFockVector = reshape(newFockVector, [], 1);
            newDensVector = reshape(newDensVector, [], 1);
            
            % push new Fock and Fock error in
            obj.fockVectors(:, 1:end-1) = obj.fockVectors(:, 2:end);
            obj.fockVectors(:, end) = newFockVector;
            
            % push new density error in
            obj.densVectors(:, 1:end-1) = obj.densVectors(:, 2:end);
            obj.densVectors(:, end) = newDensVector;
        end
        
        function newFockVector = Interpolate(obj)
            firstOrder = zeros(obj.NumVectors(), 1);
            for i = 1:obj.NumVectors()
                firstOrder(i) = (obj.oeiVector + obj.fockVectors(:, i))' ...
                    * obj.densVectors(:, i);
            end
            
            hessian = zeros(obj.NumVectors());
            for i = 1:obj.NumVectors()
                for j = 1:obj.NumVectors();
                    hessian(j, i) ...
                        = (obj.fockVectors(:, i) - obj.fockVectors(:, j))' ...
                        * (obj.densVectors(:, i) - obj.densVectors(:, j));
                end
            end
            hessian = -hessian;
            
            % reduced gradient
            coeffs = obj.ReducedGradient( ...
                hessian, firstOrder, [zeros(obj.NumVectors()-1,1); 1]);
            
            newFockVector = obj.fockVectors * coeffs;
        end
        
    end
    
    methods (Access = private)
        
        function num = NumVectors(obj)
            num = size(obj.fockVectors, 2);
        end
        
        function varFull = ReducedGradient(obj, hessian, firstOrder, iniPoint)
            indDep = 1;
            indAct = 1:obj.NumVectors();
            indAct = indAct(indAct~=indDep);
            varFull = iniPoint;
            constr = ones(1, obj.NumVectors());
            for iter = 1:2000
                varDep = varFull(indDep);
                varAct = varFull(indAct);
                constrDep = constr(indDep);
                constrAct = constr(indAct);
                
                grad = firstOrder + hessian*varFull;
                gradDep = grad(indDep);
                gradAct = grad(indAct);
                
                redGrad = gradAct - constrAct'/constrDep*gradDep;
                
                delVarAct = zeros(obj.NumVectors()-1,1);
                delVarAct(redGrad<0) = -redGrad(redGrad<0);
                delVarAct(varAct>0) = -redGrad(varAct>0);
                
                if(norm(delVarAct) < 1e-8)
                    break;
                end
                
                stepSize = 2 ./ (2 + iter);
                varActSim = varAct + stepSize .* delVarAct;
                delVarAct(varActSim<0) = 0;
                
                delVarDep = - constrDep \ constrAct * delVarAct;
                varDepSim = varDep + stepSize .* delVarDep;
                
                if(varDepSim < 0)
                    strPos = indAct(varActSim>0);
                    indDep = strPos(1);
                    indAct = 1:obj.NumVectors();
                    indAct = indAct(indAct~=indDep);
                    continue;
                end
                
                delVarFull = zeros(obj.NumVectors(), 1);
                delVarFull(indDep) = delVarDep;
                delVarFull(indAct) = delVarAct;
                
                varFull = varFull + stepSize .* delVarFull;
                
            end
        end
        
    end
    
end