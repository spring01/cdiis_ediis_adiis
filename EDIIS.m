classdef EDIIS < handle
    
    properties (Access = private)
        
        oeiVector;
        
        fockVectors;
        densVectors;
        
    end
    
    methods
        
        function obj = EDIIS(initFockVector, initDensVector, oeiVector, numVectors)
            if(nargin < 3)
                oeiVector = initFockVector;
            end
            if(nargin < 4)
                numVectors = 5;
            end
            obj.fockVectors = zeros(length(initFockVector), numVectors);
            obj.densVectors = zeros(length(initFockVector), numVectors);
            
            obj.fockVectors(:, end) = initFockVector;
            obj.densVectors(:, end) = initDensVector;
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
                    hessian(i,j) ...
                        = (obj.fockVectors(:, i) - obj.fockVectors(:, j))' ...
                        * (obj.densVectors(:, i) - obj.densVectors(:, j));
                end
            end
            hessian = -hessian;
            
            % matlab's QP solver (seems it's rgd and slow)
            coeffs = [zeros(obj.NumVectors()-1,1); 1];
            QPoptions = optimoptions('quadprog', 'Algorithm', 'active-set', 'Display', 'off');
            coeffs = quadprog(hessian, firstOrder, ...
                [], [], ...
                ones(1,obj.NumVectors()), 1, ...
                zeros(obj.NumVectors(),1), [], ...
                coeffs, QPoptions);
            
            newFockVector = obj.fockVectors * coeffs;
        end
        
    end
    
    methods (Access = private)
        
        function num = NumVectors(obj)
            num = size(obj.fockVectors, 2);
        end
        
    end
    
end