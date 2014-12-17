classdef ADIIS < handle
    
    properties (Access = private)
        
        fockVectors;
        densVectors;
        
    end
    
    methods
        
        function obj = ADIIS(initVector, numVectors)
            if(nargin < 2)
                numVectors = 5;
            end
            obj.fockVectors = zeros(numel(initVector), numVectors);
            obj.densVectors = zeros(numel(initVector), numVectors);
        end
        
        function Push(obj, newFockVector, newDensVector)
            % push in matrices or rows or columns
            newFockVector = reshape(newFockVector, [], 1);
            newDensVector = reshape(newDensVector, [], 1);
            
            % push new Fock in
            obj.fockVectors(:, 1:end-1) = obj.fockVectors(:, 2:end);
            obj.fockVectors(:, end) = newFockVector;
            
            % push new density in
            obj.densVectors(:, 1:end-1) = obj.densVectors(:, 2:end);
            obj.densVectors(:, end) = newDensVector;
        end
        
        function newFockVector = Interpolate(obj)
            % compute errors wrt. the latest Fock or density
            errorFockVectors = obj.fockVectors - ...
                repmat(obj.fockVectors(:,end), 1, obj.NumVectors);
            errorDensVectors = obj.densVectors - ...
                repmat(obj.densVectors(:,end), 1, obj.NumVectors);
            
            % first order term and Hessian
            firstOrder = 2.*(obj.fockVectors(:, end)'*errorDensVectors);
            hessian = errorFockVectors'*errorDensVectors;
            hessian = hessian + hessian'; % multiply Hessian by 2 and cancels numerical error
            
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