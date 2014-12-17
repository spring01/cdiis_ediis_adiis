classdef CDIIS < handle
    
    properties (Access = private)
        
        fockVectors;
        errorCommutatorVectors;
        
        S_Half;
        
        startError = 1e-2;
        
    end
    
    methods
        
        function obj = CDIIS(overlapMatrix, numVectors)
            if(nargin < 2)
                numVectors = 5;
            end
            obj.fockVectors = zeros(numel(overlapMatrix), numVectors);
            obj.errorCommutatorVectors = zeros(numel(overlapMatrix), numVectors);
            
            obj.S_Half = sqrtm(overlapMatrix);
        end
        
        function Push(obj, newFockVector, newDensVector)
            % push in matrices or rows or columns
            newFockVector = reshape(newFockVector, [], 1);
            newDensVector = reshape(newDensVector, [], 1);
            
            % push new Fock in
            obj.fockVectors(:, 1:end-1) = obj.fockVectors(:, 2:end);
            obj.fockVectors(:, end) = newFockVector;
            
            % push new commutator error in
            obj.errorCommutatorVectors(:, 1:end-1) = obj.errorCommutatorVectors(:, 2:end);
            
            FtDt = obj.S_Half ...
                \ reshape(newFockVector, sqrt(length(newFockVector)), []) ...
                * reshape(newDensVector, sqrt(length(newDensVector)), []) ...
                * obj.S_Half;
            obj.errorCommutatorVectors(:, end) = ...
                reshape(FtDt - FtDt', [], 1);
        end
        
        function better = IAmBetter(obj) % than EDIIS
            better = 0;
            if(sum(abs(obj.fockVectors(:,1))) ...
                    && obj.startError > max(abs(obj.errorCommutatorVectors(:, end))))
                better = 1;
            end
        end
        
        function newFockVector = Extrapolate(obj)
            % warning('off', 'MATLAB:nearlySingularMatrix');
            
            onesVec = ones(obj.NumVectors(),1);
            diisCoefficients = ...
                [obj.errorCommutatorVectors'*obj.errorCommutatorVectors, onesVec; ...
                onesVec', 0] ...
                \ [zeros(obj.NumVectors(),1); 1];
            newFockVector = obj.fockVectors * diisCoefficients(1:end-1);
        end
        
    end
    
    methods (Access = private)
        
        function num = NumVectors(obj)
            num = size(obj.fockVectors, 2);
        end
        
    end
    
end