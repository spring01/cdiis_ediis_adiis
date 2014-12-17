classdef RHF < handle
    
    properties (SetAccess = private)
        
        overlapMat;
        kineticMat;
        corePotentialMat;
        twoElecIntegrals;
        
        nuclearRepulsionEnergy;
        numElectrons;
        
    end
    
    properties (Access = private)
        
        maxSCFIter = 500;
        RMSDensityThreshold = 1e-8;
        MaxDensityThreshold = 1e-6;
        EnergyThreshold = 1e-6;
        
    end
    
    methods
        
        function obj = RHF(properties)
            obj.overlapMat = properties.overlapMat;
            obj.kineticMat = properties.kineticMat;
            obj.corePotentialMat = properties.corePotentialMat;
            obj.twoElecIntegrals = properties.twoElecIntegrals;
            obj.nuclearRepulsionEnergy = properties.nuclearRepulsionEnergy;
            obj.numElectrons = properties.numElectrons;
        end
        
    end
    
    methods (Static)
                
        function properties = MatPsiInterface(matpsi)
            properties.overlapMat = matpsi.overlap();
            properties.kineticMat = matpsi.kinetic();
            properties.corePotentialMat = matpsi.potential();
            properties.twoElecIntegrals = matpsi.tei_allfull();
            properties.nuclearRepulsionEnergy = matpsi.Enuc();
            properties.numElectrons = matpsi.nelec();
        end
    end
    
end