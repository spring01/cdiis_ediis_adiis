function [orbital, orbitalEnergies, totalEnergy, iter] = SCFDummy(obj)
oeiVec = reshape(obj.kineticMat + obj.corePotentialMat, [], 1);
teiForCoulomb = reshape(obj.twoElecIntegrals, length(oeiVec), []);
teiForExchange = reshape(permute(obj.twoElecIntegrals, [1 3 2 4]), ...
    length(oeiVec), []);
inv_S_Half = eye(size(obj.overlapMat)) / sqrtm(obj.overlapMat);

densVec = zeros(size(oeiVec));
elecEnergy = 0;
fockVec = oeiVec;
for iter = 1:obj.maxSCFIter
    oldDensVec = densVec;
    oldElecEnergy = elecEnergy;
    [densVec, elecEnergy, orbital, orbitalEnergies] ...
        = DiagonalizeFock(reshape(fockVec, sqrt(length(fockVec)), []), ...
        inv_S_Half, obj.numElectrons);
    
    if(norm(densVec - oldDensVec) < obj.RMSDensityThreshold && max(abs(densVec - oldDensVec)) < obj.MaxDensityThreshold ...
            && abs(elecEnergy - oldElecEnergy) < obj.EnergyThreshold)
        break;
    end
    
    fockVec = oeiVec + ... H
        2 .* (teiForCoulomb * densVec) ... +2J
        - (teiForExchange * densVec); ... -K
end
totalEnergy = oeiVec' * densVec + ...
    elecEnergy + obj.nuclearRepulsionEnergy;
end


function [densityVec, elecEnergy, orbital, orbitalEnergies] = DiagonalizeFock(fockMat, inv_S_Half, numElectrons)
[orbitalOtho, orbitalEnergies] = eig(inv_S_Half*fockMat*inv_S_Half);
[orbitalEnergies, ascend_order] = sort(diag(orbitalEnergies));
orbital = inv_S_Half * orbitalOtho(:, ascend_order);
densityVec = reshape( ...
    orbital(:, 1:numElectrons/2) * orbital(:, 1:numElectrons/2)', ...
    [], 1);
elecEnergy = sum(orbitalEnergies(1:numElectrons/2));
end

