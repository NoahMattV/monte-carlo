% turn energy value into Eint

function out = energyToInt(energy)
    global e;
    E = linspace(1e-255,2*e,1001); % 0 to 2*e Joules

    A = repmat(energy, [1 length(E)]);
    [minValue,EInt] = min(abs(A-E));
    out = EInt;
end