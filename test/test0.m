run('para.m');  % 載入基本參數

scenario = satelliteScenario;
tleFile = "eccentricOrbitSatellite.tle";
testSat = satellite(scenario, tleFile, 'Name', 'TestSat');
% testSat.Orbit.Keplerian = [a, e, i_deg, 0, omega_deg, 0];
% inclinationValue = sat.Orbit.Inclination;

% disp(inclinationValue);
