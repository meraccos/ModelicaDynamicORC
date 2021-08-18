package TotalModelTrial
  model ConstantPower_static
    import SI = Modelica.SIunits;
    import Modelica.Constants.pi;
    package Medium = Modelica.Media.Water.StandardWater;
    parameter Integer n = 50;
    parameter SI.Length L = 20;
    parameter SI.Length r = 0.05;
    final parameter SI.Length dx = L / n;
    final parameter SI.Volume dV = pi * r * r * dx;
    parameter SI.Pressure p0 = 1e5;
    parameter Medium.SpecificEnthalpy h0 = 250e3;
    parameter SI.Power Q = 35000;
    Medium.MassFlowRate w;
    Medium.ThermodynamicState state[n];
    Medium.SpecificEnthalpy h[n];
    Medium.Density rho[n];
    Real x[n];
  equation
    for i in 1:n - 1 loop
      w * h[i] + Q - w * h[i + 1] = 0;
    end for;
    w = Modelica.Media.Common.smoothStep(time - 5, 0.5, 0.3, 3);
    h[1] = h0;
    state = Medium.setState_phX(p0, h);
    x = Medium.vapourQuality(state);
    rho = Medium.density(state);
  end ConstantPower_static;

  model ConstantPower_dynamic
    import SI = Modelica.SIunits;
    import Modelica.Constants.pi;
    package Medium = Modelica.Media.Water.StandardWater;
    parameter Integer n = 50;
    parameter SI.Length L = 20;
    parameter SI.Length r = 0.05;
    final parameter SI.Length dx = L / n;
    final parameter SI.Volume dV = pi * r * r * dx;
    parameter SI.Pressure p0 = 1e5;
    parameter Medium.SpecificEnthalpy h0 = 250e3;
    parameter SI.Power Q = 35000;
    Medium.MassFlowRate w;
    Medium.ThermodynamicState state[n];
    Medium.SpecificEnthalpy h[n];
    Medium.Density rho[n];
    Real x[n];
  initial equation
    for i in 1:n - 1 loop
      w * h[i] + Q - w * h[i + 1] = 0;
    end for;
  equation
    for i in 1:n - 1 loop
      rho[i + 1] * dV * der(h[i + 1]) = Q + w * (h[i] - h[i + 1]);
    end for;
    w = Modelica.Media.Common.smoothStep(time - 5, 0.8, 0.5, 3);
    h[1] = h0;
    state = Medium.setState_phX(p0, h);
    x = Medium.vapourQuality(state);
    rho = Medium.density(state);
  end ConstantPower_dynamic;

  model ConstantPower_compression
    import SI = Modelica.SIunits;
    import Modelica.Constants.pi;
    package Medium = Modelica.Media.Water.StandardWater;
    parameter Integer n = 30;
    parameter SI.Length L = 20;
    parameter SI.Length r = 0.05;
    final parameter SI.Length dx = L / n;
    final parameter SI.Volume dV = pi * r * r * dx;
    parameter Medium.SpecificEnthalpy h0 = 250e3;
    parameter SI.Power Q = 1000;
    Medium.MassFlowRate w_0;
    Medium.MassFlowRate w[n](start = fill(0.5, n));
    parameter SI.Pressure p_0 = 1e5;
    Medium.ThermodynamicState state[n];
    Medium.SpecificEnthalpy h[n](start = fill(5.5*h0, n));
    Medium.Density rho[n];
    Real x[n];
  initial equation
    for i in 1:n - 1 loop
      Q + w_0 * h[i] = w_0 * h[i + 1];
    end for;
    for i in 2:n loop
      w[i] = w_0;
      der(rho[i]) = 0;
    end for;
  equation
    for i in 1:n - 1 loop
      rho[i + 1] * dV * der(h[i + 1]) = homotopy(Q + w[i] * h[i] - w[i + 1] * h[i + 1], 0);
    end for;
    dV * der(rho[2:n]) = w[1:n-1] - w[2:n];  
    w[1] = w_0;
    h[1] = h0;
    w_0 = homotopy(Modelica.Media.Common.smoothStep(time - 5, 0.8, 0.5, 3), w_0);
    state[1:n] = Medium.setState_phX(p_0, h[1:n]);
    x = Medium.vapourQuality(state);
    rho = Medium.density(state);
  end ConstantPower_compression;
end TotalModelTrial;