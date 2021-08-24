package Rankine
  import SI = Modelica.SIunits;
  import Modelica.Constants.pi;
  package Tests
    model Test_StaticHX
      ThermoPower.Water.SourceMassFlow sourceMassFlow1(T = 365, h = 2.8e6, use_T = false, w0 = 5) annotation(
        Placement(visible = true, transformation(origin = {54, -16}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
      ThermoPower.Water.SourceMassFlow sourceMassFlow(h = 250e3, w0 = 5) annotation(
        Placement(visible = true, transformation(origin = {-70, 26}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      ThermoPower.Water.SinkPressure sinkPressure1(h = 0) annotation(
        Placement(visible = true, transformation(origin = {-68, -14}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
      inner ThermoPower.System system annotation(
        Placement(visible = true, transformation(origin = {84, 90}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      ThermoPower.Water.SinkPressure sinkPressure(h = 0) annotation(
        Placement(visible = true, transformation(origin = {46, 26}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Components.StaticHX staticHX(h_t = 8000, hc_in_nom = 250e3, hh_in_nom = 2.8e6, n = 50, p_cold(displayUnit = "Pa"), p_cold_nom(displayUnit = "Pa") = 1e5, p_hot(displayUnit = "Pa"), p_hot_nom(displayUnit = "Pa") = 1e5, wc_nom = 5, wh_nom = 5)  annotation(
        Placement(visible = true, transformation(origin = {-10, 4}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    equation
  connect(sourceMassFlow.flange, staticHX.cold_inlet) annotation(
        Line(points = {{-60, 26}, {-16, 26}, {-16, 12}}, color = {0, 0, 255}));
  connect(staticHX.cold_outlet, sinkPressure.flange) annotation(
        Line(points = {{-4, 12}, {-2, 12}, {-2, 26}, {36, 26}}, color = {0, 0, 255}));
  connect(sinkPressure1.flange, staticHX.hot_outlet) annotation(
        Line(points = {{-58, -14}, {-16, -14}, {-16, -6}}, color = {0, 0, 255}));
  connect(sourceMassFlow1.flange, staticHX.hot_inlet) annotation(
        Line(points = {{44, -16}, {-4, -16}, {-4, -6}}, color = {0, 0, 255}));
    end Test_StaticHX;

    model Test_QuasiStaticHX
      ThermoPower.Water.SinkPressure sinkPressure1 annotation(
        Placement(visible = true, transformation(origin = {-68, -14}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
      ThermoPower.Water.SourceMassFlow sourceMassFlow1( h = 2.8e6, use_T = false, w0 = 5) annotation(
        Placement(visible = true, transformation(origin = {54, -16}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
      ThermoPower.Water.SourceMassFlow sourceMassFlow(h = 250e3, use_in_w0 = true, w0 = 0.5) annotation(
        Placement(visible = true, transformation(origin = {-70, 26}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      inner ThermoPower.System system annotation(
        Placement(visible = true, transformation(origin = {84, 90}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      ThermoPower.Water.SinkPressure sinkPressure annotation(
        Placement(visible = true, transformation(origin = {46, 26}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.Ramp ramp(duration = 2, height = 1, offset = 0.5, startTime = 10) annotation(
        Placement(visible = true, transformation(origin = {-112, 54}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Components.StaticHX staticHX( h_t = 8000, hc_in_nom = 250e3, hh_in_nom = 2.8e6, n = 25, p_cold(displayUnit = "Pa"), p_cold_nom(displayUnit = "Pa") = 1e5, p_hot(displayUnit = "Pa"), p_hot_nom(displayUnit = "Pa") = 1e5, wc_nom = 5, wh_nom = 5) annotation(
        Placement(visible = true, transformation(origin = {-10, 4}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    equation
      connect(ramp.y, sourceMassFlow.in_w0) annotation(
        Line(points = {{-100, 54}, {-74, 54}, {-74, 32}}, color = {0, 0, 127}));
  connect(sourceMassFlow.flange, staticHX.cold_inlet) annotation(
        Line(points = {{-60, 26}, {-40, 26}, {-40, 12}, {-16, 12}}, color = {0, 0, 255}));
  connect(staticHX.cold_outlet, sinkPressure.flange) annotation(
        Line(points = {{-4, 12}, {16, 12}, {16, 26}, {36, 26}}, color = {0, 0, 255}));
  connect(sinkPressure1.flange, staticHX.hot_outlet) annotation(
        Line(points = {{-58, -14}, {-40, -14}, {-40, -6}, {-16, -6}}, color = {0, 0, 255}));
  connect(staticHX.hot_inlet, sourceMassFlow1.flange) annotation(
        Line(points = {{-4, -6}, {16, -6}, {16, -16}, {44, -16}}, color = {0, 0, 255}));
    end Test_QuasiStaticHX;

    model Test_DynamicHX
  Modelica.Blocks.Sources.Ramp ramp(duration = 2, height = 0.2, offset = 0.5, startTime = 10) annotation(
        Placement(visible = true, transformation(origin = {-112, 54}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      ThermoPower.Water.SinkPressure sinkPressure1 annotation(
        Placement(visible = true, transformation(origin = {-68, -14}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
      ThermoPower.Water.SourceMassFlow sourceMassFlow(h = 150e3, use_in_w0 = true, w0 = 0.5) annotation(
        Placement(visible = true, transformation(origin = {-70, 26}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      ThermoPower.Water.SinkPressure sinkPressure annotation(
        Placement(visible = true, transformation(origin = {46, 26}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      ThermoPower.Water.SourceMassFlow sourceMassFlow1(h = 2.8e6, use_T = false, w0 = 0.5) annotation(
        Placement(visible = true, transformation(origin = {54, -16}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
      inner ThermoPower.System system annotation(
        Placement(visible = true, transformation(origin = {86, 86}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Components.DynamicHX dynamicHX(Q_nom = 1, h_t = 2000, hc_in_nom = 150e3, hh_in_nom = 2.8e6, n = 30, p_cold(displayUnit = "Pa"), p_cold_nom(displayUnit = "Pa") = 99999.99999999999, p_hot(displayUnit = "Pa"), p_hot_nom(displayUnit = "Pa") = 99999.99999999999, use_q_nom = true, wc_nom = 0.5, wh_nom = 0.5) annotation(
        Placement(visible = true, transformation(origin = {-10, 4}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    equation
      connect(ramp.y, sourceMassFlow.in_w0) annotation(
        Line(points = {{-100, 54}, {-74, 54}, {-74, 32}}, color = {0, 0, 127}));
      connect(sourceMassFlow.flange, dynamicHX.cold_inlet) annotation(
        Line(points = {{-60, 26}, {-42, 26}, {-42, 10}, {-20, 10}}, color = {0, 0, 255}));
      connect(sinkPressure1.flange, dynamicHX.hot_outlet) annotation(
        Line(points = {{-58, -14}, {-42, -14}, {-42, -2}, {-20, -2}}, color = {0, 0, 255}));
      connect(sinkPressure.flange, dynamicHX.cold_outlet) annotation(
        Line(points = {{36, 26}, {16, 26}, {16, 10}, {0, 10}}, color = {0, 0, 255}));
      connect(sourceMassFlow1.flange, dynamicHX.hot_inlet) annotation(
        Line(points = {{44, -16}, {18, -16}, {18, -2}, {0, -2}}, color = {0, 0, 255}));
    end Test_DynamicHX;
  end Tests;

  package Models
    model StaticHX_model
      package Medium = Modelica.Media.Water.StandardWater;
      parameter Integer n = 20 "Number of finite divisions";
      parameter Integer N = 100 "Number of tubes";
      parameter SI.Length L = 2 "Length of the HX";
      parameter SI.Length r = 0.05 "Radius of the tubes";
      final parameter SI.Length dx = L / n;
      parameter SI.CoefficientOfHeatTransfer h_t = 150 "Heat transfer coefficient";
      SI.Temperature Th[n] "Hot side temperature";
      SI.Temperature Tc[n] "Cold side temperature";
      Medium.SpecificEnthalpy hh[n](start = linspace(0.73 * hh_in, hh_in, n)) "Hot side specific enthalpy";
      Medium.SpecificEnthalpy hc[n](start = linspace(hc_in, 8.8 * hc_in, n)) "Cold side specific enthalpy";
      Medium.ThermodynamicState state_h[n];
      Medium.ThermodynamicState state_c[n];
      SI.Power Qt[n] "Heat transfer at slice i";
      Real x_h[n] "Hot side vapour quality";
      Real x_c[n] "Hot side vapour quality";
      SI.Temperature dT[n] "Temperature difference";
      //The following variable and parameters can be substituted with connectors
      parameter Medium.MassFlowRate wh = wc / N "Hot side mass flow rate";
      parameter Medium.MassFlowRate wc = 0.5 "Cold side mass flow rate";
      parameter Medium.SpecificEnthalpy hh_in = 3004772 "Hot side inlet enthalpy";
      parameter Medium.SpecificEnthalpy hc_in = 150000 "Cold side inlet enthalpy";
      parameter SI.Pressure p_hot = 1e5 "Hot side Pressure";
      parameter SI.Pressure p_cold = 1e5 "Cold side pressure";
    protected
      parameter Medium.Temperature Tc_in = Medium.temperature_ph(p_cold, hc_in);
      parameter Medium.Temperature Th_in = Medium.temperature_ph(p_hot, hh_in);  
    equation
    //Heat transfer Q_i ( the second value in homotopy is an approximation for dT_ave )
      Qt[1:n] = 2 * pi * r * dx * h_t * homotopy(Th[1:n] - Tc[1:n], (Th[n] - Tc[1]) / 2);
    //Energy conservation for cold fluid
      Qt[2:n] * N + wc * hc[1:n - 1] = wc * hc[2:n];
    //Energy conservation for hot fluid
      Qt[1:n - 1] + wh * hh[1:n - 1] = wh * hh[2:n];
    //Setting the states
      state_h = Medium.setState_phX(p_hot, hh);
      state_c = Medium.setState_phX(p_cold, hc);
    //Extracting the temperatures from the states
      Th = state_h.T;
      Tc = state_c.T;
    //Extracting the qualities from the states
      x_h = Medium.vapourQuality(state_h);
      x_c = Medium.vapourQuality(state_c);
    //Boundary conditions
      hc[1] = hc_in;
      hh[n] = hh_in;
    //Temperature difference at the inlet and outlet
      dT = Th - Tc;
      annotation(
        __OpenModelica_commandLineOptions = "--matchingAlgorithm=PFPlusExt --indexReductionMethod=dynamicStateSelection -d=initialization,NLSanalyticJacobian -d=bltdump",
        __OpenModelica_simulationFlags(lv = "LOG_STATS", s = "dassl"),
        Documentation(info = "<html><head></head><body>Static<div>Shell and tube</div><div>Constant HT coefficient</div><div>Incompressible</div><div>All phases</div><div><br></div><div>7/27/2021</div><div>Meraj Mammadov</div></body></html>"));
    end StaticHX_model;

    //Setting the states
    //state_h[1:n] = Medium.setState_phX(p_hot, hh[1:n]);
    ////state_h[1] = Medium.setState_pT(p_hot, Th_in);
    //state_c[1:n] = Medium.setState_phX(p_cold, hc[1:n]);

    model DynamicHX_model
      package Medium = Modelica.Media.Water.StandardWater;
      parameter Integer n = 50 "Number of finite divisions";
      parameter Integer N = 100 "Number of tubes";
      parameter SI.Length L = 2 "Length of the HX";
      parameter SI.Length r = 0.05 "Radius of the tubes";
      final parameter SI.Length dx = L / n;
      final parameter SI.Volume dV = pi * r * r * dx;
      parameter SI.CoefficientOfHeatTransfer h_t = 150 "Heat transfer coefficient";
      SI.Temperature Th[n] "Hot side temperature";
      SI.Temperature Tc[n] "Cold side temperature";
      Medium.SpecificEnthalpy hh[n] "Hot side specific enthalpy";
      Medium.SpecificEnthalpy hc[n] "Cold side specific enthalpy";
      Medium.ThermodynamicState state_h[n];
      Medium.ThermodynamicState state_c[n];
      SI.Power Qt[n] "Heat transfer at slice i";
      Real x_h[n] "Hot side vapour quality";
      Real x_c[n] "Hot side vapour quality";
      SI.Temperature dT[n] "Temperature difference";
      //The following variable and parameters can be substituted with connectors
      parameter Medium.MassFlowRate wh = 0.01 "Hot side mass flow rate";
      Medium.MassFlowRate wc "Cold side mass flow rate";
      parameter Medium.SpecificEnthalpy hh_in = 2804772 "Hot side inlet enthalpy";
      parameter Medium.SpecificEnthalpy hc_in = 250000 "Cold side inlet enthalpy";
      parameter SI.Pressure p_hot = 1e5 "Hot side Pressure";
      parameter SI.Pressure p_cold = 1e5 "Cold side pressure";
    initial equation
      Qt[2:n] * N = wc * (hc[2:n] - hc[1:n - 1]);
      Qt[1:n - 1] = wh * (hh[2:n] - hh[1:n - 1]);
    equation
//Dynamic change is given here, can be changed upon wish
      wc = Modelica.Media.Common.smoothStep(time - 50, 1.5, 0.5, 1);
//Heat transfer Q_i ( the second value in homotopy is an approximation for dT_ave )
      Qt[1:n] = 2 * pi * r * dx * h_t * homotopy(Th[1:n] - Tc[1:n], (Th[n] - Tc[1]) / 2);
//Energy conservation for cold fluid
      for i in 2:n loop
        Qt[i] * N + wc * (hc[i - 1] - hc[i]) = state_c[i].d * N * dV * der(hc[i]);
      end for;
//Energy conservation for hot fluid
      for i in 1:n - 1 loop
        (-Qt[i]) + wh * (hh[i + 1] - hh[i]) = state_h[i].d * dV * der(hh[i]);
      end for;
//Setting the states
      state_h = Medium.setState_phX(p_hot, hh);
      state_c = Medium.setState_phX(p_cold, hc);
//Extracting the temperatures from the states
      Th = state_h.T;
      Tc = state_c.T;
//Extracting the qualities from the states
      x_h = Medium.vapourQuality(state_h);
      x_c = Medium.vapourQuality(state_c);
//Boundary conditions
      hc[1] = hc_in;
      hh[n] = hh_in;
    //Temperature difference at the inlet and outlet
      dT = Th - Tc;
      annotation(
        __OpenModelica_commandLineOptions = "--matchingAlgorithm=PFPlusExt --indexReductionMethod=dynamicStateSelection -d=initialization,NLSanalyticJacobian -d=bltdump",
        __OpenModelica_simulationFlags(lv = "LOG_STATS", s = "dassl"));
    end DynamicHX_model;

    partial model BaseHX
      package Medium = Modelica.Media.Water.StandardWater;
      
      parameter Integer n = 25 "Number of finite divisions";
      parameter Integer N = 30 "Total number of plates";
      
      parameter SI.Length L = 0.8 "Length of the HX";
      parameter SI.Length W = 0.3 "Width of the HX";  
      parameter SI.Length t = 0.2 / N "Width of the HX";   
      final parameter SI.Length dx = L / n;
      final parameter SI.Volume dV = W * dx * t;  
      
      parameter SI.CoefficientOfHeatTransfer h_t = 400 "Heat transfer coefficient";
      
      SI.Temperature Th[n] (start = linspace(Tc_in, Th_in, n), each displayUnit = "K") "Hot side temperature";
      SI.Temperature Tc[n] (start = linspace(Tc_in, Th_in, n), each displayUnit = "K") "Cold side temperature";
      SI.Temperature dT[n] (start = fill(5, n), each displayUnit = "K")"Temperature difference";  
      
      Medium.SpecificEnthalpy hh[n] (start = linspace(hc_in_nom, hh_in_nom, n))"Hot side specific enthalpy";
      Medium.SpecificEnthalpy hc[n] (start = linspace(hc_in_nom, hh_in_nom, n))"Cold side specific enthalpy";
      
      Medium.ThermodynamicState state_h[n] annotation(HideResult = false );
      Medium.ThermodynamicState state_c[n];
      
      SI.Power Qt[n] "Heat transfer at slice i";
      parameter SI.Power Q_nom = (Th_in - Tc_in) / 2
      annotation (Dialog(tab="Initialisation"));   
      
      Real x_h[n] "Hot side vapour quality";
      Real x_c[n] "Hot side vapour quality";
      
      Medium.MassFlowRate wh"Hot side mass flow rate";
      Medium.MassFlowRate wc"Cold side mass flow rate";
      
      parameter Medium.MassFlowRate wh_nom = wc_nom "Hot side  nominal mass flow rate"
      annotation (Dialog(tab="Initialisation"));
      parameter Medium.MassFlowRate wc_nom = 0.5 "Cold side nominal mass flow rate"
      annotation (Dialog(tab="Initialisation"));  
      
      Medium.SpecificEnthalpy hh_in "Hot side inlet enthalpy";
      Medium.SpecificEnthalpy hc_in "Cold side inlet enthalpy";
    
      parameter Medium.SpecificEnthalpy hh_in_nom "Hot side inlet enthalpy"
      annotation (Dialog(tab="Initialisation"));
      parameter Medium.SpecificEnthalpy hc_in_nom "Cold side inlet enthalpy"
      annotation (Dialog(tab="Initialisation"));  
      
      SI.Pressure p_hot(start = 1e5) "Hot side Pressure";
      SI.Pressure p_cold(start = 1e5) "Cold side pressure";
      
      parameter SI.Pressure p_hot_nom = 1e5 "Hot side Pressure"
      annotation (Dialog(tab="Initialisation"));
      parameter SI.Pressure p_cold_nom = 1e5 "Cold side pressure"
      annotation (Dialog(tab="Initialisation"));    
      
      
      parameter SI.Temperature Tc_in = Medium.temperature_ph(p_cold_nom, hc_in_nom)
      annotation (Dialog(tab="Initialisation"));
      parameter SI.Temperature Th_in = Medium.temperature_ph(p_hot_nom, hh_in_nom)
      annotation (Dialog(tab="Initialisation")); 
         
    equation
    //Thermodynamic states
      state_h = Medium.setState_phX(p_hot, hh);
      state_c = Medium.setState_phX(p_cold, hc);
      
    //Temperatures
      Th = state_h.T;
      Tc = state_c.T;
      
    //Vapor qualities
      x_h = Medium.vapourQuality(state_h);
      x_c = Medium.vapourQuality(state_c);
      
    //Boundary conditions
      hc[1] = hc_in;
      hh[n] = hh_in;
      
    //Temperature difference
      dT = Th - Tc;

    annotation(
        Icon(graphics = {Bitmap(origin = {71, -1}, extent = {{-239, 101}, {95, -101}}, imageSource = "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AACEYBCXDOIAAAAAEIJBXDKIAwAAAEAIBnHJIA4AAAAAIRjEJYM4AAAAAIRgEJcM4gAAAAAQgkFcMogDAAAAQAgGcckgDgAAAAAhGMQlgzgAAAAAhGAQlwziAAAAABCCQVwyiAMAAABACAZxySAOAAAAACEYxCWDOAAAAACEYBCXDOIAAAAAEIJBXDKIAwAAAEAIBnHJIA4AAAAAIRjEJYM4AAAAAIRgEJcM4gAAAAAQgkFcMogDAAAAQAgGcckgDgAAAAAhGMQlgzgAAAAAhGAQlwziAAAAABCCQVwyiAMAAABACAZxySAOAAAAACEYxCWDOAAAAACEYBCXDOIAAAAAEIJBXDKIAwAAAEAIBnHJIA4AAAAAIRjEJYM4AAAAAIRgEJcM4gAAAAAQgkFcMogDAAAAQAgGcckgDgAAAAAhGMQlgzgAAAAAhGAQlwziAAAAABCCQVwyiAMAAABACAZxySAOAAAAACEYxCWDOAAAAACEYBCXDOIAAAAAEIJBXDKIAwAAAEAIBnHJIA4AAAAAIRjEJYM4AAAAAIRgEJcM4gAAAAAQgkFcMogDAAAAQAgGcckgDgAAAAAhGMQlgzgAAAAAhGAQlwziAAAAABCCQVwyiAMAAABACAZxySAOAAAAACEYxCWDOAAAAACEYBCXDOIAAAAAEIJBXDKIAwAAAEAIBnHJIA4AAAAAIRjEJYM4AAAAAIRgEJcM4gAAAAAQgkFcMogDAAAAQAgGcckgDgAAAAAhGMQlgzgAAAAAhGAQlwziAAAAABCCQVwyiAMAAABACAZxySAOAAAAACEYxCWDOAAAAACEYBCXDOIAAAAAEIJBXDKIAwAAAEAIBnHJIA4AAAAAIRjEJYM4AAAAAIRgEJcM4gAAAAAQgkFcMogDAAAAQAgGcckgDgAAAAAhGMQlgzgAAAAAhGAQlwziAAAAABCCQVwyiAMAAABACAZxySAOAAAAACEYxCWDOAAAAACEYBCXDOIAAAAAEIJBXDKIAwAAAEAIBnHJIA4AAAAAIRjEJYM4AAAAAIRgEJcM4gAAAAAQgkFcMogDAAAAQAgGcckgDgAAAAAhGMQlgzgAAAAAhGAQlwziAAAAABCCQVwyiAMAAABACAZxySAOAAAAACEYxCWDOAAAAACEYBCXDOIAAAAAEIJBXDKIAwAAAEAIBnHJIA4AAAAAIRjEJYM4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"), Rectangle(origin = {-65, 61}, fillPattern = FillPattern.Solid, extent = {{-15, 9}, {15, -9}}), Rectangle(origin = {57, -66}, fillPattern = FillPattern.Solid, extent = {{-15, 12}, {15, -12}}), Text(origin = {-2, 66}, extent = {{-42, 12}, {42, -12}}, textString = "cold"), Text(origin = {-1, -68}, extent = {{-39, 12}, {39, -12}}, textString = "hot")}));
    end BaseHX;
  end Models;

  package Components
    model StaticHX
      extends Models.BaseHX;
    
      ThermoPower.Water.FlangeA cold_inlet annotation(
        Placement(visible = true, transformation(origin = {-100, 64}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-66, 88}, extent = {{-18, -18}, {18, 18}}, rotation = 0)));
      ThermoPower.Water.FlangeB cold_outlet annotation(
        Placement(visible = true, transformation(origin = {100, 60}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {62, 86}, extent = {{-18, -18}, {18, 18}}, rotation = 0)));
      ThermoPower.Water.FlangeA hot_inlet annotation(
        Placement(visible = true, transformation(origin = {100, -60}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {59, -97}, extent = {{-19, -19}, {19, 19}}, rotation = 0)));
      ThermoPower.Water.FlangeB hot_outlet annotation(
        Placement(visible = true, transformation(origin = {-98, -60}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-61, -99}, extent = {{-19, -19}, {19, 19}}, rotation = 0)));
    equation
    //Heat transfer
      Qt[1:n] = W * dx * h_t * homotopy(dT[1:n], (Th[n] - Tc[1]) / 2);
      
    //Energy conservation
      Qt[2:n] * (N-1)= wc * (hc[2:n] - hc[1:n - 1]) "cold fluid";
      Qt[1:n - 1] * (N-1)= wh * (hh[2:n] - hh[1:n - 1]) "hot fluid";
    
    //Flange equations
    
    //Mass flow rates
      cold_outlet.m_flow + cold_inlet.m_flow = 0;
      wc = homotopy(cold_inlet.m_flow, wc_nom);
      hot_inlet.m_flow + hot_outlet.m_flow = 0;
      wh = homotopy(hot_inlet.m_flow, wh_nom);
      
    //Enthalpy at the flanges
      hc_in = inStream(cold_inlet.h_outflow);
      hh_in = inStream(hot_inlet.h_outflow);
      cold_inlet.h_outflow = hc[1];
      hot_inlet.h_outflow = hh[n];
      cold_outlet.h_outflow = hc[n];
      hot_outlet.h_outflow = hh[1];
      
    //Pressure at the flanges
      cold_inlet.p = cold_outlet.p;
      p_cold = cold_inlet.p;
      hot_inlet.p = hot_outlet.p;
      p_hot = hot_inlet.p;
      annotation(
        __OpenModelica_commandLineOptions = "--matchingAlgorithm=PFPlusExt --indexReductionMethod=dynamicStateSelection -d=initialization,NLSanalyticJacobian -d=bltdump",
        __OpenModelica_simulationFlags(lv = "LOG_STATS", s = "dassl"),
        Icon,
        Documentation(info = "<html><head></head><body>Static<div>Shell and tube</div><div>Constant HT coefficient</div><div>Incompressible</div><div>All phases</div><div><br></div><div>7/27/2021</div><div>Meraj Mammadov</div></body></html>"));
    end StaticHX;

    model DynamicHX
      extends Models.BaseHX;
    
      ThermoPower.Water.FlangeA cold_inlet annotation(
        Placement(visible = true, transformation(origin = {-100, 64}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-66, 88}, extent = {{-18, -18}, {18, 18}}, rotation = 0)));
      ThermoPower.Water.FlangeB cold_outlet annotation(
        Placement(visible = true, transformation(origin = {100, 60}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {61, 87}, extent = {{-19, -19}, {19, 19}}, rotation = 0)));
      ThermoPower.Water.FlangeA hot_inlet annotation(
        Placement(visible = true, transformation(origin = {100, -60}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {58, -96}, extent = {{-18, -18}, {18, 18}}, rotation = 0)));
      ThermoPower.Water.FlangeB hot_outlet annotation(
        Placement(visible = true, transformation(origin = {-98, -60}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-62, -98}, extent = {{-18, -18}, {18, 18}}, rotation = 0)));
    
    initial equation
      Qt[2:n] * (N-1)= wc * (hc[2:n] - hc[1:n - 1]) "cold fluid";
      Qt[1:n - 1] * (N-1)= wh * (hh[2:n] - hh[1:n - 1]) "cold fluid";
      
    equation
//Heat transfer
      Qt[1:n] = W * dx * h_t * homotopy(Th[1:n] - Tc[1:n], Q_nom);

//Energy equtions
      for i in 1:n - 1 loop
        (N-1) * Qt[i+1] + wc * (hc[i] - hc[i+1]) = state_c[i+1].d * dV * div(N, 2) * der(hc[i+1]) "cold fluid";
        -(N-1) * Qt[i] + wh * (hh[i+1] - hh[i]) = state_h[i].d * dV * div(N+1, 2) * der(hh[i])  "hot fluid";
      end for;
    
//Flange equations
      p_cold = cold_inlet.p;
      cold_inlet.p = cold_outlet.p;
      p_hot = hot_inlet.p;
      hot_inlet.p = hot_outlet.p;
      
      cold_inlet.m_flow + cold_outlet.m_flow = 0;
      wc = homotopy(cold_inlet.m_flow, wc_nom);
      hot_inlet.m_flow + hot_outlet.m_flow = 0;
      wh = homotopy(hot_inlet.m_flow, wh_nom);
      
      hc_in = inStream(cold_inlet.h_outflow);
      hh_in = inStream(hot_inlet.h_outflow);
      
      cold_inlet.h_outflow = hc[1];
      cold_outlet.h_outflow = hc[n];
      
      hot_inlet.h_outflow = hh[n];
      hot_outlet.h_outflow = hh[1];
      
      annotation(
        __OpenModelica_commandLineOptions = "--matchingAlgorithm=PFPlusExt --indexReductionMethod=dynamicStateSelection -d=initialization,NLSanalyticJacobian -d=bltdump",
        __OpenModelica_simulationFlags(lv = "LOG_STATS", s = "dassl"),
        Icon);
    end DynamicHX;

    model Pump
  package Medium = Modelica.Media.Water.StandardWater;
      parameter Real eta_iso = 0.85;
      parameter SI.Pressure P_in = 0.08e5;
      parameter SI.Pressure P_out = 8e5;
      Medium.ThermodynamicState state_in;
      Medium.ThermodynamicState state_iso;
      Medium.SpecificEnthalpy h_in;
      Medium.SpecificEntropy s_in;
      Medium.SpecificEnthalpy h_out;
      Medium.SpecificEnthalpy h_iso;
      Medium.SpecificInternalEnergy u_in;
      SI.SpecificVolume v_in;
      ThermoPower.Water.FlangeA inlet annotation(
        Placement(visible = true, transformation(origin = {-92, 88}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-95, -3}, extent = {{-19, -19}, {19, 19}}, rotation = 0)));
      ThermoPower.Water.FlangeB outlet annotation(
        Placement(visible = true, transformation(origin = {90, 88}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {97, -1}, extent = {{-19, -19}, {19, 19}}, rotation = 0)));
    equation
      abs(h_iso - h_in) = v_in * (outlet.p - inlet.p);
      s_in = Medium.specificEntropy(state_in);
      state_iso = Medium.setState_ps(outlet.p, s_in);
      h_iso = Medium.specificEnthalpy(state_iso);
      eta_iso = (h_iso - h_in) / (h_out - h_in);
      inlet.m_flow + outlet.m_flow = 0;
      state_in = Medium.setState_phX(inlet.p, h_in);
      h_in = inStream(inlet.h_outflow);
      h_out = outlet.h_outflow;
      inlet.h_outflow = outlet.h_outflow;
      u_in = Medium.specificInternalEnergy(state_in);
      v_in = (h_in - u_in) / inlet.p;
      annotation(
        Icon(graphics = {Text(extent = {{-100, -64}, {100, -90}}, textString = "%name"), Ellipse(origin = {1, 4}, fillPattern = FillPattern.HorizontalCylinder, extent = {{-77, 78}, {77, -78}}), Polygon(origin = {-11, -2}, fillColor = {255, 255, 127}, fillPattern = FillPattern.Solid, points = {{-19, 46}, {-19, -38}, {65, 6}, {-19, 46}})}));
  end Pump;

    model StodolaTurbine "Steam turbine"
      package Medium = ThermoPower.Water.StandardWater;
      
      parameter SI.PerUnit eta_iso_nom=0.92 "Nominal isentropic efficiency";
      parameter SI.Area Kt "Kt coefficient of Stodola's law";
      Medium.Density rho "Inlet density"; 
        
      parameter Boolean explicitIsentropicEnthalpy=true
        "Outlet enthalpy computed by isentropicEnthalpy function";
      parameter Medium.MassFlowRate wstart=wnom "Mass flow rate start value"
        annotation (Dialog(tab="Initialisation"));
      parameter SI.PerUnit PRstart "Pressure ratio start value"
        annotation (Dialog(tab="Initialisation"));
      parameter Medium.MassFlowRate wnom "Inlet nominal flowrate";
      parameter Medium.AbsolutePressure pnom "Nominal inlet pressure";
      parameter Real eta_mech=0.98 "Mechanical efficiency";
      parameter Boolean usePartialArcInput = false
        "Use the input connector for the partial arc opening";
    
      outer ThermoPower.System system "System wide properties";
    
      Medium.ThermodynamicState steamState_in;
      Medium.ThermodynamicState steamState_iso;
    
      SI.Angle phi(start=0) "shaft rotation angle";
      SI.Torque tau "net torque acting on the turbine";
      SI.AngularVelocity omega "shaft angular velocity";
      Medium.MassFlowRate w "Mass flow rate";
      Medium.SpecificEnthalpy hin "Inlet enthalpy";
      Medium.SpecificEnthalpy hout "Outlet enthalpy";
      Medium.SpecificEnthalpy hiso "Isentropic outlet enthalpy";
      Medium.SpecificEntropy sin "Inlet entropy";
      Medium.AbsolutePressure pin(start=pnom) "Outlet pressure";
      Medium.AbsolutePressure pout(start=pnom/PRstart) "Outlet pressure";
      SI.PerUnit PR "pressure ratio";
      SI.Power Pm "Mechanical power input";
      SI.PerUnit eta_iso "Isentropic efficiency";
      SI.PerUnit theta "Partial arc opening in p.u."; 
    
      Modelica.Blocks.Interfaces.RealInput partialArc if usePartialArcInput
        "Partial arc opening in p.u." annotation (Placement(
            transformation(extent={{-60,-50},{-40,-30}}, rotation=0)));
      Modelica.Mechanics.Rotational.Interfaces.Flange_a shaft_a annotation (
          Placement(transformation(extent={{-76,-10},{-56,10}}, rotation=0)));
      Modelica.Mechanics.Rotational.Interfaces.Flange_b shaft_b annotation (
          Placement(transformation(extent={{54,-10},{74,10}}, rotation=0)));
      ThermoPower.Water.FlangeA inlet(redeclare package Medium = Medium, m_flow(min=0))
        annotation (Placement(transformation(extent={{-100,60},{-60,100}},
              rotation=0)));
      ThermoPower.Water.FlangeB outlet(redeclare package Medium = Medium, m_flow(max=0))
        annotation (Placement(transformation(extent={{60,60},{100,100}},
              rotation=0)));
    protected
      Modelica.Blocks.Interfaces.RealInput partialArc_int "Internal connector for partial arc input";
    equation
      PR = pin/pout "Pressure ratio";
      theta = partialArc_int;
  if not usePartialArcInput then
        partialArc_int = 1 "Default value if not connector input is disabled";
// otherwise partialArc_int is connected to partialArc input connector
      end if;
  if explicitIsentropicEnthalpy then
        hiso = Medium.isentropicEnthalpy(pout, steamState_in) "Isentropic enthalpy";
//dummy assignments
        sin = 0;
        steamState_iso = Medium.setState_phX(Medium.p_default, Medium.h_default);
      else
        sin = Medium.specificEntropy(steamState_in);
        steamState_iso = Medium.setState_ps(pout, sin);
        hiso = Medium.specificEnthalpy(steamState_iso);
      end if;
      hin - hout = eta_iso*(hin - hiso) "Computation of outlet enthalpy";
      Pm = eta_mech*w*(hin - hout) "Mechanical power from the steam";
      Pm = -tau*omega "Mechanical power balance";
    
      inlet.m_flow + outlet.m_flow = 0 "Mass balance";
    
      assert(w >= -wnom/100, "The turbine model does not support flow reversal");
// Mechanical boundary conditions
      shaft_a.phi = phi;
      shaft_b.phi = phi;
      shaft_a.tau + shaft_b.tau = tau;
      der(phi) = omega;
// steam boundary conditions and inlet steam properties
      steamState_in = Medium.setState_phX(pin, inStream(inlet.h_outflow));
      hin = inStream(inlet.h_outflow);
      hout = outlet.h_outflow;
      pin = inlet.p;
      pout = outlet.p;
      w = inlet.m_flow;
    
      rho =  Medium.density(steamState_in);
      w = homotopy(Kt*theta*sqrt(pin*rho)*ThermoPower.Functions.sqrtReg(1 - (1/PR)^2),
                   theta*wnom/pnom*pin) "Stodola's law";
      eta_iso = eta_iso_nom "Constant efficiency";
// The next equation is provided to close the balance but never actually used
      inlet.h_outflow = outlet.h_outflow;
      connect(partialArc, partialArc_int);
      annotation (
        Icon(graphics={Polygon( fillPattern = FillPattern.Solid, lineThickness = 0.5, points = {{-28, 76}, {-28, 28}, {-22, 28}, {-22, 82}, {-60, 82}, {-60, 76}, {-28, 76}}),
            Polygon( fillPattern = FillPattern.Solid, lineThickness = 0.5, points = {{26, 56}, {32, 56}, {32, 76}, {60, 76}, {60, 82}, {26, 82}, {26, 56}}),
            Rectangle(fillColor = {160, 160, 164}, fillPattern = FillPattern.Sphere, extent = {{-60, 8}, {60, -8}}),
            Polygon(fillPattern = FillPattern.Solid, lineThickness = 0.5, points = {{-28, 28}, {-28, -26}, {32, -60}, {32, 60}, {-28, 28}}),
            Text(extent={{-130,-60},{128,-100}}, textString="%name")}),
        Diagram(graphics),
        Documentation(info="<html>
    <p>This base model contains the basic interface, parameters and definitions for steam turbine models. It lacks the actual performance characteristics, i.e. two more equations to determine the flow rate and the efficiency.
    <p>This model does not include any shaft inertia by itself; if that is needed, connect a <tt>Modelica.Mechanics.Rotational.Inertia</tt> model to one of the shaft connectors.
    <p><b>Modelling options</b></p>
    <p>The following options are available to calculate the enthalpy of the outgoing steam:
    <ul><li><tt>explicitIsentropicEnthalpy = true</tt>: the isentropic enthalpy <tt>hout_iso</tt> is calculated by the <tt>Medium.isentropicEnthalpy</tt> function. <li><tt>explicitIsentropicEnthalpy = false</tt>: the isentropic enthalpy is given equating the specific entropy of the inlet steam <tt>steam_in</tt> and of a fictional steam state <tt>steam_iso</tt>, which has the same pressure of the outgoing steam, both computed with the function <tt>Medium.specificEntropy</tt>.</pp></ul>
    </html>", revisions="<html>
    <ul>
    <li><i>20 Apr 2005</i>
    by <a href=\"mailto:francesco.casella@polimi.it\">Francesco Casella</a>:<br>
       First release.</li>
    <li><i>5 Oct 2011</i>
    by <a href=\"mailto:francesco.casella@polimi.it\">Francesco Casella</a>:<br>
       Small changes in alias variables.</li>
    </ul>
    </html>"));
    end StodolaTurbine;
  end Components;

  package Functions
    function sigmoid
      import Modelica.Math.exp;
      input Real x;
      input Real y1;
      input Real y2;
      input Real c = 0;
      input Real k = 1;
      output Real y;
    algorithm
      y := y1 + (y2 - y1) / (1 + exp(-k * (x - c)));
    end sigmoid;

    function ht_formula
      package Medium = Modelica.Media.Water.StandardWater;
      input Medium.ThermodynamicState state;
      input SI.Velocity v;
      output SI.CoefficientOfHeatTransfer h;
    protected
      Real C "Nusselt's number coefficient";
      parameter SI.Length t = 0.01 "Thickness";
      parameter Real n_t = 1 / 3;
      parameter Real m_t = 0.5;
      Real Nu "Nusselt's number";
      Real Re "Reynold's number";
      Real Pr "Prandtl's number";
      Medium.ThermalConductivity k;
      Medium.Density rho;
      Medium.DynamicViscosity mu;
      Real x;
    algorithm
      rho := Medium.density(state);
      mu := Medium.dynamicViscosity(state);
      k := Medium.thermalConductivity(state);
      Re := rho * v * 2 * 0.01 / mu;
      Pr := Medium.prandtlNumber(state);
      x := Medium.vapourQuality(state);
      if x == 0 then
        C := 1.29;
      elseif x == 1 then
        C := 0.063;
      else
        C := 1.29 - x * (1.29 - 0.063);
      end if;
      Nu := C * Re ^ m_t * Pr ^ n_t;
      h := k * Nu / t;
    end ht_formula;

    function ht_formula_smooth
      package Medium = Modelica.Media.Water.StandardWater;
      input Medium.ThermodynamicState state "State of the fluid";
      input Real w_a "Mass flow rate per unit area";
      output SI.CoefficientOfHeatTransfer ht_res "result ht value";
    protected
      Medium.SaturationProperties sat;
      Medium.SpecificEnthalpy h;
      Medium.SpecificEnthalpy hf;
      Medium.SpecificEnthalpy hg;
      SI.CoefficientOfHeatTransfer ht;
      SI.CoefficientOfHeatTransfer ht_f;
      SI.CoefficientOfHeatTransfer ht_g;
      SI.Velocity v;
      SI.Velocity v_0;
      SI.Velocity v_1;
      Medium.ThermodynamicState state_f;
      Medium.ThermodynamicState state_g;
      Real x;
    algorithm
      sat.psat := Medium.pressure(state);
      sat.Tsat := Medium.saturationPressure(sat.psat);
      h := Medium.specificEnthalpy(state);
      hf := Medium.bubbleEnthalpy(sat);
      hg := Medium.dewEnthalpy(sat);
      x := (h - hf) / (hg - hf);
      state_f := Medium.setState_ph(Medium.pressure(state), hf);
      state_g := Medium.setState_ph(Medium.pressure(state), hg);
      v := w_a / Medium.density(state);
      v_0 := v * Medium.density(state) / Medium.bubbleDensity(sat);
      v_1 := v * Medium.density(state) / Medium.dewDensity(sat);
      ht := ht_coef(state, v);
      ht_f := ht_coef(state_f, v_0);
      ht_g := ht_coef(state_g, v_1);
      if x > 0.3 and x < 1.001 then
        ht_res := ht_f + x * (ht_g - ht_f);
      else
        ht_res := ht;
      end if;
    end ht_formula_smooth;

    function ht_formula_sigmoid
      //take two ref points and use sigmoid
      package Medium = Modelica.Media.Water.StandardWater;
      input Medium.ThermodynamicState state "State of the fluid";
      input Real VxRho "Mass flow rate per unit area";
      output SI.CoefficientOfHeatTransfer ht_res;
    protected
      parameter Real k;
      parameter Real d;
      Medium.SaturationProperties sat;
      Medium.SpecificEnthalpy h;
      Medium.SpecificEnthalpy hf;
      Medium.SpecificEnthalpy hg;
      SI.CoefficientOfHeatTransfer ht;
      SI.CoefficientOfHeatTransfer ht_f;
      SI.CoefficientOfHeatTransfer ht_g;
      SI.Velocity v;
      SI.Velocity v_0;
      SI.Velocity v_1;
      Medium.ThermodynamicState state_f;
      Medium.ThermodynamicState state_g;
      Real x;
    algorithm
      sat.psat := Medium.pressure(state);
      sat.Tsat := Medium.saturationPressure(sat.psat);
      h := Medium.specificEnthalpy(state);
      hf := Medium.bubbleEnthalpy(sat);
      hg := Medium.dewEnthalpy(sat);
      x := (h - hf) / (hg - hf);
      state_f := Medium.setState_ph(Medium.pressure(state), hf);
      state_g := Medium.setState_ph(Medium.pressure(state), hg);
      v := VxRho / Medium.density(state);
      v_0 := v * Medium.density(state) / Medium.bubbleDensity(sat);
      v_1 := v * Medium.density(state) / Medium.dewDensity(sat);
      ht_f := ht_coef(state_f, v_0);
      ht_g := ht_coef(state_g, v_1);
      ht_res := sigmoid(x, ht_f, ht_g, 0.1, 35);
    end ht_formula_sigmoid;

    function ht_ref_sigmoid
      package Medium = Modelica.Media.Water.StandardWater;
      input Medium.ThermodynamicState state "State of the fluid";
      input SI.CoefficientOfHeatTransfer ht_f;
      input SI.CoefficientOfHeatTransfer ht_g;
      output SI.CoefficientOfHeatTransfer ht_res "result ht value";
    protected
      parameter Real k;
      parameter Real d;
      Medium.SaturationProperties sat;
      Medium.SpecificEnthalpy h;
      Medium.SpecificEnthalpy hf;
      Medium.SpecificEnthalpy hg;
      SI.CoefficientOfHeatTransfer ht;
      Real x;
    algorithm
      sat.psat := Medium.pressure(state);
      sat.Tsat := Medium.saturationPressure(sat.psat);
      h := Medium.specificEnthalpy(state);
      hf := Medium.bubbleEnthalpy(sat);
      hg := Medium.dewEnthalpy(sat);
      x := (h - hf) / (hg - hf);
      ht_res := sigmoid(x, ht_f, ht_g, 0.1, 35);
    end ht_ref_sigmoid;

    package tests
      model ht_formula_test
        package Medium = Modelica.Media.Water.StandardWater;
        Medium.SpecificEnthalpy h(start = 10e3);
        Medium.ThermodynamicState state;
        parameter SI.Pressure p = 1e5;
        parameter SI.MassFlowRate w = 0.5;
        parameter SI.Area A = Modelica.Constants.pi * 0.01 ^ 2;
        SI.CoefficientOfHeatTransfer ht;
        Real x;
      equation
        der(h) = 500;
        state = Medium.setState_phX(p, h);
        ht = ht_formula(state, w / A);
        x = Medium.vapourQuality(state);
      end ht_formula_test;

      model ht_formula_sigmoid_test
        package Medium = Modelica.Media.Water.StandardWater;
        Medium.SpecificEnthalpy h(start = 10e3);
        Medium.ThermodynamicState state;
        parameter SI.Pressure p = 1e5;
        parameter SI.MassFlowRate w = 0.5;
        parameter SI.Area A = Modelica.Constants.pi * 0.05 ^ 2;
        SI.CoefficientOfHeatTransfer ht;
        Real x;
      equation
        der(h) = 500;
        state = Medium.setState_phX(p, h);
        ht = ht_sigmoid(state, w / A);
        x = Medium.vapourQuality(state);
      end ht_formula_sigmoid_test;

      model ht_formula_smooth_test
        package Medium = Modelica.Media.Water.StandardWater;
        Medium.SpecificEnthalpy h(start = 10e3);
        Medium.ThermodynamicState state;
        parameter SI.Pressure p = 1e5;
        parameter SI.MassFlowRate w = 0.5;
        parameter SI.Area A = Modelica.Constants.pi * 0.01 ^ 2;
        SI.CoefficientOfHeatTransfer ht;
        Real x;
      equation
        der(h) = 500;
        state = Medium.setState_phX(p, h);
        ht = ht_smooth(state, w / A);
        x = Medium.vapourQuality(state);
      end ht_formula_smooth_test;

      model ht_ref_sigmoid_test
        package Medium = Modelica.Media.Water.StandardWater;
        Medium.SpecificEnthalpy h(start = 10e3);
        Medium.ThermodynamicState state;
        parameter SI.Pressure p = 1e5;
        parameter SI.MassFlowRate w = 0.5;
        parameter SI.Area A = Modelica.Constants.pi * 0.01 ^ 2;
        SI.CoefficientOfHeatTransfer ht;
        Real x;
      equation
        der(h) = 500;
        state = Medium.setState_phX(p, h);
        ht = HT_coeff.functions.ht_ref_sigmoid(state, 120e3, 200e3);
        x = Medium.vapourQuality(state);
      end ht_ref_sigmoid_test;
    end tests;
  end Functions;

  package Utilities
    model EnthalpySensor
      Modelica.Blocks.Interfaces.RealOutput h(unit = "J/kg") "Absolute temperature as output signal" annotation(
        Placement(visible = true, transformation(extent = {{82, 42}, {102, 62}}, rotation = 0), iconTransformation(extent = {{90, -10}, {110, 10}}, rotation = 0)));
      ThermoPower.Water.FlangeA inlet annotation(
        Placement(visible = true, transformation(origin = {-100, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-100, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      ThermoPower.Water.FlangeB outlet annotation(
        Placement(visible = true, transformation(origin = {100, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {118, -38}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    equation
      h = outlet.h_outflow;
      inlet.m_flow + outlet.m_flow = 0;
      inlet.p = outlet.p;
      inlet.h_outflow = inStream(outlet.h_outflow);
      inStream(inlet.h_outflow) = outlet.h_outflow;
      annotation(
        Diagram(coordinateSystem(preserveAspectRatio = true, extent = {{-100, -100}, {100, 100}}), graphics = {Ellipse(fillColor = {191, 0, 0}, fillPattern = FillPattern.Solid, lineThickness = 0.5, extent = {{-20, -98}, {20, -60}}), Rectangle(lineColor = {191, 0, 0}, fillColor = {191, 0, 0}, fillPattern = FillPattern.Solid, extent = {{-12, 40}, {12, -68}}), Line(points = {{12, 0}, {90, 0}}, color = {0, 0, 255}), Line(points = {{-94, 0}, {-14, 0}}, color = {191, 0, 0}), Polygon(lineThickness = 0.5, points = {{-12, 40}, {-12, 80}, {-10, 86}, {-6, 88}, {0, 90}, {6, 88}, {10, 86}, {12, 80}, {12, 40}, {-12, 40}}), Line(points = {{-12, 40}, {-12, -64}}, thickness = 0.5), Line(points = {{12, 40}, {12, -64}}, thickness = 0.5), Line(points = {{-40, -20}, {-12, -20}}), Line(points = {{-40, 20}, {-12, 20}}), Line(points = {{-40, 60}, {-12, 60}}), Text(extent = {{102, -28}, {60, -78}}, textString = "h")}),
        Icon(coordinateSystem(preserveAspectRatio = true, extent = {{-100, -100}, {100, 100}}), graphics = {Ellipse(fillColor = {191, 0, 0}, fillPattern = FillPattern.Solid, lineThickness = 0.5, extent = {{-20, -98}, {20, -60}}), Rectangle(lineColor = {191, 0, 0}, fillColor = {191, 0, 0}, fillPattern = FillPattern.Solid, extent = {{-12, 40}, {12, -68}}), Line(points = {{12, 0}, {90, 0}}, color = {0, 0, 255}), Line(points = {{-90, 0}, {-12, 0}}, color = {191, 0, 0}), Polygon(lineThickness = 0.5, points = {{-12, 40}, {-12, 80}, {-10, 86}, {-6, 88}, {0, 90}, {6, 88}, {10, 86}, {12, 80}, {12, 40}, {-12, 40}}), Line(points = {{-12, 40}, {-12, -64}}, thickness = 0.5), Line(points = {{12, 40}, {12, -64}}, thickness = 0.5), Line(points = {{-40, -20}, {-12, -20}}), Line(points = {{-40, 20}, {-12, 20}}), Line(points = {{-40, 60}, {-12, 60}}), Text(extent = {{126, -20}, {26, -120}}, textString = "h"), Text(lineColor = {0, 0, 255}, extent = {{-150, 130}, {150, 90}}, textString = "%name")}),
        Documentation(info = "<html>
    <p>
    This is an ideal absolute temperature sensor which returns
    the temperature of the connected port in Kelvin as an output
    signal.  The sensor itself has no thermal interaction with
    whatever it is connected to.  Furthermore, no
    thermocouple-like lags are associated with this
    sensor model.
    </p>
    </html>"));
    end EnthalpySensor;

    model PID
      Modelica.Blocks.Interfaces.RealInput y annotation(
        Placement(visible = true, transformation(origin = {-100, -50}, extent = {{-20, -20}, {20, 20}}, rotation = 0), iconTransformation(origin = {-80, -40}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
      Modelica.Blocks.Interfaces.RealInput r annotation(
        Placement(visible = true, transformation(origin = {-100, 50}, extent = {{-20, -20}, {20, 20}}, rotation = 0), iconTransformation(origin = {-80, 40}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
      Modelica.Blocks.Interfaces.RealOutput x annotation(
        Placement(visible = true, transformation(origin = {99, -1}, extent = {{-19, -19}, {19, 19}}, rotation = 0), iconTransformation(origin = {90, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      parameter Boolean proportional = true "if false, inversely proportial";
      parameter Real x_0;
      parameter Real Kp;
      parameter Real b;
      parameter Real Ki;
      Real e;
      Real P;
      Real I;
    initial equation
      x = x_0;
    equation
      if proportional then
        e = r - y;
      else
        e = y - r;
      end if;
      P = Kp * e + b;
      der(I) = Ki * e;
      x = P + I;
      annotation(
        Diagram,
        Icon(graphics = {Rectangle(origin = {0, -1}, fillPattern = FillPattern.VerticalCylinder, extent = {{-70, 67}, {70, -67}}), Text(origin = {1, 1}, lineColor = {255, 255, 255}, fillColor = {255, 255, 255}, extent = {{-39, 27}, {39, -27}}, textString = "PID")}));
    end PID;

    model FixedP
      ThermoPower.Water.FlangeB flange annotation(
        Placement(visible = true, transformation(extent = {{80, -20}, {120, 20}}, rotation = 0), iconTransformation(extent = {{80, -20}, {120, 20}}, rotation = 0)));
      parameter SI.Pressure P;
    equation
      flange.m_flow = 0;
      flange.p = P;
      annotation(
        Icon(graphics = {Ellipse(fillPattern = FillPattern.Solid, extent = {{-80, 80}, {80, -80}}), Text(origin = {5, 0}, lineColor = {255, 255, 255}, fillColor = {255, 255, 255}, extent = {{-53, 42}, {53, -42}}, textString = "FixedP"), Text(origin = {2, -88}, extent = {{-72, 14}, {72, -14}}, textString = "%name")}));
    end FixedP;
  end Utilities;

  package CombinationSteps
    model static_01
      package Medium = Modelica.Media.Water.StandardWater;
      parameter Integer n = 50 "Number of finite divisions";
      parameter Integer N = 30 "Number of tubes";
      parameter SI.Length L = 5 "Length of the HX";
      parameter SI.Length r = 0.05 "Radius of the tubes";
      final parameter SI.Length dx = L / n;
      //  parameter SI.CoefficientOfHeatTransfer h_t = 170 "Heat transfer coefficient";
      SI.CoefficientOfHeatTransfer h_t[n] "Heat transfer coefficient";
      SI.Temperature Th[n] "Hot side temperature";
      SI.Temperature Tc[n] "Cold side temperature";
      Real dT[n] "Cold side temperature";
      Medium.SpecificEnthalpy hh[n](each start = hh_in * 0.9) "Hot side specific enthalpy";
      Medium.SpecificEnthalpy hc[n](each start = hc_in * 2) "Cold side specific enthalpy";
      Medium.ThermodynamicState state_h[n];
      Medium.ThermodynamicState state_c[n];
      SI.Power Qt[n](min = fill(0, n)) "Heat transfer at slice i";
      Real x_h[n] "Hot side vapour quality";
      Real x_c[n] "Hot side vapour quality";
      //The following variable and parameters can be substituted with connectors
      parameter Medium.MassFlowRate wh = wc "Hot side mass flow rate (total)";
      parameter Medium.MassFlowRate wc = 0.5 "Cold side mass flow rate";
      parameter Medium.SpecificEnthalpy hh_in = 3004772 "Hot side inlet enthalpy";
      parameter Medium.SpecificEnthalpy hc_in = 150000 "Cold side inlet enthalpy";
      parameter SI.Pressure p_hot = 1e5 "Hot side Pressure";
      parameter SI.Pressure p_cold = 1e5 "Cold side pressure";
    equation
//Heat transfer Q_i ( the second value in homotopy is an approximation for dT_ave )
      for i in 1:n loop
        if time < 1 then
          Qt[i] = homotopy(2 * pi * r * dx * h_t[i] * 50, 0);
          h_t[i] = 16e3;
        else
          Qt[i] = homotopy(2 * pi * r * dx * h_t[i] * N * dT[i], 2 * pi * r * dx * h_t[i] * dT[i]);
          h_t[i] = homotopy(Functions.ht_ref_sigmoid(state_c[i], 120e3, 200e3), 160e3);
        end if;
      end for;
//Energy conservation for cold fluid
      Qt[2:n] = wc * (hc[2:n] - hc[1:n - 1]);
//Energy conservation for hot fluid
      Qt[1:n - 1] = wh * (hh[2:n] - hh[1:n - 1]);
//Setting the states
      state_h = Medium.setState_phX(p_hot, hh);
      state_c = Medium.setState_phX(p_cold, hc);
//Extracting the temperatures from the states
      Th = Medium.temperature(state_h);
      Tc = Medium.temperature(state_c);
      dT = Th - Tc;
//Extracting the qualities from the states
      x_h = Medium.vapourQuality(state_h);
      x_c = Medium.vapourQuality(state_c);
//Boundary conditions
      hc[1] = hc_in;
      hh[n] = hh_in;
      annotation(
        __OpenModelica_commandLineOptions = "--matchingAlgorithm=PFPlusExt --indexReductionMethod=dynamicStateSelection -d=initialization,NLSanalyticJacobian -d=bltdump",
        __OpenModelica_simulationFlags(lv = "LOG_STATS", s = "dassl"),
        Documentation(info = "<html><head></head><body>Static<div>Shell and tube</div><div>Constant HT coefficient</div><div>Incompressible</div><div>All phases</div><div><br></div><div>7/27/2021</div><div>Meraj Mammadov</div></body></html>"));
    end static_01;

    model dynamic_01
      ThermoPower.Water.SinkPressure sinkPressure1(p0 = 80e5) annotation(
        Placement(visible = true, transformation(origin = {-68, -14}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
      ThermoPower.Water.SourceMassFlow sourceMassFlow(h = 183.36e3, use_in_w0 = false, w0 = 123.58) annotation(
        Placement(visible = true, transformation(origin = {-70, 26}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      ThermoPower.Water.SinkPressure sinkPressure(p0 = 80e5) annotation(
        Placement(visible = true, transformation(origin = {76, 36}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      ThermoPower.Water.SourceMassFlow sourceMassFlow1(h = 3004772, use_T = false, use_in_w0 = true, w0 = 123) annotation(
        Placement(visible = true, transformation(origin = {26, -36}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
      inner ThermoPower.System system annotation(
        Placement(visible = true, transformation(origin = {86, 86}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Components.DynamicHX dynamicHX(N = 100, Q_nom = 1, h_t = 100000, n = 50, use_q_nom = true, wc_nom = 123, wh_nom = 123) annotation(
        Placement(visible = true, transformation(origin = {-4, 4}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Utilities.EnthalpySensor enthalpySensor annotation(
        Placement(visible = true, transformation(origin = {36, 54}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Utilities.PID pid(Ki = 0.0001, Kp = 0.01, b = 122, x_0 = 123) annotation(
        Placement(visible = true, transformation(origin = {93, -1}, extent = {{17, 17}, {-17, -17}}, rotation = 0)));
      Modelica.Blocks.Sources.Ramp ramp(duration = 100, height = 1638e3, offset = 1.12e6, startTime = 10) annotation(
        Placement(visible = true, transformation(origin = {154, -6}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
    equation
      connect(sourceMassFlow.flange, dynamicHX.cold_inlet) annotation(
        Line(points = {{-60, 26}, {-42, 26}, {-42, 10}, {-14, 10}}, color = {0, 0, 255}));
      connect(sinkPressure1.flange, dynamicHX.hot_outlet) annotation(
        Line(points = {{-58, -14}, {-42, -14}, {-42, -2}, {-14, -2}}, color = {0, 0, 255}));
      connect(sourceMassFlow1.flange, dynamicHX.hot_inlet) annotation(
        Line(points = {{16, -36}, {16, -2}, {6, -2}}, color = {0, 0, 255}));
      connect(enthalpySensor.inlet, dynamicHX.cold_outlet) annotation(
        Line(points = {{26, 54}, {6, 54}, {6, 10}}, color = {0, 0, 255}));
      connect(enthalpySensor.h, pid.y) annotation(
        Line(points = {{46, 54}, {106, 54}, {106, 6}}, color = {0, 0, 127}));
      connect(pid.x, sourceMassFlow1.in_w0) annotation(
        Line(points = {{78, 0}, {30, 0}, {30, -30}}, color = {0, 0, 127}));
      connect(ramp.y, pid.r) annotation(
        Line(points = {{144, -6}, {126, -6}, {126, -8}, {106, -8}}, color = {0, 0, 127}));
      connect(enthalpySensor.outlet, sinkPressure.flange) annotation(
        Line(points = {{48, 50}, {56, 50}, {56, 36}, {66, 36}}, color = {0, 0, 255}));
    protected
    end dynamic_01;

    model boiler_test
      ThermoPower.Water.SinkPressure sinkPressure1(p0 = 80e5) annotation(
        Placement(visible = true, transformation(origin = {-74, 72}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
      ThermoPower.Water.SourceMassFlow sourceMassFlow(h = 183.36e3, use_in_w0 = false, w0 = 123.58) annotation(
        Placement(visible = true, transformation(origin = {-74, 26}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      ThermoPower.Water.SinkPressure sinkPressure(p0 = 80e5) annotation(
        Placement(visible = true, transformation(origin = {46, 30}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      ThermoPower.Water.SourceMassFlow sourceMassFlow1(h = 3004772, use_T = false, w0 = 1350) annotation(
        Placement(visible = true, transformation(origin = {44, 76}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
      inner ThermoPower.System system annotation(
        Placement(visible = true, transformation(origin = {86, 86}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Components.DynamicHX dynamicHX(N = 100, Q_nom = 1, h_t = 100000, n = 50, use_q_nom = true, wc_nom = 123, wh_nom = 123) annotation(
        Placement(visible = true, transformation(origin = {-6, 50}, extent = {{-10, 10}, {10, -10}}, rotation = 0)));
    equation
      connect(sourceMassFlow.flange, dynamicHX.cold_inlet) annotation(
        Line(points = {{-64, 26}, {-42, 26}, {-42, 44}, {-16, 44}}, color = {0, 0, 255}));
      connect(sinkPressure1.flange, dynamicHX.hot_outlet) annotation(
        Line(points = {{-64, 72}, {-36, 72}, {-36, 56}, {-16, 56}}, color = {0, 0, 255}));
      connect(sinkPressure.flange, dynamicHX.cold_outlet) annotation(
        Line(points = {{36, 30}, {16, 30}, {16, 44}, {4, 44}}, color = {0, 0, 255}));
      connect(sourceMassFlow1.flange, dynamicHX.hot_inlet) annotation(
        Line(points = {{34, 76}, {18, 76}, {18, 56}, {4, 56}}, color = {0, 0, 255}));
    protected
    end boiler_test;

    model turbine_test
      inner ThermoPower.System system annotation(
        Placement(visible = true, transformation(origin = {86, 86}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      ThermoPower.Water.SourcePressure sourcePressure(h = 2758e3, p0 = 80e5) annotation(
        Placement(visible = true, transformation(origin = {-68, 22}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      ThermoPower.Water.SinkPressure sinkPressure(h = 1939.3e3, p0 = 0.08e5) annotation(
        Placement(visible = true, transformation(origin = {52, 22}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Mechanics.Rotational.Sources.ConstantSpeed constantSpeed(w_fixed = 3000 / 60 * 3.14159) annotation(
        Placement(visible = true, transformation(extent = {{46, -34}, {26, -14}}, rotation = 0)));
  Rankine.Components.StodolaTurbine stodolaTurbine(Kt = 0.0067, eta_iso_nom = 0.85, phi(displayUnit = "rad"), pin(start = 80e5), pnom = 80e5, pout(start = 0.08e5), wnom = 123)  annotation(
        Placement(visible = true, transformation(origin = {-12, 14}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    equation
  connect(sourcePressure.flange, stodolaTurbine.inlet) annotation(
        Line(points = {{-58, 22}, {-20, 22}}, color = {0, 0, 255}));
  connect(stodolaTurbine.outlet, sinkPressure.flange) annotation(
        Line(points = {{-4, 22}, {42, 22}}, color = {0, 0, 255}));
  connect(stodolaTurbine.shaft_b, constantSpeed.flange) annotation(
        Line(points = {{-6, 14}, {12, 14}, {12, -24}, {26, -24}}));
    end turbine_test;

    model boiler_and_turbine_test
      ThermoPower.Water.SourceMassFlow sourceMassFlow(h = 183.36e3, use_in_w0 = false, w0 = 123.58) annotation(
        Placement(visible = true, transformation(origin = {-74, 26}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      inner ThermoPower.System system annotation(
        Placement(visible = true, transformation(origin = {86, 86}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Components.DynamicHX dynamicHX(N = 100, Q_nom = 1, h_t = 100000, n = 50, wc_nom = 123, wh_nom = 123) annotation(
        Placement(visible = true, transformation(origin = {-6, 50}, extent = {{-10, 10}, {10, -10}}, rotation = 0)));
      ThermoPower.Water.SinkPressure sinkPressure1(p0 = 80e5) annotation(
        Placement(visible = true, transformation(origin = {-74, 72}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
      ThermoPower.Water.SourceMassFlow sourceMassFlow1(h = 3004772, use_T = false, w0 = 1350) annotation(
        Placement(visible = true, transformation(origin = {44, 76}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
      ThermoPower.Water.SteamTurbineStodola steamTurbineStodola(Kt = 0.0067, eta_iso_nom = 0.85, phi(displayUnit = "rad"), pin(start = 80e5), pnom = 80e5, pout(start = 0.08e5), wnom = 123.5) annotation(
        Placement(visible = true, transformation(origin = {48, 14}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Mechanics.Rotational.Sources.ConstantSpeed constantSpeed(w_fixed = 3000 / 60 * 3.14159) annotation(
        Placement(visible = true, transformation(extent = {{96, -28}, {76, -8}}, rotation = 0)));
      ThermoPower.Water.SinkPressure sinkPressure2(h = 1939.3e3, p0 = 0.08e5) annotation(
        Placement(visible = true, transformation(origin = {112, 28}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    equation
      connect(sinkPressure1.flange, dynamicHX.hot_outlet) annotation(
        Line(points = {{-64, 72}, {-36, 72}, {-36, 56}, {-16, 56}}, color = {0, 0, 255}));
      connect(sourceMassFlow.flange, dynamicHX.cold_inlet) annotation(
        Line(points = {{-64, 26}, {-42, 26}, {-42, 44}, {-16, 44}}, color = {0, 0, 255}));
      connect(sourceMassFlow1.flange, dynamicHX.hot_inlet) annotation(
        Line(points = {{34, 76}, {18, 76}, {18, 56}, {4, 56}}, color = {0, 0, 255}));
      connect(sinkPressure2.flange, steamTurbineStodola.outlet) annotation(
        Line(points = {{102, 28}, {56, 28}, {56, 22}}, color = {0, 0, 255}));
      connect(constantSpeed.flange, steamTurbineStodola.shaft_b) annotation(
        Line(points = {{76, -18}, {54, -18}, {54, 14}}));
      connect(dynamicHX.cold_outlet, steamTurbineStodola.inlet) annotation(
        Line(points = {{4, 44}, {40, 44}, {40, 22}}, color = {0, 0, 255}));
    end boiler_and_turbine_test;

    model condenser_test_pid
      ThermoPower.Water.SourceMassFlow sourceMassFlow(T = 288.15, use_T = true, use_in_w0 = true, w0 = 123.58) annotation(
        Placement(visible = true, transformation(origin = {-92, -84}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      inner ThermoPower.System system annotation(
        Placement(visible = true, transformation(origin = {86, 82}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Components.DynamicHX dynamicHX(L = 5, N = 100, Q_nom = 1, h_t = 100000, n = 50, use_q_nom = true, wc_nom = 123, wh_nom = 123) annotation(
        Placement(visible = true, transformation(origin = {-6, -44}, extent = {{-10, 10}, {10, -10}}, rotation = 0)));
      ThermoPower.Water.SinkPressure sinkPressure1(p0 = 0.08e5) annotation(
        Placement(visible = true, transformation(origin = {-110, -34}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
      ThermoPower.Water.SinkPressure sinkPressure(p0 = 1e5) annotation(
        Placement(visible = true, transformation(origin = {68, -88}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      ThermoPower.Water.SourceMassFlow sourceMassFlow1(h = 1939.3e3, use_T = false, w0 = 123) annotation(
        Placement(visible = true, transformation(origin = {72, -18}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.Ramp ramp(duration = 100, height = -1.62612e6, offset = 1.8e6, startTime = 10) annotation(
        Placement(visible = true, transformation(origin = {-248, -38}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Utilities.PID pid(Ki = 0.001, Kp = 0.001, b = 122, proportional = false, x_0 = 123) annotation(
        Placement(visible = true, transformation(origin = {-155, -69}, extent = {{-17, 17}, {17, -17}}, rotation = 0)));
      Utilities.EnthalpySensor enthalpySensor annotation(
        Placement(visible = true, transformation(origin = {-68, 6}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
    equation
      connect(sourceMassFlow.flange, dynamicHX.cold_inlet) annotation(
        Line(points = {{-82, -84}, {-49, -84}, {-49, -50}, {-16, -50}}, color = {0, 0, 255}));
      connect(sourceMassFlow1.flange, dynamicHX.hot_inlet) annotation(
        Line(points = {{62, -18}, {13, -18}, {13, -38}, {4, -38}}, color = {0, 0, 255}));
      connect(sinkPressure.flange, dynamicHX.cold_outlet) annotation(
        Line(points = {{58, -88}, {16, -88}, {16, -50}, {4, -50}}, color = {0, 0, 255}));
      connect(enthalpySensor.h, pid.y) annotation(
        Line(points = {{-78, 6}, {-169, 6}, {-169, -62}}, color = {0, 0, 127}));
      connect(ramp.y, pid.r) annotation(
        Line(points = {{-237, -38}, {-203, -38}, {-203, -76}, {-169, -76}}, color = {0, 0, 127}));
      connect(sinkPressure1.flange, enthalpySensor.outlet) annotation(
        Line(points = {{-100, -34}, {-80, -34}, {-80, 2}}, color = {0, 0, 255}));
      connect(enthalpySensor.inlet, dynamicHX.hot_outlet) annotation(
        Line(points = {{-58, 6}, {-16, 6}, {-16, -38}}, color = {0, 0, 255}));
      connect(pid.x, sourceMassFlow.in_w0) annotation(
        Line(points = {{-140, -68}, {-96, -68}, {-96, -78}}, color = {0, 0, 127}));
    end condenser_test_pid;

    model condenser_test
      ThermoPower.Water.SourceMassFlow sourceMassFlow1(h = 1939.3e3, use_T = false, w0 = 123) annotation(
        Placement(visible = true, transformation(origin = {72, -18}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
      ThermoPower.Water.SinkPressure sinkPressure1(p0 = 0.08e5) annotation(
        Placement(visible = true, transformation(origin = {-92, -10}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
      ThermoPower.Water.SourceMassFlow sourceMassFlow(T = 288.15, use_T = true, use_in_w0 = false, w0 = 2472) annotation(
        Placement(visible = true, transformation(origin = {-92, -84}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      ThermoPower.Water.SinkPressure sinkPressure(p0 = 1e5) annotation(
        Placement(visible = true, transformation(origin = {68, -88}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Components.DynamicHX dynamicHX(L = 5, N = 100, Q_nom = 1, h_t = 100000, n = 50, use_q_nom = true, wc_nom = 123, wh_nom = 123) annotation(
        Placement(visible = true, transformation(origin = {-6, -44}, extent = {{-10, 10}, {10, -10}}, rotation = 0)));
      inner ThermoPower.System system annotation(
        Placement(visible = true, transformation(origin = {86, 82}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    equation
      connect(sinkPressure.flange, dynamicHX.cold_outlet) annotation(
        Line(points = {{58, -88}, {16, -88}, {16, -50}, {4, -50}}, color = {0, 0, 255}));
      connect(sourceMassFlow1.flange, dynamicHX.hot_inlet) annotation(
        Line(points = {{62, -18}, {13, -18}, {13, -38}, {4, -38}}, color = {0, 0, 255}));
      connect(sourceMassFlow.flange, dynamicHX.cold_inlet) annotation(
        Line(points = {{-82, -84}, {-49, -84}, {-49, -50}, {-16, -50}}, color = {0, 0, 255}));
      connect(sinkPressure1.flange, dynamicHX.hot_outlet) annotation(
        Line(points = {{-82, -10}, {-16, -10}, {-16, -38}}, color = {0, 0, 255}));
    end condenser_test;

    model boiler_and_turbine_and_condenser_test
      ThermoPower.Water.SinkPressure sinkPressure1(p0 = 80e5) annotation(
        Placement(visible = true, transformation(origin = {-72, 74}, extent = {{8, -8}, {-8, 8}}, rotation = 0)));
      ThermoPower.Water.SourceMassFlow sourceMassFlow1(h = 3004772, use_T = false, w0 = 1350) annotation(
        Placement(visible = true, transformation(origin = {42, 74}, extent = {{8, -8}, {-8, 8}}, rotation = 0)));
      Components.DynamicHX dynamicHX(N = 100, Q_nom = 1, h_t = 100000, n = 50, wc_nom = 123, wh_nom = 123) annotation(
        Placement(visible = true, transformation(origin = {-6, 50}, extent = {{-10, 10}, {10, -10}}, rotation = 0)));
      inner ThermoPower.System system annotation(
        Placement(visible = true, transformation(origin = {86, 86}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Mechanics.Rotational.Sources.ConstantSpeed constantSpeed(w_fixed = 3000 / 60 * 3.14159) annotation(
        Placement(visible = true, transformation(origin = {20, 8}, extent = {{-6, -6}, {6, 6}}, rotation = 0)));
      ThermoPower.Water.SourceMassFlow sourceMassFlow(h = 183.36e3, use_in_w0 = false, w0 = 123.58) annotation(
        Placement(visible = true, transformation(origin = {-73, 27}, extent = {{-9, -9}, {9, 9}}, rotation = 0)));
      Components.DynamicHX dynamicHX1(L = 5, N = 100, Q_nom = 1, h_t = 100000, n = 50, wc_nom = 123, wh_nom = 123) annotation(
        Placement(visible = true, transformation(origin = {8, -50}, extent = {{-10, 10}, {10, -10}}, rotation = 0)));
      ThermoPower.Water.SinkPressure sinkPressure(p0 = 1e5) annotation(
        Placement(visible = true, transformation(origin = {75, -91}, extent = {{-7, -7}, {7, 7}}, rotation = 0)));
      ThermoPower.Water.SourceMassFlow sourceMassFlow3(T = 288.15, use_T = true, use_in_w0 = false, w0 = 2472) annotation(
        Placement(visible = true, transformation(origin = {-80, -88}, extent = {{-8, -8}, {8, 8}}, rotation = 0)));
      ThermoPower.Water.SinkPressure sinkPressure3(p0 = 0.08e5) annotation(
        Placement(visible = true, transformation(origin = {-80, -14}, extent = {{8, -8}, {-8, 8}}, rotation = 0)));
    Components.StodolaTurbine stodolaTurbine(Kt = 0.0067, eta_iso_nom = 0.85, phi(displayUnit = "rad"), pin(start = 80e5), pnom = 80e5, pout(start = 0.08e5), wnom = 123) annotation(
        Placement(visible = true, transformation(origin = {65, 7}, extent = {{-17, -17}, {17, 17}}, rotation = -90)));
    equation
      connect(sourceMassFlow1.flange, dynamicHX.hot_inlet) annotation(
        Line(points = {{34, 74}, {18, 74}, {18, 56}, {4, 56}}, color = {0, 0, 255}));
      connect(sinkPressure1.flange, dynamicHX.hot_outlet) annotation(
        Line(points = {{-64, 74}, {-36, 74}, {-36, 56}, {-16, 56}}, color = {0, 0, 255}));
      connect(sourceMassFlow.flange, dynamicHX.cold_inlet) annotation(
        Line(points = {{-64, 27}, {-42, 27}, {-42, 44}, {-16, 44}}, color = {0, 0, 255}));
      connect(sinkPressure3.flange, dynamicHX1.hot_outlet) annotation(
        Line(points = {{-72, -14}, {-2, -14}, {-2, -44}}, color = {0, 0, 255}));
      connect(sourceMassFlow3.flange, dynamicHX1.cold_inlet) annotation(
        Line(points = {{-72, -88}, {-39, -88}, {-39, -56}, {-2, -56}}, color = {0, 0, 255}));
      connect(sinkPressure.flange, dynamicHX1.cold_outlet) annotation(
        Line(points = {{68, -91}, {26, -91}, {26, -56}, {18, -56}}, color = {0, 0, 255}));
  connect(dynamicHX.cold_outlet, stodolaTurbine.inlet) annotation(
        Line(points = {{0, 42}, {78, 42}, {78, 20}}, color = {0, 0, 255}));
  connect(dynamicHX1.hot_inlet, stodolaTurbine.outlet) annotation(
        Line(points = {{14, -40}, {78, -40}, {78, -6}}, color = {0, 0, 255}));
  connect(stodolaTurbine.shaft_b, constantSpeed.flange) annotation(
        Line(points = {{66, -4}, {34, -4}, {34, 8}, {26, 8}}));
    end boiler_and_turbine_and_condenser_test;

    model pump_test
      inner ThermoPower.System system annotation(
        Placement(visible = true, transformation(origin = {84, 80}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      ThermoPower.Water.SourceMassFlow sourceMassFlow(h = 173.88e3, p0 = 0.08e5, w0 = 123) annotation(
        Placement(visible = true, transformation(origin = {-88, 12}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      ThermoPower.Water.SinkPressure sinkPressure(h = 183.405e3, p0 = 80e5) annotation(
        Placement(visible = true, transformation(origin = {74, 4}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Components.Pump Pump(P_in = 7999.999999999999, P_out = 7999999.999999999) annotation(
        Placement(visible = true, transformation(origin = {-6, -10}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    equation
      connect(sourceMassFlow.flange, Pump.inlet) annotation(
        Line(points = {{-78, 12}, {-14, 12}, {-14, -8}}, color = {0, 0, 255}));
  connect(Pump.outlet, sinkPressure.flange) annotation(
        Line(points = {{4, -10}, {64, -10}, {64, 4}}, color = {0, 0, 255}));
    protected
    end pump_test;

    model boiler_and_turbine_and_condenser_and_pump_test
      Components.DynamicHX dynamicHX1(L = 5, N = 100, Q_nom = 1, h_t = 100000, n = 50, use_q_nom = true, wc_nom = 123, wh_nom = 123) annotation(
        Placement(visible = true, transformation(origin = {8, -50}, extent = {{-10, 10}, {10, -10}}, rotation = 0)));
      ThermoPower.Water.SteamTurbineStodola steamTurbineStodola(Kt = 0.0067, eta_iso_nom = 0.85, phi(displayUnit = "rad"), pin(start = 80e5), pnom = 80e5, pout(start = 0.08e5), wnom = 123.5) annotation(
        Placement(visible = true, transformation(origin = {52, 10}, extent = {{-14, -14}, {14, 14}}, rotation = -90)));
      inner ThermoPower.System system annotation(
        Placement(visible = true, transformation(origin = {86, 86}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Mechanics.Rotational.Sources.ConstantSpeed constantSpeed(w_fixed = 3000 / 60 * 3.14159) annotation(
        Placement(visible = true, transformation(origin = {20, 8}, extent = {{-6, -6}, {6, 6}}, rotation = 0)));
      ThermoPower.Water.SourceMassFlow sourceMassFlow1(h = 3004772, use_T = false, w0 = 1350) annotation(
        Placement(visible = true, transformation(origin = {42, 74}, extent = {{8, -8}, {-8, 8}}, rotation = 0)));
      ThermoPower.Water.SourceMassFlow sourceMassFlow(h = 183.36e3, use_in_w0 = false, w0 = 123.58) annotation(
        Placement(visible = true, transformation(origin = {-73, 27}, extent = {{-9, -9}, {9, 9}}, rotation = 0)));
      Components.DynamicHX dynamicHX(N = 100, Q_nom = 1, h_t = 100000, n = 50, use_q_nom = true, wc_nom = 123, wh_nom = 123) annotation(
        Placement(visible = true, transformation(origin = {-6, 50}, extent = {{-10, 10}, {10, -10}}, rotation = 0)));
      ThermoPower.Water.SinkPressure sinkPressure(p0 = 1e5) annotation(
        Placement(visible = true, transformation(origin = {75, -91}, extent = {{-7, -7}, {7, 7}}, rotation = 0)));
      ThermoPower.Water.SinkPressure sinkPressure1(p0 = 80e5) annotation(
        Placement(visible = true, transformation(origin = {-72, 74}, extent = {{8, -8}, {-8, 8}}, rotation = 0)));
      ThermoPower.Water.SourceMassFlow sourceMassFlow3(T = 288.15, use_T = true, use_in_w0 = false, w0 = 2472) annotation(
        Placement(visible = true, transformation(origin = {-80, -88}, extent = {{-8, -8}, {8, 8}}, rotation = 0)));
      Components.Pump Pump(P_in = 7999.999999999999, P_out = 7999999.999999999) annotation(
        Placement(visible = true, transformation(origin = {-114, -20}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
      ThermoPower.Water.SinkPressure sinkPressure2(h = 183.405e3, p0 = 80e5) annotation(
        Placement(visible = true, transformation(origin = {-160, 34}, extent = {{-10, 10}, {10, -10}}, rotation = 180)));
    equation
      connect(sourceMassFlow.flange, dynamicHX.cold_inlet) annotation(
        Line(points = {{-64, 27}, {-42, 27}, {-42, 44}, {-16, 44}}, color = {0, 0, 255}));
      connect(sinkPressure1.flange, dynamicHX.hot_outlet) annotation(
        Line(points = {{-64, 74}, {-36, 74}, {-36, 56}, {-16, 56}}, color = {0, 0, 255}));
      connect(sourceMassFlow3.flange, dynamicHX1.cold_inlet) annotation(
        Line(points = {{-72, -88}, {-39, -88}, {-39, -56}, {-2, -56}}, color = {0, 0, 255}));
      connect(steamTurbineStodola.outlet, dynamicHX1.hot_inlet) annotation(
        Line(points = {{64, -2}, {64, -44}, {18, -44}}, color = {0, 0, 255}));
      connect(dynamicHX.cold_outlet, steamTurbineStodola.inlet) annotation(
        Line(points = {{4, 44}, {63, 44}, {63, 21}}, color = {0, 0, 255}));
      connect(sourceMassFlow1.flange, dynamicHX.hot_inlet) annotation(
        Line(points = {{34, 74}, {18, 74}, {18, 56}, {4, 56}}, color = {0, 0, 255}));
      connect(sinkPressure.flange, dynamicHX1.cold_outlet) annotation(
        Line(points = {{68, -91}, {26, -91}, {26, -56}, {18, -56}}, color = {0, 0, 255}));
      connect(constantSpeed.flange, steamTurbineStodola.shaft_b) annotation(
        Line(points = {{26, 8}, {26, 5}, {52, 5}, {52, 1}}));
      connect(Pump.outlet, sinkPressure2.flange) annotation(
        Line(points = {{-120, -13}, {-120, 6.5}, {-150, 6.5}, {-150, 34}}, color = {0, 127, 255}));
      connect(Pump.inlet, dynamicHX1.hot_outlet) annotation(
        Line(points = {{-106, -18}, {-2, -18}, {-2, -44}}, color = {0, 127, 255}));
    end boiler_and_turbine_and_condenser_and_pump_test;

    model pump_tretrfsf
      ThermoPower.Water.SourceMassFlow sourceMassFlow3(T = 288.15, use_T = true, use_in_w0 = false, w0 = 2472) annotation(
        Placement(visible = true, transformation(origin = {-80, -88}, extent = {{-8, -8}, {8, 8}}, rotation = 0)));
      ThermoPower.Water.SinkPressure sinkPressure(p0 = 1e5) annotation(
        Placement(visible = true, transformation(origin = {75, -91}, extent = {{-7, -7}, {7, 7}}, rotation = 0)));
      Components.DynamicHX dynamicHX(N = 100, Q_nom = 1, h_t = 100000, n = 50, use_q_nom = true, wc_nom = 123, wh_nom = 123) annotation(
        Placement(visible = true, transformation(origin = {-6, 50}, extent = {{-10, 10}, {10, -10}}, rotation = 0)));
      ThermoPower.Water.SourceMassFlow sourceMassFlow1(h = 3004772, use_T = false, w0 = 1350) annotation(
        Placement(visible = true, transformation(origin = {42, 74}, extent = {{8, -8}, {-8, 8}}, rotation = 0)));
      inner ThermoPower.System system annotation(
        Placement(visible = true, transformation(origin = {86, 86}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      ThermoPower.Water.SinkPressure sinkPressure1(p0 = 80e5) annotation(
        Placement(visible = true, transformation(origin = {-72, 74}, extent = {{8, -8}, {-8, 8}}, rotation = 0)));
      ThermoPower.Water.SinkPressure sinkPressure3(p0 = 0.08e5) annotation(
        Placement(visible = true, transformation(origin = {-80, -14}, extent = {{8, -8}, {-8, 8}}, rotation = 0)));
      ThermoPower.Water.SteamTurbineStodola steamTurbineStodola(Kt = 0.0067, eta_iso_nom = 0.85, phi(displayUnit = "rad"), pin(start = 80e5), pnom = 80e5, pout(start = 0.08e5), wnom = 123.5) annotation(
        Placement(visible = true, transformation(origin = {52, 10}, extent = {{-14, -14}, {14, 14}}, rotation = -90)));
      Components.DynamicHX dynamicHX1(L = 5, N = 100, Q_nom = 1, h_t = 100000, n = 50, use_q_nom = true, wc_nom = 123, wh_nom = 123) annotation(
        Placement(visible = true, transformation(origin = {8, -50}, extent = {{-10, 10}, {10, -10}}, rotation = 0)));
      Modelica.Mechanics.Rotational.Sources.ConstantSpeed constantSpeed(w_fixed = 3000 / 60 * 3.14159) annotation(
        Placement(visible = true, transformation(origin = {20, 8}, extent = {{-6, -6}, {6, 6}}, rotation = 0)));
      Components.Pump Pump(P_in = 7999.999999999999, P_out = 7999999.999999999) annotation(
        Placement(visible = true, transformation(origin = {-106, 26}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      ThermoPower.Water.SourceMassFlow sourceMassFlow2(h = 173.88e3, p0 = 0.08e5, w0 = 123) annotation(
        Placement(visible = true, transformation(origin = {-188, 48}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    equation
      connect(steamTurbineStodola.outlet, dynamicHX1.hot_inlet) annotation(
        Line(points = {{64, -2}, {64, -44}, {18, -44}}, color = {0, 0, 255}));
      connect(sinkPressure.flange, dynamicHX1.cold_outlet) annotation(
        Line(points = {{68, -91}, {26, -91}, {26, -56}, {18, -56}}, color = {0, 0, 255}));
      connect(sourceMassFlow3.flange, dynamicHX1.cold_inlet) annotation(
        Line(points = {{-72, -88}, {-39, -88}, {-39, -56}, {-2, -56}}, color = {0, 0, 255}));
      connect(sinkPressure1.flange, dynamicHX.hot_outlet) annotation(
        Line(points = {{-64, 74}, {-36, 74}, {-36, 56}, {-16, 56}}, color = {0, 0, 255}));
      connect(dynamicHX.cold_outlet, steamTurbineStodola.inlet) annotation(
        Line(points = {{4, 44}, {63, 44}, {63, 21}}, color = {0, 0, 255}));
      connect(sinkPressure3.flange, dynamicHX1.hot_outlet) annotation(
        Line(points = {{-72, -14}, {-2, -14}, {-2, -44}}, color = {0, 0, 255}));
      connect(constantSpeed.flange, steamTurbineStodola.shaft_b) annotation(
        Line(points = {{26, 8}, {26, 5}, {52, 5}, {52, 1}}));
      connect(sourceMassFlow1.flange, dynamicHX.hot_inlet) annotation(
        Line(points = {{34, 74}, {18, 74}, {18, 56}, {4, 56}}, color = {0, 0, 255}));
      connect(sourceMassFlow2.flange, Pump.inlet) annotation(
        Line(points = {{-178, 48}, {-114, 48}, {-114, 28}}, color = {0, 0, 255}));
      connect(Pump.outlet, dynamicHX.cold_inlet) annotation(
        Line(points = {{-100, 34}, {-16, 34}, {-16, 44}}, color = {0, 127, 255}));
    end pump_tretrfsf;
  equation

  end CombinationSteps;
  annotation(
    uses(ThermoPower(version = "3.1"), Modelica(version = "3.2.3")));
end Rankine;
