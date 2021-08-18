package ThermoPower_steps
  model test01
  
    ThermoSysPro.WaterSteam.HeatExchangers.DynamicWaterWaterExchanger
      echangeurAPlaques1D(
      modec=1,
      modef=1,
      N=5) annotation (Placement(transformation(extent={{-20,40},{0,60}},
            rotation=0)));
    ThermoSysPro.WaterSteam.BoundaryConditions.SourceP sourceP(
                                             T0=340)
      annotation (Placement(transformation(extent={{-80,40},{-60,60}}, rotation=0)));
    ThermoSysPro.WaterSteam.BoundaryConditions.SourceP sourceP1
                                              annotation (Placement(
          transformation(extent={{-60,20},{-40,40}}, rotation=0)));
    ThermoSysPro.WaterSteam.BoundaryConditions.SinkP puitsP
                                           annotation (Placement(transformation(
            extent={{40,40},{60,60}}, rotation=0)));
    ThermoSysPro.WaterSteam.BoundaryConditions.SinkP puitsP1
                                            annotation (Placement(transformation(
            extent={{20,20},{40,40}}, rotation=0)));
  ThermoSysPro.InstrumentationAndControl.Blocks.Sources.Rampe rampe(Duration = 400,Finalvalue = 5e5, Initialvalue = 3e5, Starttime = 300)  annotation(
      Placement(visible = true, transformation(origin = {-96, 50}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  equation
    connect(sourceP.C, echangeurAPlaques1D.Ec)
      annotation (Line(points={{-60,50},{-20,50}}, color={0,0,255}));
    connect(sourceP1.C, echangeurAPlaques1D.Ef) annotation (Line(points={{-40,30},
            {-15,30},{-15,44}}, color={0,0,255}));
    connect(echangeurAPlaques1D.Sc, puitsP.C) annotation (Line(points={{0,50},{20,
            50},{20,50},{40,50}},      color={0,0,255}));
    connect(echangeurAPlaques1D.Sf, puitsP1.C) annotation (Line(points={{-5,44},{
            -6,44},{-6,30},{20,30}}, color={0,0,255}));
  connect(rampe.y, sourceP.IPressure) annotation(
      Line(points = {{-84, 50}, {-74, 50}}, color = {0, 0, 255}));
    annotation (experiment(StopTime=1000), Diagram(graphics),
      Icon(graphics={
          Rectangle(
            lineColor={200,200,200},
            fillColor={248,248,248},
            fillPattern=FillPattern.HorizontalCylinder,
            extent={{-100.0,-100.0},{100.0,100.0}},
            radius=25.0),
          Rectangle(
            lineColor={128,128,128},
            extent={{-100.0,-100.0},{100.0,100.0}},
            radius=25.0),
          Polygon(
            origin={8.0,14.0},
            lineColor={78,138,73},
            fillColor={78,138,73},
            pattern=LinePattern.None,
            fillPattern=FillPattern.Solid,
            points={{-58.0,46.0},{42.0,-14.0},{-58.0,-74.0},{-58.0,46.0}})}),
      Documentation(info="<html>
  <p><b>Copyright &copy; EDF 2002 - 2019 </p>
  <p><b>ThermoSysPro Version 3.2 </h4>
  </html>"));
  
  end test01;

  model test02
    ThermoPower.Water.SteamTurbineStodola steamTurbine(Kt = 0.0104, PRstart = 1000, phi(displayUnit = "rad"), pin(fixed = false, start = 80e5), pnom = 80e5, pout(fixed = true, start = 0.08e5), wnom = 151.278, wstart = 55) annotation(
      Placement(visible = true, transformation(extent = {{-26, 26}, {24, 76}}, rotation = 0)));
    Modelica.Mechanics.Rotational.Sources.ConstantSpeed constantSpeed(phi(fixed = true, start = 0), w_fixed = 157) annotation(
      Placement(visible = true, transformation(origin = {38, 54}, extent = {{6, -6}, {-6, 6}}, rotation = 0)));
    ThermoPower.Water.SourcePressure sourcePressure(h = 2578e3, p0 = 80e5, use_in_p0 = true) annotation(
      Placement(visible = true, transformation(origin = {-66, 26}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Blocks.Sources.Step step(height = 20e5, offset = 60e5, startTime = 5) annotation(
      Placement(visible = true, transformation(origin = {-106, 52}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    inner ThermoPower.System system annotation(
      Placement(visible = true, transformation(origin = {-88, 84}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    ThermoPower.Examples.RankineCycle.Models.PrescribedPressureCondenser condenser(initOpt = ThermoPower.Choices.Init.Options.fixedState, p (displayUnit = "Pa") = 8000) annotation(
      Placement(visible = true, transformation(extent = {{58, -8}, {98, 32}}, rotation = 0)));
  equation
    connect(constantSpeed.flange, steamTurbine.shaft_b) annotation(
      Line(points = {{32, 54}, {25, 54}, {25, 51}, {15, 51}}));
    connect(sourcePressure.flange, steamTurbine.inlet) annotation(
      Line(points = {{-56, 26}, {-21, 26}, {-21, 71}}, color = {0, 0, 255}));
    connect(step.y, sourcePressure.in_p0) annotation(
      Line(points = {{-94, 52}, {-70, 52}, {-70, 34}}, color = {0, 0, 127}));
    connect(steamTurbine.outlet, condenser.steamIn) annotation(
      Line(points = {{19, 71}, {78, 71}, {78, 32}}, color = {0, 0, 255}));
    connect(condenser.steamIn, condenser.steamIn) annotation(
      Line(points = {{78, 32}, {78, 32}}));
  protected
  end test02;

  model step01
  ThermoPower.Water.SteamTurbineStodola steamTurbine(Kt = 0.0104, PRstart = 551.521, phi(displayUnit = "rad"), pin(start = 29.727e5), pnom = 3000000, pout(start = 0.0539e5), wnom = 56.3196, wstart = 56.3196) annotation(
      Placement(visible = true, transformation(origin = {-29, 25}, extent = {{-17, -17}, {17, 17}}, rotation = 0)));
  ThermoPower.Water.SourcePressure sourcePressure(h = 3.245e6, p0 = 29.727e5, use_in_p0 = false) annotation(
      Placement(visible = true, transformation(origin = {-72, 72}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  ThermoPower.Water.SinkPressure sinkPressure(h = 2.216e6, p0 = 0.0539e5) annotation(
      Placement(visible = true, transformation(origin = {31, 69}, extent = {{-11, -11}, {11, 11}}, rotation = 0)));
  Modelica.Mechanics.Rotational.Sources.ConstantSpeed constantSpeed(phi(fixed = true, start = 0), w_fixed = 157) annotation(
      Placement(visible = true, transformation(origin = {9, 17}, extent = {{5, -5}, {-5, 5}}, rotation = 0)));
  inner ThermoPower.System system annotation(
      Placement(visible = true, transformation(origin = {86, 88}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  equation
  connect(constantSpeed.flange, steamTurbine.shaft_b) annotation(
      Line(points = {{4, 17}, {-7, 17}, {-7, 25}, {-18, 25}}));
  connect(sourcePressure.flange, steamTurbine.inlet) annotation(
      Line(points = {{-62, 72}, {-43, 72}, {-43, 39}}, color = {0, 0, 255}));
  connect(steamTurbine.outlet, sinkPressure.flange) annotation(
      Line(points = {{-15, 39}, {12.5, 39}, {12.5, 69}, {20, 69}}, color = {0, 0, 255}));
  end step01;

  model step02
ThermoPower.Water.SteamTurbineStodola steamTurbine(Kt = 0.0104, PRstart = 551.521, phi(displayUnit = "rad"), pin(start = 29.727e5), pnom = 29.727e5, pout(start = 0.0539e5), wnom = 56.3196, wstart = 56.3196) annotation(
      Placement(visible = true, transformation(origin = {-25, 59}, extent = {{-17, -17}, {17, 17}}, rotation = 0)));
  ThermoPower.Water.SourcePressure sourcePressure(h = 3.2455e6, p0 = 29.727e5, use_in_p0 = false) annotation(
      Placement(visible = true, transformation(origin = {-72, 72}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Mechanics.Rotational.Sources.ConstantSpeed constantSpeed(phi(fixed = true, start = 0), w_fixed = 157) annotation(
      Placement(visible = true, transformation(origin = {3, 57}, extent = {{5, -5}, {-5, 5}}, rotation = 0)));
  inner ThermoPower.System system annotation(
      Placement(visible = true, transformation(origin = {86, 88}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  ThermoPower.Examples.RankineCycle.Models.PrescribedPressureCondenser condenser(initOpt = ThermoPower.Choices.Init.Options.fixedState, p(displayUnit = "Pa") = 0.0539e5) annotation(
      Placement(visible = true, transformation(origin = {48, 4}, extent = {{-16, -16}, {16, 16}}, rotation = 0)));
  equation
    connect(constantSpeed.flange, steamTurbine.shaft_b) annotation(
      Line(points = {{-2, 57}, {-8, 57}, {-8, 59}, {-14, 59}}));
    connect(sourcePressure.flange, steamTurbine.inlet) annotation(
      Line(points = {{-62, 72}, {-39, 72}, {-39, 73}}, color = {0, 0, 255}));
  connect(steamTurbine.outlet, condenser.steamIn) annotation(
      Line(points = {{-12, 72}, {48, 72}, {48, 20}}, color = {0, 0, 255}));
  end step02;

  model step03
ThermoPower.Water.SteamTurbineStodola steamTurbine(Kt = 0.0104, PRstart = 551.521, phi(displayUnit = "rad"), pin(start = 29.727e5), pnom = 29.727e5, pout(start = 0.0539e5), wnom = 56.3196, wstart = 56.3196) annotation(
      Placement(visible = true, transformation(origin = {-25, 59}, extent = {{-17, -17}, {17, 17}}, rotation = 0)));
  ThermoPower.Water.SourcePressure sourcePressure(h = 3.2455e6, p0 = 29.727e5, use_in_p0 = false) annotation(
      Placement(visible = true, transformation(origin = {-72, 72}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Mechanics.Rotational.Sources.ConstantSpeed constantSpeed(phi(fixed = true, start = 0), w_fixed = 157) annotation(
      Placement(visible = true, transformation(origin = {3, 57}, extent = {{5, -5}, {-5, 5}}, rotation = 0)));
  inner ThermoPower.System system annotation(
      Placement(visible = true, transformation(origin = {86, 88}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  ThermoPower.Examples.RankineCycle.Models.PrescribedPressureCondenser condenser(initOpt = ThermoPower.Choices.Init.Options.fixedState, p(displayUnit = "Pa") = 0.0539e5) annotation(
      Placement(visible = true, transformation(origin = {70, 18}, extent = {{-12, -12}, {12, 12}}, rotation = 0)));
  ThermoPower.Water.SinkPressure sinkPressure(h = 147163, p0 = 30.1462e5) annotation(
      Placement(visible = true, transformation(origin = {-45, -41}, extent = {{11, -11}, {-11, 11}}, rotation = 0)));
    ThermoPower.Examples.RankineCycle.Models.PrescribedSpeedPump prescribedSpeedPump(
      
      redeclare package FluidMedium = Modelica.Media.Water.StandardWater,
      head_nom={450,300,0}, hstart = 147163,n0=1500,
      nominalInletPressure(displayUnit = "Pa") = 5390,
      nominalMassFlowRate= 56.32,
      nominalOutletPressure(displayUnit = "Pa") =3000000,
      q_nom={0,0.055,0.1},
      rho0(displayUnit = "kg/m3") =1000) annotation (Placement(visible = true, transformation(origin = {20, -36}, extent = {{16, -16}, {-16, 16}}, rotation = 0)));
    Modelica.Blocks.Sources.Ramp gasTemperature(duration = 0, height = 0, offset = 1500) annotation(
      Placement(visible = true, transformation(origin = {27, -5}, extent = {{-7, -7}, {7, 7}}, rotation = 0)));
  equation
    connect(constantSpeed.flange, steamTurbine.shaft_b) annotation(
      Line(points = {{-2, 57}, {-8, 57}, {-8, 59}, {-14, 59}}));
    connect(sourcePressure.flange, steamTurbine.inlet) annotation(
      Line(points = {{-62, 72}, {-39, 72}, {-39, 73}}, color = {0, 0, 255}));
  connect(steamTurbine.outlet, condenser.steamIn) annotation(
      Line(points = {{-12, 72}, {70, 72}, {70, 30}}, color = {0, 0, 255}));
    connect(prescribedSpeedPump.inlet, condenser.waterOut) annotation(
      Line(points = {{36, -36}, {70, -36}, {70, 6}}, color = {0, 0, 255}));
  connect(sinkPressure.flange, prescribedSpeedPump.outlet) annotation(
      Line(points = {{-34, -40}, {4, -40}, {4, -36}}, color = {0, 0, 255}));
  connect(gasTemperature.y, prescribedSpeedPump.nPump) annotation(
      Line(points = {{34, -6}, {48, -6}, {48, -26}, {32, -26}}, color = {0, 0, 127}));
  end step03;

  model step04
ThermoPower.Water.SteamTurbineStodola steamTurbine(Kt = 0.0104, PRstart = 551.521, phi(displayUnit = "rad"), pin(start = 29.727e5), pnom = 29.727e5, pout(start = 0.0539e5), wnom = 56.3196, wstart = 56.3196) annotation(
      Placement(visible = true, transformation(origin = {-25, 59}, extent = {{-17, -17}, {17, 17}}, rotation = 0)));
  ThermoPower.Water.SourcePressure sourcePressure(h = 3.2455e6, p0 = 29.727e5, use_in_p0 = false) annotation(
      Placement(visible = true, transformation(origin = {-72, 72}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Mechanics.Rotational.Sources.ConstantSpeed constantSpeed(phi(fixed = true, start = 0), w_fixed = 157) annotation(
      Placement(visible = true, transformation(origin = {3, 57}, extent = {{5, -5}, {-5, 5}}, rotation = 0)));
  inner ThermoPower.System system annotation(
      Placement(visible = true, transformation(origin = {86, 88}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  ThermoPower.Examples.RankineCycle.Models.PrescribedPressureCondenser condenser(initOpt = ThermoPower.Choices.Init.Options.fixedState, p(displayUnit = "Pa") = 0.0539e5) annotation(
      Placement(visible = true, transformation(origin = {70, 18}, extent = {{-12, -12}, {12, 12}}, rotation = 0)));
  ThermoPower.Water.SinkPressure sinkPressure(h = 147163, p0 = 30.1462e5) annotation(
      Placement(visible = true, transformation(origin = {-45, -41}, extent = {{11, -11}, {-11, 11}}, rotation = 0)));
    ThermoPower.Examples.RankineCycle.Models.PrescribedSpeedPump prescribedSpeedPump(
      
      redeclare package FluidMedium = Modelica.Media.Water.StandardWater,
      head_nom={450,300,0}, hstart = 147163,n0=1500,
      nominalInletPressure(displayUnit = "Pa") = 5390,
      nominalMassFlowRate= 56.32,
      nominalOutletPressure(displayUnit = "Pa") =3000000,
      q_nom={0,0.055,0.1},
      rho0(displayUnit = "kg/m3") =1000) annotation (Placement(visible = true, transformation(origin = {20, -36}, extent = {{16, -16}, {-16, 16}}, rotation = 0)));
    Modelica.Blocks.Sources.Ramp gasTemperature(duration = 0, height = 0, offset = 1500) annotation(
      Placement(visible = true, transformation(origin = {29, 5}, extent = {{-7, -7}, {7, 7}}, rotation = 0)));
  equation
    connect(constantSpeed.flange, steamTurbine.shaft_b) annotation(
      Line(points = {{-2, 57}, {-8, 57}, {-8, 59}, {-14, 59}}));
    connect(sourcePressure.flange, steamTurbine.inlet) annotation(
      Line(points = {{-62, 72}, {-39, 72}, {-39, 73}}, color = {0, 0, 255}));
  connect(steamTurbine.outlet, condenser.steamIn) annotation(
      Line(points = {{-12, 72}, {70, 72}, {70, 30}}, color = {0, 0, 255}));
    connect(prescribedSpeedPump.inlet, condenser.waterOut) annotation(
      Line(points = {{36, -36}, {70, -36}, {70, 6}}, color = {0, 0, 255}));
  connect(sinkPressure.flange, prescribedSpeedPump.outlet) annotation(
      Line(points = {{-34, -40}, {4, -40}, {4, -36}}, color = {0, 0, 255}));
  connect(gasTemperature.y, prescribedSpeedPump.nPump) annotation(
      Line(points = {{37, 5}, {48, 5}, {48, -26}, {32, -26}}, color = {0, 0, 127}));
  
  end step04;
  annotation(
    uses(ThermoPower(version = "3.1"), ThermoSysPro(version = "3.2"), Modelica(version = "3.2.3")));
end ThermoPower_steps;
