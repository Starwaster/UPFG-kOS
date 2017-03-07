parameter apo.
parameter per.
parameter inc.
parameter lan_d.
parameter slip.

//runoncepath("newupfg").
runoncepath("newpeg").
runoncepath("trgt_app").
runoncepath("vehicle_build").
//runoncepath("CSEroutine").
runoncepath("inst_az").
runoncepath("misc_lib").

global g0 is 9.80655.
global bigR is body:radius.
global mu is body:mu.
local vehicle is falcon9().
local trgt is launchTargeting(ship:latitude,ship:longitude,per,apo,inc,lan_d,slip).
global t0 is time:seconds + trgt["eta"].
//global t is time:seconds - t0.
local state is lexicon().
global g_pitch is 90.
global g_yaw is inst_az(inc,trgt["velocity"]).
lock steering to heading(g_yaw,g_pitch).
global throt is 0.
lock throttle to throt.
local kick_pitch is 0.
SET WARPMODE TO "RAILS".
warpto(t0 - 10).
//Countdown
local ignition is 0.
//sas on.
until 0 {
  set state to update(state).
  print "Time since launch : " + round(state["time"],2) + "s".
  if state["time"] > -3 and ignition = 0{
    set throt to 1.
    stage.
    set ignition to 1.
  }
  if state["time"] > 0 and ignition = 1 and state["acc"] > g0 {
    stage.
    set ignition to 2.
    set t0 to time:seconds.
  }
  if ship:velocity:surface:mag>60 {
    set kick_pitch to 90 - ship:velocity:surface:mag/state["time"].
    break.
  }
  wait 0.2.
}
//sas off.

//Open-loop control
until 0 {
  set state to update(state).
  print "Time               : " + round(state["time"],2) + " s".
  print "Target Apoapsis    : " + apo + " km".
  print "Current Apoapsis   : " + round(ship:apoapsis/1000,2) + " km".
  print "Target Periapsis   : " + per + " km".
  print "Current Periapsis  : " + round(ship:Periapsis/1000,2) + " km".
  print "Target Inclination : " + inc + " degrees".
  print "Inclination        : " + round(ship:orbit:inclination,2) +" degrees".
  print "Target LAN         : " + lan_d + " degrees".
  print "LAN                : " + round(ship:obt:lan,2) +" degrees".
  print "Heading            : " + round(g_yaw,2) + " , " + round(g_pitch,2).
  print "Current Thrust     : " + round(state["thrust"],2) + " kN".
  local turn_pitch is 90 - vang(ship:up:vector,ship:srfprograde:vector).
  if g_pitch > kick_pitch {
    set g_pitch to g_pitch - 0.2.
  }
  else if turn_pitch < g_pitch {
    set g_pitch to turn_pitch.
  }
  if ship:maxthrust < 1 {
    stage.
    break.
  }
  wait 0.2.
}

print "Staging ! Ullage for 3 Seconds".
rcs on.
set ship:control:fore to 0.9.
wait 3.
stage.
print "Ignition ! Spool up for 5 seconds".
until 0 {
  wait 5.
  break.
}
set ship:control:fore to 0.

//PEG-loop
local tc is 0.1.
local good is false.
local guided is lexicon("A",0,"B",0,"T",vehicle["maxT"]).
local tcall is state["time"].
until 0 {
  set state to update(state).
  if state["time"]-tcall > tc {
    set guided to peg(tc,trgt,state,guided).
    local t1 is guided["T"].
    set guided to peg(tc,trgt,state,guided).
    local t2 is guided["T"].
    if abs(t2-t1)/t1 < 0.05 {
      set A to guided["A"].
      set B to guided["B"].
      set C to guided["C"].
    }
    set tcall to state["time"].
  }
  set g_pitch to A - B*(state["time"] - tcall) + C.
  set g_pitch to arcsin(min(1,max(-1,g_pitch))).
  //set g_yaw to inst_az(inc,trgt["velocity"]).
  if guided["T"] < tc {
    break.
  }
  if state["velocity"] > trgt["velocity"] {
    break.
  }
  print "Time               : " + round(state["time"],2) + " s".
  print "Target Apoapsis    : " + apo + " km".
  print "Current Apoapsis   : " + round(ship:apoapsis/1000,2) + " km".
  print "Target Periapsis   : " + per + " km".
  print "Current Periapsis  : " + round(ship:Periapsis/1000,2) + " km".
  print "Target Inclination : " + inc + " degrees".
  print "Inclination        : " + round(ship:orbit:inclination,2) +" degrees".
  print "Target LAN         : " + lan_d + " degrees".
  print "LAN                : " + round(ship:obt:lan,2) +" degrees".
  print "Heading            : " + round(g_yaw,2) + " , " + round(g_pitch,2).
  print "Current Thrust     : " + round(state["thrust"],2) + " kN".
  print "Velocity to go     : " + round(guided["vgo"],2) + " m/s".
  print "Time to go         : " + round(guided["T"],2) + " s".
  wait 0.
}

unlock steering.
unlock throttle.
set ship:control:pilotmainthrottle to 0.
