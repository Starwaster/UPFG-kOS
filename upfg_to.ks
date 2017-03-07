parameter apo.
parameter per.
parameter inc.
parameter lan_d.
parameter slip.

runoncepath("newupfg").
//runoncepath("newpeg").
runoncepath("trgt_app").
runoncepath("vehicle_build").
runoncepath("CSEroutine").
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
global g_yaw is trgt["azimuth"].
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

//UPFG loop (closed-loop)
local good is 0.
local upfg_debug is lexicon().
local upfg_internal is lexicon().
local upfg_guidance is lexicon().
local out is lexicon().
local tc is 0.1.
until 0 {
  set state to update(state).
  local cser is lexicon().
  if good = 0 {
    set cser["dtcp"] to 0.
    set cser["xcp"] to 0.
    set cser["A"] to 0.
    set cser["D"] to 0.
    set cser["E"] to 0.
    local rdinit is rodrigues(state["radius"]:normalized, -trgt["normal"], 20).
    set rdinit to rdinit * trgt["radius"].
    local vdinit is trgt["velocity"]*vcrs(-trgt["normal"], rdinit):normalized.
    set vdinit to vdinit - state["velocity"].
    set upfg_internal["cser"] to cser.
    set upfg_internal["rbias"] to v(0,0,0).
    set upfg_internal["rd"] to rdinit.
    set upfg_internal["rgrav"] to -(mu/2)*state["radius"]/state["radius"]:mag^3.
    set upfg_internal["tb"] to 0.
    set upfg_internal["time"] to state["time"].
    set upfg_internal["tgo"] to vehicle["maxT"].
    set upfg_internal["v"] to state["velocity"].
    set upfg_internal["vgo"] to vdinit.
    set out to convergeUPFG(vehicle,trgt,state,upfg_internal,state["time"],50).
    if out[2] = 1{
      set upfg_internal to out[1].
      set upfg_guidance to out[0].
      set good to 1.
    }
  }
  if (state["time"] - upfg_internal["time"]) > tc {
    set out to upfg(vehicle,trgt,state,upfg_internal).
    set upfg_internal to out[1].
    set upfg_guidance to out[0].
    set upfg_debug to out[2].
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
    print "Velocity to go     : " + round(upfg_internal["vgo"]:mag,2) + " m/s".
    print "Time to go         : " + round(upfg_internal["tgo"],2) + " s".
  }
  if upfg_guidance["tgo"] < 0.05 {
    break.
  }
  if state["velocity"]:mag >= trgt["velocity"] {
    break.
  }
  set g_yaw to  upfg_guidance["yaw"].
  set g_pitch to upfg_guidance["pitch"].
  wait 0.2.
}


unlock steering.
unlock throttle.
set ship:control:pilotmainthrottle to 0.
