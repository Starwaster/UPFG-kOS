function list_of {
  parameter i,n.
  //parameter n.
  local out is list().
  from {local x is -1.}
  until x = n
  step { set x to x+1. }
  do {
    out:add(i).
  }
  return out.
}

function getAngleFromFrame {
  parameter vectors.
  parameter frame.
  parameter type.

  set vectors to vectors:normalized.
  if type = "pitch" {
    local angle is 90 - vang(frame[0],vectors).
    return angle.
  }
  else if type = "yaw" {
    local angle is arctan2(vdot(frame[2],vectors),vdot(frame[0],vectors)).
    if angle < 0 {
      return 360+angle.
    }
    else return angle.
  }
}

function safeAcosd {
  parameter angle.

  set angle to min(angle, 1).
  set angle to max(angle, -1).
  return arccos(angle).
}

function swapYZ {
  parameter vectors.
  return v(vectors:x,vectors:z,vectors:y).
}

function update {
  parameter state.
  clearscreen.
  local ve is 0.
  local fT is 0.
  list engines in eng.
  for e in eng {
    if e:ignition {
      set fT to fT + e:thrust.
      set ve to e:isp * g0.
    }
  }
  set state["thrust"] to fT.
  set state["ve"] to ve.
  set state["acc"] to fT/ship:mass.
  set state["time"] to time:seconds - t0.
  set state["mass"] to ship:mass.
  set state["radius"] to getECIVecs(ship:orbit)[0].
  set state["velocity"] to getECIVecs(ship:orbit)[1].
  return state.
}

function rodrigues {
  parameter vectors.
  parameter axis.
  parameter angle.

  set axis to axis:normalized.
  local rotated is vectors*cos(angle).
  set rotated to rotated + vcrs(axis,vectors)*sin(angle).
  set rotated to rotated + axis*vdot(axis,vectors)*(1-cos(angle)).
  return rotated.
}

function convergeUPFG {
  parameter vehicle is list().
  parameter trgt is lexicon().
  parameter state is lexicon().
  parameter internal is lexicon().
  parameter times is 0.
  parameter maxIters is 50.

  local conv is 0.05.
  local fail is 1.
  local out is upfg(vehicle,trgt,state,internal).
  set internal to out[1].
  local guidance is out[0].

  local i is 1.
  until i = maxIters
  {
    local t1 is internal["tgo"].
    set out to upfg(vehicle,trgt,state,internal).
    set internal to out[1].
    set guidance to out[0].
    local t2 is internal["tgo"].
    if abs((t1-t2)/t1) < conv {
      if times > 0 {
        print "UPFG converged after : " + i + " iterations.".
      }
      set fail to 0.
      return list(guidance,internal,1).
    }
    set i to i+1.
  }
  return list(guidance,internal,0).
}

FUNCTION east {
	RETURN VCRS(SHIP:UP:VECTOR, SHIP:NORTH:VECTOR).
}

//	navPitch	gets current compass pitch (degrees above horizon)
FUNCTION navPitchV {
  DECLARE PARAMETER pointing.
	RETURN 90-VANG(SHIP:UP:VECTOR, pointing).
}

//	navHdgV		gets compass heading for any vector
FUNCTION navHdgV {
	DECLARE PARAMETER pointing.
	LOCAL result IS ARCTAN2( VDOT(east(), pointing), VDOT(SHIP:NORTH:VECTOR, pointing) ).
	IF result < 0 { RETURN result+360. } ELSE RETURN result.
}

//	navHeading	shortcut to get current compass heading (degrees from north towards east)
FUNCTION navHeading {
	RETURN navHdgV(SHIP:FACING:FOREVECTOR).
}

//	navRoll		gets current roll angle
FUNCTION navRoll {
	IF VANG(SHIP:FACING:VECTOR, SHIP:UP:VECTOR) < 0.2 { RETURN 0. } //	deadzone against gimbal lock (when vehicle is too vertical, roll angle becomes indeterminate)
	ELSE {
		LOCAL raw IS VANG(VXCL(SHIP:FACING:VECTOR, SHIP:UP:VECTOR), SHIP:FACING:STARVECTOR).
		IF VANG(SHIP:UP:VECTOR, SHIP:FACING:TOPVECTOR) > 90 {
			RETURN 270-raw.
		} ELSE {
			RETURN -90-raw.
		}
	}.
}

//liborbitalstate

function eciVecsToKepElem {
  parameter mu, r, v.

  // Semi-major Axis
  local a to 1 / ((2 / r:mag) - (v:sqrmagnitude / mu)).

  // Eccentricity
  local h to vcrs(r, v).
  local eVec to vcrs(v, h) / mu - r:normalized.
  local e to eVec:mag.

  // Inclination
  local i to arccos(h:z / h:mag).

  // Longitude of Ascending Node
  local n to V(-h:y, h:x, 0).
  local lan to 0.
  if n:mag <> 0 {
    set lan to arccos(n:x / n:mag).
    if n:y < 0 { set lan to 360 - lan. }
  }

  // Argument of Periapsis
  local aop to 0.
  if n:mag <> 0 and e <> 0 {
    set aop to arccos(vdot(n, eVec) / (n:mag * e)).
    if eVec:z < 0 { set aop to 360 - aop. }
  }

  // True Anomaly
  local ta to 0.
  if e <> 0 { // TODO: do something reasonable when the orbit is circular
    set ta to arccos(vdot(eVec, r) / (e * r:mag)).
    if vdot(r, v) < 0 { set ta to 360 - ta. }
  }

//  // Eccentric Anomaly
//  local ea to arccos((e + cos(ta)) / (1 + e * cos(ta))).
//  if 180 < ta and ta < 360 { set ea to 360 - ea. }
//
//  // Mean Anomaly
//  local ma to ea - e * sin(ea).

  // TODO: change ta to ma?
  return list(a, e, i, lan, aop, ta).
}

function kepElemToEciVecs {
  // Best reference I've found:
  // http://ccar.colorado.edu/ASEN5070/handouts/kep2cart_2002.doc
  parameter mu, elems.

  local a to elems[0].
  local e to elems[1].
  local i to elems[2].
  local lan to elems[3].
  local aop to elems[4].
  local ta to elems[5].

  local p to a * (1 - e^2).
  local r to p / (1 + e * cos(ta)).
  local h to sqrt(mu * p).

  local coslan to cos(lan).
  local sinlan to sin(lan).
  local cosaopta to cos(aop + ta).
  local sinaopta to sin(aop + ta).
  local cosi to cos(i).
  local sini to sin(i).

  local x to r * (coslan * cosaopta - sinlan * sinaopta * cosi).
  local y to r * (sinlan * cosaopta + coslan * sinaopta * cosi).
  local z to r * (sini * sinaopta).

  local sinta to sin(ta).
  local herpsinta to h * e * sinta / (r * p).
  local hr to h / r.

  local vx to x * herpsinta - hr * (coslan * sinaopta + sinlan * cosaopta * cosi).
  local vy to y * herpsinta - hr * (sinlan * sinaopta - coslan * cosaopta * cosi).
  local vz to z * herpsinta + hr * sini * cosaopta.

  return list(V(x, y, z), V(vx, vy, vz)).
}

function getECIVecs
{
	parameter p_obt.
	return kepElemToEciVecs(p_obt:body:mu, list(p_obt:semimajoraxis,
	                                            p_obt:eccentricity,
									      	    p_obt:inclination,
										          p_obt:longitudeofascendingnode,
      										    p_obt:argumentofperiapsis,
	      									    p_obt:trueanomaly)).
}
