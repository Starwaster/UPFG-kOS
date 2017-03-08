function peg {
  parameter cycle is 0.1.
  parameter trgt is lexicon().//cycle,alti,vt,vr,tgt,acc,ve,oldA,oldB,oldT.
  parameter state is lexicon().
  parameter previous is lexicon().

  set alti to ship:obt:body:distance.
  set vt to sqrt(ship:velocity:orbit:sqrmagnitude - ship:verticalspeed^2).
  set vr to ship:verticalspeed.
  set tgt to trgt["radius"].
  set acc to state["acc"].
  set ve to state["ve"].
  set oldA to previous["A"].
  set oldB to previous["B"].
  set oldT to previous["T"].

  set tau to ve/acc.
  if oldA = 0 and oldB = 0 {
    set b0 to ve*ln(tau / (tau - oldT)).
    set b1 to b0*tau - ve*oldT.
    set c0 to b0*oldT - b1.
    set c1 to c0*tau - ve*oldT^2/2.

    set z0 to -vr.
    set z1 to tgt - alti - vr*oldT.

    set oldB to (z1/c0 - z0/b0) / (c1/c0 - b1/b0).
    set oldA to (z0 - b1*oldB)/b0.
  }
  set angM to vcrs(v(alti,0,0),v(vr,vt,0)):mag.
  set tgtV to trgt["velocity"].
  set tgtM to vcrs(v(tgt,0,0),v(0,tgtV,0)):mag.
  set dMom to tgtM - angM.

  set C to (mu/tgt^2 - tgtV^2/tgt) / (acc / (1- oldt/tau)).
  set f_r_T to oldA + oldB*oldT + C.
  set C to (mu/alti^2 - vt^2/alti) / acc.
  set f_r to oldA + C.
  set f_r_dot to (f_r_T - f_r)/oldT.

  set f_theta to 1 - f_r^2/2.
  set f_theta_dot to -f_r*f_r_dot.
  set f_theta_2dot to -f_r_dot^2/2.

  set avgR to (alti + tgt)/2.
  set deltaV to dMom/avgR + ve*(oldT - cycle)*(f_theta_dot+f_theta_2dot*tau) + f_theta_2dot*ve*(oldT-cycle)^2/2.
  set deltaV to deltaV/(f_theta + f_theta_dot*tau + f_theta_2dot*tau^2).

  set T to tau*(1- constant:e^(-deltaV/ve)).

  if T > 5 {
    set b0 to ve*ln(tau / (tau - T)).
    set b1 to b0*tau - ve*T.
    set c0 to b0*T - b1.
    set c1 to c0*tau - ve*T^2/2.

    set z0 to -vr.
    set z1 to tgt - alti - vr*T.

    set B to (z1/c0 - z0/b0) / (c1/c0 - b1/b0).
    set A to (z0 - b1*B)/b0.
  }
  else {
    set A to oldA.
    set B to oldB.
  }
  set out to lexicon().
  set out["A"] to A.
  set out["B"] to B.
  set out["C"] to C.
  set out["T"] to T.
  set out["Vgo"] to b0.

  return out.
}
