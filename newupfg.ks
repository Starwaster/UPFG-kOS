function upfg {
  parameter vehicle is lexicon().
  parameter trgt is lexicon().
  parameter state is lexicon().
  parameter previous is lexicon().

  local gamma is trgt["angle"].
  local iy is trgt["normal"].
  local rdval is trgt["radius"].
  local vdval is trgt["velocity"].
  local t_ is state["time"].
  local m_ is state["mass"].
  local r_ is state["radius"].
  local v_ is state["velocity"].
  local cser is previous["cser"].
  local rbias is previous["rbias"].
  local rd is previous["rd"].
  local rgrav is previous["rgrav"].
  local tp is previous["time"].
  local vprev is previous["v"].
  local vgo is previous["vgo"].
  set my_proc to processor("second").
  set my_conn to my_proc:connection.

  //BLOCK 1
  //local n_ is vehicle:length - 1.
  local SM is 1.
  local aL is 0.
  local md is 0.
  local ve is 0.
  local fT is 0.
  local acc is 0.
  local tau is 0.
  local tb is 0.

  list engines in eng.
  for e in eng {
    if e:ignition {
      set fT to fT + e:thrust.
      set ve to e:isp * g0.
    }
  }
  set md to fT/ve.
  set acc to fT/m_.
  set tau to ve/acc.
  set tb to vehicle["maxT"].

  //BLOCK 2
  local dt is t_ - tp.
  local dvsensed is v_ - vprev.
  set vgo to vgo - dvsensed.
  local vgo1 is vgo.
  set tb to tb - previous["tb"].

  //BLOCK 3
  local L1 is vgo:mag.
  set tgo to tau * (1 - constant:e^(-L1/ve)).

  //BLOCK 4
  local L_ is L1.                         //b0
  local J_ is L_*tau - ve * tgo.          //b1
  local S_ is L_*tgo - J_.                //c0
  local Q_ is S_*tau - ve * tgo^2 / 2.    //c1
  local P_ is Q_*tau - ve * tgo^3 / 6.
  local H_ is J_*tgo - Q_.

  //BLOCK 5
  local lambda is vgo:normalized.
  local rgrav1 is rgrav.
  if previous["tgo"] <> 0 {
    set rgrav to (tgo/previous["tgo"])^2 * rgrav.
  }
  local rgo is rd - (r_ + v_*tgo + rgrav). local rgo1 is rgo.
  local iz is vcrs(rd,iy).
  set iz to iz:normalized.
  local iz1 is iz.
  local rgoxy is rgo - vdot(iz,rgo)*iz.
  local rgoz is (S_ - vdot(lambda,rgoxy)) / vdot(lambda,iz).
  set rgo to rgoxy + rgoz*iz  + rbias.
  local lambdade is Q_ - S_*J_/L_.
  local lambdadot is (rgo - S_ * lambda) / lambdade.
  local uF is lambda - lambdadot*J_/L_.
  set uF to uF:normalized.
  local phi is arccos(vdot(uF,lambda)) * constant:degtorad.
  local phidot is -phi*L_/J_.
  local vthrust is (L_ - 0.5*L_*phi^2 - J_*phi*phidot - 0.5*H_*phidot^2)*lambda.
  set vthrust to vthrust - (L_*phi + J_*phidot)*lambdadot:normalized.
  local rthrust is (S_ - 0.5*S_*phi^2 - Q_*phi*phidot - 0.5*P_*phidot^2)*lambda.
  set rthrust to rthrust - (S_*phi + Q_*phidot)*lambdadot:normalized.
  local vbias is vgo - vthrust.
  set rbias to rgo - rthrust.

  //BLOCK 6
  local up_ is r_:normalized.
  local east_ is vcrs(v(0,0,1),up_).
  set east_ to east_:normalized.
  local frame is list(up_ , v(0,0,1), east_).
  local pitch_ is getAngleFromFrame(uf, frame, "pitch").
  local yaw_ is getAngleFromFrame(uf, frame, "yaw").

  //BLOCK 7
  local rc1 is r_ - 0.1*rthrust - vthrust*tgo/30.
  local vc1 is v_ + 1.2*rthrust/tgo - 0.1*vthrust.
  //local send_me is list(rc1,vc1,tgo,cser).
  //my_conn:sendmessage(send_me).
  //wait until core:messages:length > 0.
  //local CSEOut is core:messages:pop:content.
  local CSEOut is CSEroutine(rc1,vc1,tgo,cser).
  local rc2 is CSEOut["r"].
  local vc2 is CSEOut["v"].
  set cser to CSEOut["last"].
  //local rc2 is swapYZ(positionat(ship,time:seconds + tgo) - body:position).
  //local vc2 is swapYZ(velocityat(ship,time:seconds + tgo):orbit).
  local vgrav is vc2 - vc1.
  set rgrav to rc2 - rc1 - vc1*tgo.

  //BLOCK 8
  local rp is r_ + v_*tgo + rthrust + rgrav.
  set rp to rp - vdot(rp,iy)*iy.
  set rd to rdval * rp:normalized.
  local ix is rd:normalized.
  set iz to  vcrs(ix,iy).
  local vd is vdval * iz.
  local vgop is vd - vgrav + vbias.
  local dvgo is 0.0*(vgop - vgo).
  set vgo to  vgop + dvgo.

  //Returns
  local current is previous.
  set current["cser"] to cser.
  set current["rbias"] to rbias.
  set current["rd"] to rd.
  set current["rgrav"] to rgrav.
  set current["tb"] to current["tb"] + dt.
  set current["time"] to t_.
  set current["tgo"] to tgo.
  set current["v"] to v_.
  set current["vgo"] to vgo.

  local guidance is lexicon().
  set guidance["pitch"] to pitch_.
  set guidance["yaw"] to yaw_.
  set guidance["tgo"] to tgo.

  local debug is lexicon().
  set debug["dvsensed"] to dvsensed.
  set debug["vgo1"] to vgo1.
  set debug["L1"] to L1.
  set debug["tgo"] to tgo.
  set debug["L"] to L_.
  set debug["J"] to J_.
  set debug["S"] to S_.
  set debug["Q"] to Q_.
  set debug["P"] to P_.
  set debug["H"] to H_.
  set debug["lambda"] to lambda.
  set debug["rgrav1"] to rgrav1.
  set debug["rgo1"] to rgo1.
  set debug["iz1"] to iz1.
  set debug["rgoxy"] to rgoxy.
  set debug["rgoz"] to rgoz.
  set debug["rgo2"] to rgo.
  set debug["lambdade"] to lambdade.
  set debug["lambdadot"] to lambdadot.
  set debug["uF"] to uF.
  set debug["phi"] to phi.
  set debug["phidot"] to phidot.
  set debug["vthrust"] to vthrust.
  set debug["rthrust"] to rthrust.
  set debug["vbias"] to vbias.
  set debug["rbias"] to rbias.
  set debug["guidance"] to guidance.
  set debug["rc1"] to rc1.
  set debug["vc1"] to vc1.
  set debug["rc2"] to rc2.
  set debug["vc2"] to vc2.
  set debug["cser"] to cser.
  set debug["vgrav"] to vgrav.
  set debug["rgrav2"] to rgrav.
  set debug["rp"] to rp.
  set debug["rd"] to rd.
  set debug["ix"] to ix.
  set debug["iz2"] to iz.
  set debug["vd"] to vd.
  set debug["vgop"] to vgop.
  set debug["dvgo"] to dvgo.
  set debug["vgo2"] to vgo.
  log debug + "Time : " + t_ to debug.txt.

  return list(guidance,current,debug).
}
