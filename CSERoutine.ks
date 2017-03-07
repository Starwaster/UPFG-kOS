function CSEroutine {
  parameter r0.
  parameter v0.
  parameter dt.
  parameter last.

  if last["dtcp"] = 0 {
    set dtcp to dt.
  }
  else {
    set dtcp to last["dtcp"].
  }
  set xcp to last["xcp"].
  set x to xcp.
  set A to last["A"].
  set D to last["D"].
  set E to last["E"].

  set kmax to 10.
  set imax to 10.

  //PLATE 5-2
  if dt >= 0 {
    set f0 to 1.
  }
  else {
    set f0 to -1.
  }

  set n to 0.
  set r0m to r0:mag.

  set f1 to f0*sqrt(r0m/mu).
  set f2 to 1/f1.
  set f3 to f2/r0m.
  set f4 to f1*r0m.
  set f5 to f0/sqrt(r0m).
  set f6 to f0*sqrt(r0m).

  set ir0 to r0/r0m.
  set v0s to f1*v0.
  set sigma0s to vdot(ir0,v0s).
  set b0 to vdot(v0s,v0s) - 1.
  set alphas to 1-b0.

  //PLATE 5-3
  set xguess to f5*x.
  set xlast to  f5*xcp.
  set xmin to 0.
  set dts to f3*dt.
  set dtlast to f3*dtcp.
  set dtmin to 0.

  set xmax to 2*constant:pi/sqrt(abs(alphas)).

  if alphas > 0 {
    set dtmax to xmax/alphas.
    set xP to xmax.
    set Ps to dtmax.
    until dts < Ps {
      set n to n+1.
      set dts to dts-Ps.
      set dtlast to dtlast-Ps.
      set xguess to xguess-xP.
      set xlast to xlast-xP.
    }
  }
  else {
    set dtmax to KTTI(xmax,sigma0s,alphas,kmax)["t"].
    if dtmax < dts {
      until dtmax < dts {
        set dtmin to dtmax.
        set xmin to xmax.
        set xmax to 2*xmax.
        set dtmax to KTTI(xmax,sigma0s,alphas,kmax)["t"].
      }
    }
  }

  //PLATE 5-4
  if xmin >= xguess or xguess >=xmax {
    set xguess to 0.5*(xmin+xmax).
  }

  set dtguess to KTTI(xguess,sigma0s,alphas,kmax)["t"].

  if dts < dtguess {
    if xguess<xlast and xlast<xmax and dtguess<dtlast and dtlast<dtmax {
      set xmax to xlast.
      set dtmax to dtlast.
    }
  }
  else {
    if xmin<xlast and xlast<xguess and dtmin<dtlast and dtlast<dtguess {
      set xmin to xlast.
      set dtmin to dtlast.
    }
  }

  set KILout to KIL(imax,dts,xguess,dtguess,xmin,dtmin,xmax,dtmax,sigma0s,alphas,kmax,A,D,E).
  set xguess to KILout["xguess"].
  set dtguess to KILout["dtguess"].
  set A to KILout["A"].
  set D to KILout["D"].
  set E to KILout["E"].

  //PLATE 5-5
  set rs to 1 + 2*(b0*A + sigma0s*D*E).
  set b4 to 1/rs.

  if n > 0 {
    set xc to f6*(xguess+n*xP).
    set dtc to f4*(dtguess+n*Ps).
  }
  else {
    set xc to f6*xguess.
    set dtc to f4*dtguess.
  }

  set last["dtcp"] to dtc.
  set last["xcp"] to xc.
  set last["A"] to A.
  set last["D"] to D.
  set last["E"] to E.

  set F to 1 - 2*A.
  set Gs to 2*(D*E + sigma0s*A).
  set Fts to -2*b4*D*E.
  set Gt to 1 - 2*b4*A.

  set r to r0m*(F*ir0 + Gs*v0s).
  set v to f2*(Fts*ir0 + Gt*v0s).

  set out to lexicon().
  set out["r"] to r.
  set out["v"] to v.
  set out["last"] to last.

  return out.
}

function KTTI {
  parameter xarg,s0s,a,kmax.
  set u1 to USS(xarg,a,kmax).

  set zs to 2*u1.
  set E to 1 - 0.5*a*zs^2.
  set W to sqrt(max(0.5+E/2,0)).
  set D to W*zs.
  set A to D^2.
  set B to 2*(E+s0s*D).
  set Q to QCF(W).
  set t to D*(B+A*Q).

  set out to lexicon().
  set out["t"] to t.
  set out["A"] to A.
  set out["D"] to D.
  set out["E"] to E.

  return out.
}

function USS {
  parameter xarg,a,kmax.

  set du1 to xarg/4.
  set u1 to du1.
  set f7 to -a*du1^2.
  set k to 1.
  until k >= kmax {
    set du1 to (f7*du1) / (k*(k-1)).
    set u1old to u1.
    set u1 to u1 + du1.
    if u1 = u1old {
      break.
    }
    set k to k+1.
  }
  return u1.
}

function QCF {
  parameter W.

  if W<1 set xq to 21.04 - 13.04*W.
  else if W<4.625 set xq to (5/3)*(2*W+5).
  else if W<13.846 set xq to (10/7)*(W+12).
  else if W<44 set xq to 0.5*(W+60).
  else if W<100 set xq to 0.25*(W+164).
  else set xq to 70.

  set b to 0.
  set y to (W-1)/(W+1).
  set j to floor(xq).
  set b to y/(1+(j-1)/(j+2)*(1-b)).
  until j <= 2 {
    set j to j-1.
    set b to y/(1+(j-1)/(j+2)*(1-b)).
  }

  set out to 1/W^2 * (1 + (2-b/2) / (3*W*(W+1))).
  return out.
}

function KIL {
  parameter imax,dts,xguess,dtguess,xmin,dtmin,xmax,dtmax,s0s,a,kmax,A,D,E.

  set i to 1.
  until i >= imax {
    set dterror to dts - dtguess.
    if abs(dterror) <  1e-6 {
        break.
    }

    set SIout to SI(dterror,xguess,dtguess,xmin,dtmin,xmax,dtmax).
    set dxs to SIout["dxs"].
    set xmin to SIout["xmin"].
    set dtmin to SIout["dtmin"].
    set xmax to SIout["xmax"].
    set dtmax to SIout["dtmax"].

    set xold to xguess.
    set xguess to xguess + dxs.

    if xguess = xold break.

    set dtold to dtguess.

    set KTTIout to KTTI(xguess,s0s,a,kmax).
    set dtguess to KTTIout["t"].
    set A to KTTIout["A"].
    set D to KTTIout["D"].
    set E to KTTIout["E"].

    if dtguess = dtold break.

    set i to i+1.
  }
  set out to lexicon().
  set out["xguess"] to xguess.
  set out["dtguess"] to dtguess.
  set out["A"] to A.
  set out["D"] to D.
  set out["E"] to E.

  return out.
}

function SI {
  parameter dterror,xguess,dtguess,xmin,dtmin,xmax,dtmax.

  set etp to 1e-6.

  set dtminp to dtguess-dtmin.
  set dtmaxp to dtguess-dtmax.

  if abs(dtminp) < etp or abs(dtmaxp) < etp {
    set dxs to 0.
  }
  else {
    if dterror < 0 {
      set dxs to (xguess-xmax)*(dterror/dtmaxp).
      if (xguess+dxs) <= xmin {
        set dxs to (xguess-xmin)*(dterror/dtminp).
      }
      set xmax to xguess.
      set dtmax to dtguess.
    }
    else {
      set dxs to (xguess-xmin)*(dterror/dtminp).
      if (xguess+dxs) >= xmax {
        set dxs to (xguess-xmax)*(dterror/dtmaxp).
      }
      set xmin to xguess.
      set dtmin to dtguess.
    }
  }
  set out to lexicon().
  set out["dxs"] to dxs.
  set out["xmin"] to xmin.
  set out["dtmin"] to dtmin.
  set out["xmax"] to xmax.
  set out["dtmax"] to dtmax.
  return out.
}
