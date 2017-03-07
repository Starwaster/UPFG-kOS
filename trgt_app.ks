function launchTargeting {
  parameter lat.
  parameter lon.
  parameter per.
  parameter apo.
  parameter inc.
  parameter lan.
  parameter slip.
  local azimuth is 90.
  if inc < lat {
    print "inclination below launch site latitude".
    set azimuth to 90.
  }
  else {
    local Binertial is arcsin(cos(inc)/cos(lat)).
    local vorbit is sqrt(body:mu/(bigR+per*1000)).
    local vEarthrot is (2*constant:pi*bigR/body:rotationperiod)*cos(lat).
    local vrotx is vorbit * sin(Binertial) - vEarthrot.
    local vroty is vorbit * cos(Binertial).
    set azimuth to arctan2(vrotx,vroty).
  }

  local rel_lng is arcsin(tan(lat)/tan(inc)).
  local geo_lng is lan + rel_lng - body:rotationangle.
  set geo_lng to mod(360+geo_lng,360).
  local node_angle is geo_lng - lon.
  set node_angle to mod(node_angle+360+slip,360).
  local launchtime is 0.
  set launchtime to (node_angle/360)*body:rotationperiod.


  local target_plane_normal is v(-sin(lan)*sin(inc),cos(lan)*sin(inc),-cos(inc)).
  //local target_plane_normal is v(0,0,-1).
  local target_velocity is 2/(bigR + per*1000).
  set target_velocity to target_velocity - 1/(bigR + (per+apo)*500).
  set target_velocity to sqrt(body:mu*target_velocity).
  local a_ is (bigR + per*1000 + bigR + apo*1000)/2.
  //local e_ is 0.
  //local i_ is inc.
  local orbit_parameter is list(a_,0,inc,lan,90,1).
  local trgt is lexicon().
  set trgt["parameter"] to orbit_parameter.
  set trgt["azimuth"] to azimuth.
  set trgt["radius"] to bigR+per*1000.
  set trgt["normal"] to target_plane_normal.
  set trgt["angle"] to 0.
  set trgt["velocity"] to target_velocity.
  set trgt["eta"] to launchtime.
  return trgt.
}
