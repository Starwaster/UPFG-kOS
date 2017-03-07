function falcon9 {
  local vehicle is list(0).
  local stage_ is lexicon().
  set stage_["maxT"] to 60*6 + 33. 
  return stage_.
}
