neat_symmetric_cbar_breaks <- function(max_abs, n.breaks.approx=9) {
  possible_increments <- c(1, 2, 5, 10, 25, 50, 100, 250, 500, 1000, 2500, 5000, 10000, 25000, 50000)
  
  # Identify the increment which is closest to being able to neatly fit the
  # required number of breaks
  n.breaks.oneside <- floor(n.breaks.approx / 2) - 1
  increment <- possible_increments[which.min(abs(floor(max_abs/possible_increments - 0.3) - n.breaks.oneside))]
  print(max_abs/increment)
  
  # Generate sequence using the increment. Start at last multiple of increment before max_abs,
  # but add adjustment of 0.3 to avoid this last multiple being too close to the max_abs value
  last_break <- floor(max_abs/increment - 0.3) * increment
  breaks <- c(-max_abs, seq(-last_break, last_break, increment), max_abs)
  
  #return(length(breaks))
  return(breaks)
}

neat_symmetric_cbar_breaks_log <- function(max_abs, n.breaks.approx=9) {
  all_breaks <- c(1, 10, 100, 1000, 10000, 100000)
  
  n.breaks.oneside <- floor(n.breaks.approx / 2) - 1
  last_break_index <- which(all_breaks > max_abs)[1]
                               
  last_break_index <- valid_break_indices[length(valid_break_indices)]
  log_diff_to_max <- (log(all_breaks[last_break_index]) - log(max_abs))/log(max_abs)
  if (abs(log_diff_to_max) < 1/n.breaks.oneside) {
    last_break_index <- last_break_index - 1
  }
  
  first_break_index <- max(1, last_break_index - n.breaks.oneside + 1)
  
  positive_breaks <- all_breaks[seq(first_break_index, last_break_index)]
  breaks <- c(-max_abs, rev(-positive_breaks), 0, positive_breaks, max_abs)
  
  #return(length(breaks))
  return(breaks)
}



