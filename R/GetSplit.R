# ------------------------------------------------------------------------------
# GetSplit
# 2020/05/17
# Ryan D. Crawford
# ------------------------------------------------------------------------------
# Calculates the difference between the input time and the current time and
# prints the difference between them. The current time is returned.
# ------------------------------------------------------------------------------

GetSplit = function( inTime )
{
  curTime = Sys.time()
  diff = round( difftime( curTime, inTime, units = "mins"), 2 )
  cat( "  -- Finished in ", diff, " minutes\n", sep = '' )
  return( curTime )
}

# ------------------------------------------------------------------------------