# ------------------------------------------------------------------------------
# GetFilePaths
# 2019/05/21
# Ryan D. Crawford
# ------------------------------------------------------------------------------
# Get the files in the input directory with the provided extension. If no 
# file extension is input, all files in the directory are returned. 
# ------------------------------------------------------------------------------

GetFilePaths = function( dir, ext )
{
  # Check that there is a '/' on the input directory
  if ( !grepl("/$", dir) && dir != '' ) dir = paste0( dir, '/' )
  
  # It no file extension was input, set the default extension to an 
  # empty string
  if ( missing(ext) ) ext = ''
  
  # Run the command to import the paths into R
  cmd   = paste0( "ls ", dir, '*', ext )
  files = system( cmd, intern = TRUE )
  
  if ( length(files) == 0 )
  {
    stop("No files withere found in directiory ", dir, "with extension '",
      ext, "'\n"   
      )
  }
  return( files )
}

# ------------------------------------------------------------------------------

