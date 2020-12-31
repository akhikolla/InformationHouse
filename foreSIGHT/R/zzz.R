.onAttach <- function(lib, pkg)
{
  # unlockBinding(".ga.default", asNamespace("exposure")) 
  version <- read.dcf(file.path(lib, pkg, "DESCRIPTION"), "Version")
  
  if(interactive())
    {
      packageStartupMessage(
" ___   ___   ____      
 \\  \\ /  /  |    \\   foreSIGHT - The Exposure 
  \\     /   |  |> |  Space
   |    |   |  _ /   Generator
  /  /\\  \\  | |
 /__/  \\ _\\ |_| version ", version)
}
else
  { packageStartupMessage("Package 'foreSIGHT' version ", version) } 

  packageStartupMessage("Type 'citation(\"foreSIGHT\")' for citing this R package in publications.")
  invisible()
}
  