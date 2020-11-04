resistance_neutral <- function(ws_mean=ws_mean, pl_height = pl_height){

  
  # Some physical constants
  k = 0.41 #Von Karman constant
  z = pl_height+2 #wind measure height, we assume that the wind measure is 2 m over plant height
  b =  0.6
  a = 0.1
  #Aerodynamic conductance
  gaM = k^2 * ws_mean/
    (log((z-pl_height*b)/(pl_height*a))*log((z-pl_height*b)/(0.135*pl_height*a))) #aerodynamic conductance for momentum; Tan et al.2019 (eq. 14) + Chu et al. 2018 (Supporting information)
  
  #boundary layer conductance
  kB1 = 2 #Tan et al. 2019 eq. 20 
  u_friction = (k*ws_mean)/(log((z-pl_height*b)/(pl_height*a))+log(1.25)) #friction velocity
  gbN = (k*u_friction)/kB1 #boundary layer
  r.a = (1/gaM + 1/gbN) # total aerodynamic conductance [m s-1]

  return(r.a)
  
}