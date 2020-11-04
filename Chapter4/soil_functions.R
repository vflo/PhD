cosby <- function(percent.clay=18, percent.sand=43) {
  
  # Calculate soil hydraulic parameters based on Cosby et al (1984)
  # (Coefficients are derived from Cosby et al Table 4)
  
  # vwc is volumetric water content (m^3 of water per m^3 of soil)
  # suction is matric.potential (soil suction)
  
  # Density of water
  rho.water <- 1000. # kg per m^3
  
  # Acceleration of gravity
  gravity <- 9.81 # meters per second
  
  # Make sure texture percentages sum to 100
  percent.silt <- 100. - (percent.clay + percent.sand) 
  
  ##################
  
  # Coefficients from Cosby et al (1984) Table 4:
  bee.intercept <- 3.10
  d.bee.d.clay <- 0.1570
  d.bee.d.sand <- -0.003
  
  suction.sat.intercept <- 1.54
  dLog.suction.sat.d.sand <- -0.0095
  dLog.suction.sat.d.silt <- 0.0063
  
  k.sat.intercept <- -0.60
  dLog.k.sat.d.sand <- 0.0126
  dLog.k.sat.d.clay <- -0.0064
  
  porosity.intercept <- 50.5
  d.porosity.d.clay <- -0.037
  d.porosity.d.sand <- -0.142
  
  ##################
  
  # Clapp and Hornberger "b" exponent
  bee = bee.intercept + d.bee.d.clay * percent.clay + d.bee.d.sand * percent.sand
  
  # Saturated soil water suction psi_s (mm)
  suction.sat = 10 * (10^(suction.sat.intercept + 
                            dLog.suction.sat.d.sand * percent.sand +
                            dLog.suction.sat.d.silt * percent.silt))
  
  # Saturated soil hydraulic conductivity K_s (mm per hour)
  hydraulic.cond.sat = 25.4 * 10^(k.sat.intercept +
                                    dLog.k.sat.d.sand * percent.sand +
                                    dLog.k.sat.d.clay * percent.clay)
  
  # Volumetric soil water concentration at saturation point θs (porosity)
  vwc.sat = 0.01 * (porosity.intercept + 
                      d.porosity.d.clay * percent.clay +
                      d.porosity.d.sand * percent.sand)
  
  # Volumetric soil moisture concentration at Field Capacity (Bonan Fig 9.8)
  suction.field.capacity <- 1000.   # mm
  vwc.field.capacity <- vwc.sat * (suction.field.capacity / suction.sat) ^ (-1/bee)
  
  # Volumetric soil moisture concentration at wilting point (Bonan Fig 9.8)
  suction.wilt <- 150000   # mm
  vwc.wilt <- vwc.sat * (suction.wilt / suction.sat) ^ (-1/bee)
  
  return(list(bee=bee, suction.sat=suction.sat, hydraulic.cond.sat=hydraulic.cond.sat,
              vwc.sat=vwc.sat, vwc.field.capacity=vwc.field.capacity, vwc.wilt=vwc.wilt))
  
}



clapp.and.hornberger <- function(vwc=.4, percent.sand=43, percent.clay=18){ 
  
  # Use the method of Clapp and Hornberger (1978) to estimate the soil
  # matric potential ("suction") and hydraulic conductivity of a soil 
  # given the texture (%sand and %clay) and the volumetric water content (vwc)
  
  # First compute soil properties from texture using Cosby et al (1984)
  soil <- cosby(percent.sand=percent.sand, percent.clay=percent.clay)
  
  # soil suction as a function of vwc
  suction <- soil$suction.sat * (vwc / soil$vwc.sat)^-soil$bee
  
  # hydraulic conductivity as a function of vwc
  hydraulic.cond <- soil$hydraulic.cond.sat * 
    (vwc / soil$vwc.sat)^(2 * soil$bee + 3)
  
  return(list(vwc=vwc, percent.sand=percent.sand,
              percent.silt=(100-(percent.sand+percent.clay)),
              percent.clay=percent.clay,
              suction=suction, hydraulic.cond=hydraulic.cond,
              bee=soil$bee, suction.sat=soil$suction.sat, 
              hydraulic.cond.sat=soil$hydraulic.cond.sat,
              vwc.sat=soil$vwc.sat, vwc.field.capacity=soil$vwc.field.capacity, 
              vwc.wilt=soil$vwc.wilt))
  
}