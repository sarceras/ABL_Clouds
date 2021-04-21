# Define constants.
rhow = 1000. # kg/m^3 water density
lambda_ = 2.45*10**6 # J/kg
rho = 1.225 # kg/m^3
g = 9.8
cp = 1005. # J/kg/K
P0 = 101.325*10**3 # Pa
Rd = 287.
Rv = 461.
z_wind = 2. # m
u_wind = 2. # m/s
kvc = 0.41 # von Karman constant
C0 = 400 # ppm C02 mean concentration
R_s0 = 1367.0 # Solar constant W/m^2
R_earth = 6371. # km
beta = 0.2
t_noon = 12
y_atm = 9
r = R_earth/y_atm

# Soil parameters.
mpsi_s = {"Sandy loam": -.7*10**-3, "Loamy sand": -.17*10.**-3., "Loam": -1.43*10.**-3., "Clay": -1.82*10**-3.} # (MPa)
b = {"Sandy loam": 4.9, "Loamy sand": 4.38, "Loam": 5.39, "Clay": 11.4} # (unitless)
Ks = {"Sandy loam": 80., "Loamy sand": 100., "Loam": 20., "Clay": 1.} # Saturated hydraulic condutivity (cm/day)

# Dictionaries per species.
alfas = {"Grass": 0.20, "Forest": 0.15, "White": 1.0, "Black": 0.0}
ds = {"Grass": 0.4, "Forest": 0.1, "Black": 1.}
height = {"Grass": 1., "Forest": 10.}
Zr = {"Grass": .6, "Forest": .65} # Rooting depth (m)
LAI = {"Grass": 1.71, "Forest": 5.} # Leaf area index (m2/m2)
ga = {"Grass": 0.01, "Forest": 0.02} # Atmospheric Conductance, per unit ground area (m/s)
gpmax = {"Grass": 11.7, "Forest": 0.056} # Plant conductance, per unit leaf area(um/(MPa s))
RAIW = {"Grass": 5.6, "Forest": 10.} # Root area index (m2/m2)