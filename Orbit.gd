@tool
class_name Orbit
extends Kepler

# A Class to represent an Orbit around a particular AstroBody
# @author Rotiahn / https://github.com/Rotiahn/
# @class
# classdesc Orbit is a class for defining orbit parameters and characteristics useful to all AstroBodys
# @param {AstroBody}   primary     - The AstroBody around which this orbit exists
# @param {number}      a           - The semi-major axis of the orbit in meters
# @param {number}      ecc         - The eccentricity of the orbit
# @param {number}      mAnomaly    - The Mean Anomaly of the oribt
# @param {number}      rotI        - The inclination of the orbit (rotation of plane from horizontal)
# @param {number}      rotW        - The Argument of perifocus (rotation of orbit around normal of its inclined plane)
# @param {number}      rotOmeg     - The Longitude of Ascending Node (rotation of orbital plane around vertical axis)
# @example
# //returns circular orbit of earth
# var earthOrbit = new KEPLER.Orbit({mass:KEPLER.SOL_MASS},KEPLER.AU,0,0,0,0);
# @module kepler

# signals

# enums

# constants

# @export variables
## PROPERTIES IN INSPECTOR ##

##TODO: require an Astrobody object here instead of Node3D
@export var primary: Node3D #Astrobody

#Orbital Elements
@export_group("Orbital Elements")
@export var a: int = 0  							## Semi-major axis, a (m)
@export var ecc: float = 0 							## Eccentricity, e, (circular: ecc=0, eliptical: 0<ecc<1, parabolic: ecc=1, hyperbolic: e>1)
@export_range(0,360) var mAnomalyDeg: float = 0		## Mean anomaly, M
@export_range(0,360) var rotIDeg: float = 0 		## Inclination w.r.t xz-plane, i       
@export_range(0,360) var rotWDeg: float = 0			## Argument of Perifocus, w
@export_range(0,360) var rotOmegDeg: float = 0 		## Longitude of Ascending Node, OMEGA  
#@export var mu : float=0 							## Standard Gravitational Parameter  (m^3/s^2) 
#@export var peri: int = 0 							## Periapsis distance, q (m)
#@export var apo: int = 0  							## Apoapsis distance (m)
#@export var T: float = 0							## Sidereal orbit period (s)
#@export var meanMotion: float = 0 					## Mean motion, n (rad/s)
#@export var periT: float = 0 						## Time of periapsis (s)
@export_group("")


# public variables
var mAnomaly: float = mAnomalyDeg * PI / 180		## Mean anomaly, M
var rotI: float = rotIDeg * PI / 180 		## Inclination w.r.t xz-plane, i       
var rotW: float = rotWDeg * PI / 180			## Argument of Perifocus, w
var rotOmeg: float = rotOmegDeg * PI / 180 		## Longitude of Ascending Node, OMEGA  
var mu : float=0 							## Standard Gravitational Parameter  (m^3/s^2) 
var peri: int = 0 							## Periapsis distance, q (m)
var apo: int = 0  							## Apoapsis distance (m)
var T: float = 0							## Sidereal orbit period (s)
var meanMotion: float = 0 					## Mean motion, n (rad/s)
var tAnomaly: float = 0						## True anomaly, nu (rad)
var periT: float = 0 						## Time of periapsis (s)
var E: float = 0							## Eccentric anomaly (rad)

# private variables

# @onready variables

# optional built-in virtual _init method

# optional built-in virtual _enter_tree() method

# built-in virtual _ready method
# Called when the node enters the scene tree for the first time.
func _ready():
	pass # Replace with function body.

##Part II: calculation functions:
	##Source: http://www.bogan.ca/orbits/kepler/orbteqtn.html

# Calculate Periapsis
# * @function calcualtePeri
# * @returns {number} peri - (m) the periapsis of the orbit based on its orbital elements
# * @private
func _calculatePeri():
	if (ecc < 1): # circular or eliptical
		return (1 - ecc) * a
	elif (ecc == 1): # parabolic
		return a
	elif (ecc > 1): # hyperbolic
		return (1 - ecc) * a


# Calculate Apoapsis
# * @function calculateApo
# * @returns {number} apo - (m) the apoapsis of the orbit based on its orbital elements
# * Infinity for parabolic and hyperbolic orbits
# * @private
func _calculateApo():
	if (ecc < 1): # circular or eliptical
		return (1 + ecc) * a
	elif (ecc >= 1): # parabolic or hyperbolic
		return INF


# Calculate Period (Sidereal Year)
# * @function calculateT
# * @returns {number} T - (s) The period (sidereal Year) for the orbit.
# * @private
func _calculateT():
	if (ecc < 1): # circular or eliptical
		return 2 * PI * pow( ( (a*a*a)/(mu) ) , 0.5)
	elif (ecc >= 1): # parabolic or hyperbolic
		return INF


# Calculate Mean Motion
# * @function calculateMeanMotion
# * @returns {number} meanMotion - (rad/s) the Mean motion of the orbiting body
# * @private
func _calculateMeanMotion():
	if (ecc < 1): # circular or eliptical
		return pow( ( (a*a*a)/(mu) ) , -0.5)
	elif (ecc >= 1): #parabolic or hyperbolic
		return INF


# Calculate True Anomaly
# * @function calculateTAnomaly
# * @returns {number} tAnomaly - (rad) True anomaly, nu
# * @private
func _calculateTAnomaly():
	if (ecc == 0): # circular
		return mAnomaly
	elif (ecc < 1): # eliptical
		E = _calculateE()
		tAnomaly = 2 * atan2( 
				 sqrt( (1+ecc) ) * sin( E/2 ) 
				,sqrt( (1-ecc) ) * cos( E/2 ) 
		) # //https://en.wikipedia.org/wiki/True_anomaly
		return tAnomaly
	elif (ecc == 1): #parabolic
		#var periT = Math.pow( ( (2*a*a*a)/(mu) ) , 0.5) * mAnomaly        //http://www.bogan.ca/orbits/kepler/orbteqtn.html
		#var peri  = a;                                                    //see Note 1 above
		#var A = (3/2) * Math.sqrt(mu / (2 * peri*peri*peri ) ) * (periT)  //https://en.wikipedia.org/wiki/Parabolic_trajectory#Barker.27s_equation
		#var B = Math.cbrt( A + Math.sqrt( (A*A) + 1 ) )                   //https://en.wikipedia.org/wiki/Parabolic_trajectory#Barker.27s_equation
		#var tAnomaly = 2 * arctan (B - 1/B)                               //https://en.wikipedia.org/wiki/Parabolic_trajectory#Barker.27s_equation

		#var mu2p3 = mu/(2*peri*peri*per)

		#var periT = Math.sqrt( 1/mu2p3 ) * mAnomaly
		#var periT = mAnomaly / Math.sqrt( mu2p3 )

		#var A = (3/2) * Math.sqrt( mu2p3 ) * (mAnomaly / Math.sqrt( mu2p3 ) )
		#var A = (3/2) * mAnomaly
		@warning_ignore("integer_division")
		var A = (3.0/2.0) * mAnomaly;
		var B = pow(  A + sqrt( (A*A) + 1 )  , (1.0 / 3.0) );
		tAnomaly = 2 * atan( B - (1/B) );

		return tAnomaly
	elif (ecc >= 1):  #hyperbolic
		# cosh(F) = (ecc + cos(tAnomaly)) / (1+ecc*cos(tAnomaly))          //http://www.bogan.ca/orbits/kepler/orbteqtn.html
		# Using analogous method as elliptical solution
		E = _calculateE()
		var tanh_tAnomaly =  sqrt( (1+ecc) )* sin( E/2 ) / sqrt( (1-ecc) )* cos( E/2 )
		tAnomaly = 2 * atanh(tanh_tAnomaly)
		return tAnomaly


# Calculate Time of Periapsis
# * @function calculatePeriT
# * @returns {number} periT - (s) Time of periapsis
# * @private
func _calculatePeriT():
	if (ecc < 1): #circular or eliptical
		return pow( ( (a*a*a)/(mu) ) , 0.5) * mAnomaly
	elif (ecc == 1): #parabolic
		return pow( ( (2*a*a*a)/(mu) ) , 0.5) * mAnomaly
	elif (ecc > 1): #hyperbolic
		return pow( ( ((-a)*(-a)*(-a))/(mu) ) , 0.5) * mAnomaly;


# Calculate Eccentric anomaly
# * @function calculateE
# * @returns {number} E - (rad)   Eccentric anomaly
# * @private
func _calculateE():
	if (ecc == 0): #circular
		return mAnomaly
	elif (ecc < 1): # eliptical, per guidance from Markus
		var M = mAnomaly
		E = M
		var i = 0
		while (abs(E - (ecc * sin(E)) - M) > 0.0000000001):
			E -= (E - ecc * sin(E) - M) / (1 - ecc * cos(E))
			i += 1
			if (i>=1000):
				push_error('took too long to determine E for '+(str(get_instance_id())))
		return E
	elif (ecc == 1): #parabolic
		#D = tan(tAnomaly/2);
		tAnomaly = _calculateTAnomaly()
		var D = tan(tAnomaly/2)
		return D
	elif (ecc > 1): # hyperbolic
		var M = mAnomaly;
		#M = ecc * sinh(F) - F
		#0 = ecc * sinh(F) - F - M
		var F = M
		var i = 0
		while (abs((ecc * sinh(F)) - F - M) > 0.0000000001):
			F-= (ecc * sinh(F) - F - M) / (ecc * cosh(F) - 1)
			i += 1
			if (i>=1000):
				push_error('took too long to determine E for '+(str(get_instance_id())))
		return F

# remaining built-in virtual methods

# public methods

# private methods

# subclasses
