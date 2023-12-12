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

@export var primary: AstroBody

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
func _init(
    param_primary:AstroBody,
    param_a:int,
    param_ecc:float,
    param_mAnomaly:float,
    param_rotI:float,
    param_rotW:float,
    param_rotOmeg:float,
):
    primary = param_primary
    a = param_a
    ecc = param_ecc
    mAnomaly = param_mAnomaly
    rotI = param_rotI
    rotW = param_rotW
    rotOmeg = param_rotOmeg

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


# update Orbital Elements based on Cartesian Position and Velocity
# * @function keplerize
# * @param {Vector3} position - Vector of current position
# * @param {Vector3} velocity - Vector of current velocity
# * @see {@link http://microsat.sm.bmstu.ru/e-library/Ballistics/kepler.pdf}
# * @private
func _keplerize(new_mu,position,velocity):
    #var r = Vector3(position) # m
    #var v = Vector3(velocity) # m/s

    #This Function is meant to be used for determining orbits from points and velocities.
    #It is primarily used for applying velocity changes and calculating resulting new orbital elements
    var ang_momentum = position.cross(velocity) # m^2/s

    var rot_omeg = atan( ang_momentum.x / (-ang_momentum.z) ) #radians
    var rot_i = atan( sqrt(ang_momentum.x*ang_momentum.x + ang_momentum.z*ang_momentum.z) / ang_momentum.y ) #radians

    #rotate position into orbital frame:
    var r = Vector3(position); # m

    var axis_omeg = Vector3(0,1,0) # Y+ = North
    r.rotated(axis_omeg, rot_omeg)

    var axis_i = Vector3(1,0,0)
    r.rotated(axis_i,rot_i)

    #determine argument of latitude u, where tan(u) = tan(w+v) = p2/p1
    var arg_lat = atan(r.z/r.x)  #radians

    # a = (GM * r) / (2GM - r*velocity_scalar^2)
    # e = SQRT(1 - (h^2 / (GM*a))

    var mu2 = new_mu;  #  m^3/s^2
    var r_scalar = position.length()  # m
    var v_scalar = velocity.length()  # m
    var h_scalar = ang_momentum.length(); # m^2/s
    var a = (mu2 * r_scalar) / (2*mu2 - r_scalar*(v_scalar*v_scalar));  #  m^4/s^2  / (m^3/s^2 - m*m/s*m/s) = m^4/s^2  / m^3/s^2 = m
    var ecc = sqrt(1 - (h_scalar*h_scalar / (mu2*a)));  #  m^2/s * m^2/s  / (m^3/s^2 * m) = m^4/s^2 / m^4/s^2 = no units

    # radial velocity = position DOT velocity / position_scalar
    var rad_v = position.dot(velocity) / r_scalar;  #  m*m/s / m = m/s

    # cos(E) = (a-r)/(ae)
    # sin(E) = (r_scalar*rad_v)/(ecc*sqrt(mu*a))
    var sin_E = (a-r_scalar)/(a*ecc);
    var cos_E = (r_scalar*rad_v)/(ecc*sqrt(mu2*a));

    # tan(v) = ( sqrt(1-ecc*ecc)*sin_E )/( cos_E - ecc )
    var v = atan( (sqrt(1-ecc*ecc)*sin_E )/( cos_E - ecc ) );  #radians

    # u = w+v = arg_lat
    # w = arg_lat-v
    var rot_w = arg_lat - v; # radians

    # E - ecc*sin(E) = M

    var E = asin(sin_E)  #radians
    var M = E - ecc*sin_E;  #radians

    var elements = {	 'a'       :a
                        ,'ecc'     :ecc
                        ,'mAnomaly':mAnomaly
                        ,'rotI'    :rotI
                        ,'rotW'    :rotW
                        ,'rotOmeg' :rotOmeg
                        }
    return elements


# Apply Reverse Rotations for Kepler elements -> Cartesian Elements (x,y,z)
# * Used by this.getPosition() and this.getVelocity()
# * NOTE: XZ plane is the plane of reference with X+ axis = reference direction and Y+ axis = "north"
# * @function reverseRotations
# * @param {KEPLER.Vector3} vector - the vector (relative to the orbital plane) to be rotated to match the world reference
# * @returns {KEPLER.Vector3} - Returns a KEPLER.Vector3 which defines the position in the orbit in world reference frame (RELATIVE TO PRIMARY)
# * @see {@link http://microsat.sm.bmstu.ru/e-library/Ballistics/kepler.pdf}
# * @private
func _reverseRotations(vector):
    #NOTE: XZ plane is the (typical) plane of reference with X+ axis = reference direction and Y+ axis = "north"

    #Part I: Rotate orbit around y world-axis by angle -rotW so that periapsis lines up with reference direction
    var axisW = Vector3(0,1,0)
    vector.rotated(axisW,-rotW)

    #Part II: Rotate orbital plane around x world-axis by angle -rotI so that orbital plane lines up with reference plane
    var axisI = Vector3(1,0,0)
    vector.rotated(axisI,-rotI)

    #Part III: Rotate orbital plane around y world-axis by angle -rotOmeg so that ascending node lines up with reference direction
    var axisOmeg = Vector3(0,1,0)
    vector.rotated(axisOmeg,-rotOmeg);

    return vector;


# Get Cartesian position (x,y,z)
# * @function getPosition
# * @returns {KEPLER.Vector3} - Returns a KEPLER.Vector3 which defines the position in the orbit (INCORPORATES PRIMARY)
# * @see {@link http://microsat.sm.bmstu.ru/e-library/Ballistics/kepler.pdf}
# * @public
func getPosition():

    #Part I: Update Orbital Elements
    #this.updateAllElements();
    updateElement("E")

    #Part II: Create initial elipse
    var position = Vector3(
         a*cos(E) -a*ecc
        ,a*sqrt(1-(ecc*ecc))*sin(E)
        ,0
    )

    #Part III: Conduct rotations (reversed):
    var positionFinal = _reverseRotations(position)

    #Part IV: Add position vector of primary:
    positionFinal.add(primary.getPosition())

    return positionFinal


# Get Cartesian velocity (x,y,z)
# * @function getVelocity
# * @returns {KEPLER.Vector3} - Returns a KEPLER.Vector3 which defines the position in the orbit (INCORPORATES PRIMARY)
# * @see {@link http://microsat.sm.bmstu.ru/e-library/Ballistics/kepler.pdf}
# * @public
func getVelocity():

    #Part I: Update Orbital Elements
    #this.updateAllElements();
    updateElement("E")
    updateElement("meanMotion")

    #Part II: Create initial elipse
    var velocity = Vector3(
         ( (meanMotion*a)/( 1-(ecc*cos(E)) ) )*( -sin(E) )
        ,( (meanMotion*a)/( 1-(ecc*cos(E)) ) )*( sqrt(1-(ecc*ecc))*cos(E) )
        ,0
    )

    #Part III: Conduct rotations (reversed):
    var velocityFinal = _reverseRotations(velocity);

    #Part IV: Add position vector of primary:
    velocityFinal.add(primary.getVelocity())

    return velocityFinal;


# remaining built-in virtual methods

# public methods

# Update single element
# * @member {object} updateElement
# * @example
# * //updates E (eccentric anomaly)
# * this.updateElement['E']()
# * @public
func updateElement(element):
    match element:
        "peri"        : peri = _calculatePeri();				return peri
        "apo"         : apo = _calculateApo();					return apo
        "T"           : T = _calculateT();						return T
        "meanMotion"  : meanMotion = _calculateMeanMotion(); 	return meanMotion
        "tAnomaly"    : tAnomaly = _calculateTAnomaly();		return tAnomaly
        "periT"       : periT = _calculatePeriT();				return periT
        "E"           : E = _calculateE(); 						return E


# Update all derivable elements
# * @function updateAllElements
# * @public
func updateAllElements():
    updateElement("peri")
    updateElement("apo")
    updateElement("T")
    updateElement("meanMotion")
    updateElement("tAnomaly")
    updateElement("periT")
    updateElement("E")
    
    
# Get all orbital Elements
# * @function getElements
# * @returns {Object} - Returns an object which includes all orbital elements.
# * @public
func getElements():
    updateAllElements()
    var retObject = {
            a          :a
        ,ecc        :ecc
        ,mAnomaly   :mAnomaly
        ,rotI       :rotI
        ,rotW       :rotW
        ,rotOmeg    :rotOmeg
        ,mu         :mu
        ,peri       :peri
        ,apo        :apo
        ,T          :T
        ,meanMotion :meanMotion
        ,tAnomaly   :tAnomaly
        ,periT      :periT
        ,E          :E
    }
    return retObject


# Create a new (clone) Orbit with the same parameters at this one, and return it.
# * @function clone
# * @returns {KEPLER.Orbit} - Returns a new KEPLER.Orbit with the same parameters as this one.
# * @example
# * //Returns orbitB is a copy of orbitB, but different objects
# * orbitA = new KEPLER.Orbit({mass:KEPLER.SOL},100e3,0,0,0,0);
# * orbitB = orbitA.clone();
# * orbitA === orbitB; //false
# * //All True:
# * for (key in Object.keys(orbitA.getElements())) {console.log(key,':',(orbitA.getElements()[key]===orbitB.getElements()[key]));};
# * @public
func clone():
    #Part I: gather Orbital Elements
    updateAllElements()
    var elements = getElements();

    #Part II: Create clone
    var clone = Orbit.new(
            primary
        ,a
        ,ecc
        ,mAnomaly
        ,rotI
        ,rotW
        ,rotOmeg
    )
    return clone

 # Add Time: revolve object forward in time
# * @function addTime
# * @param {number} time - the time (in seconds) to adjust the object's movement
# * @returns {KEPLER.Orbit} - Returns this KEPLER.Orbit in it's new state after the transition
# * @public
func addTime(deltaTime):
    #Adding Time can be completely accomplished with updating the mean Anomaly (mAnomaly) with a new value.
    #deltaMAnomaly = (deltaTime*meanMotion)%(2PI)
    updateElement("meanMotion"); # (rad/s)
    var deltaMAnomaly = ( deltaTime * meanMotion ) % (2 * PI);  # ( (s) * (rad/s) ) % (rad)
    mAnomaly = ( (mAnomaly+deltaMAnomaly)%(2*PI) + (2*PI) )%(2*PI);  # (rad), forces to always be between 0 and 2PI
    updateAllElements();
    return self;


# Subtract Time: revolve object backwards in time
# * @function subTime
# * @param {number} time - the time (in seconds) to adjust the object's movement
# * @returns {KEPLER.Orbit} - Returns this KEPLER.Orbit in it's new state after the transition
# * @public
func subTime(deltaTime):
    return addTime(-deltaTime)


# Add Cartesian velocity to this orbit to cause a change in orbital functions
# * @function addVelocity
# * @param {KEPLER.Vector3} deltaV - Vector to be added to the current velocity and adjust orbital elements
# * @returns {KEPLER.Orbit} - Returns a KEPLER.Vector3 which defines the position in the orbit (INCORPORATES PRIMARY)
# * @see {@link http://microsat.sm.bmstu.ru/e-library/Ballistics/kepler.pdf}
# * @public
func addVelocity(deltaV):
    #Part I: Get Cartesian elements
    var position = getPosition();
    var velocity = getVelocity();

    #Part II: Add deltaV to velocity;
    velocity.add(deltaV);

    #Part III: Update orbital elements
    var result = _keplerize(mu,position,velocity);

    a           = result.a;
    ecc         = result.ecc;
    mAnomaly    = result.mAnomaly;
    rotI        = result.rotI;
    rotW        = result.rotW;
    rotOmeg     = result.rotOmeg;

    updateAllElements();

    return result;


# private methods

# subclasses
