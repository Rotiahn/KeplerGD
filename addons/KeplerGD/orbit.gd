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
@export_group("OrbitalElements")
@export_range(0.0,1.79769e308,1.0e-14) var a: float = 0  							## Semi-major axis, a (m)
@export var a_units: LengthEnum = 1                 ##Units for distance (m)
@export_range(0.0,1.79769e10,1.0e-14) var ecc: float = 0 							## Eccentricity, e, (circular: ecc=0, eliptical: 0<ecc<1, parabolic: ecc=1, hyperbolic: e>1)
@export_range(0,360,1.0e-14) var mAnomalyDeg: float = 0		## Mean anomaly, M
@export_range(0,360,1.0e-14) var rotIDeg: float = 0 		## Inclination w.r.t xz-plane, i       
@export_range(0,360,1.0e-14) var rotWDeg: float = 0			## Argument of Perifocus, w
@export_range(0,360,1.0e-14) var rotOmegDeg: float = 0 		## Longitude of Ascending Node, OMEGA  
#@export var peri: int = 0 							## Periapsis distance, q (m)
#@export var mu : float=0 							## Standard Gravitational Parameter  (m^3/s^2) 
#@export var apo: int = 0  							## Apoapsis distance (m)
#@export var T: float = 0							## Sidereal orbit period (s)
#@export var meanMotion: float = 0 					## Mean motion, n (rad/s)
#@export var periT: float = 0 						## Time of periapsis (s)
#@export_group("")


# public variables
var mAnomaly: float = mAnomalyDeg * PI / 180		## Mean anomaly, M
var rotI: float = rotIDeg * PI / 180 		## Inclination w.r.t xz-plane, i       
var rotW: float = rotWDeg * PI / 180		## Argument of Perifocus, w
var rotOmeg: float = rotOmegDeg * PI / 180 	## Longitude of Ascending Node, OMEGA  
var mu : float=0 							## Standard Gravitational Parameter  (m^3/s^2) 
var peri: float = 0 						## Periapsis distance, q (m)
var apo: float = 0  						## Apoapsis distance (m)
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
func _ready() -> void:
    print('a')
    print(a)
    print(a_units)

    match a_units:
        LengthEnum.METER:
            pass
        LengthEnum.KILOMETER:
            a *= KILOMETER
            a_units=LengthEnum.METER
        LengthEnum.KM:
            a *= KM
            a_units=LengthEnum.METER
        LengthEnum.AU: 
            a*= AU
            a_units=LengthEnum.METER
        LengthEnum.LIGHT_SECOND:
            a*= LIGHT_SECOND
            a_units=LengthEnum.METER
        LengthEnum.LIGHT_MINUTE:
            a*= LIGHT_MINUTE
            a_units=LengthEnum.METER
        LengthEnum.LIGHT_HOUR: 
            a*= LIGHT_HOUR
            a_units=LengthEnum.METER
        LengthEnum.LIGHT_DAY: 
            a*= LIGHT_DAY
            a_units=LengthEnum.METER
        LengthEnum.LIGHT_YEAR:
            a*=LIGHT_YEAR
            a_units=LengthEnum.METER
        _:
            var error_text:String = 'Error: invalid "A" calculations for '
            error_text += get_parent().name
            error_text += ' did not have a valid length_units'
            print(error_text)
            a = -INF

    print(a)

# Called when the node enters the scene tree for the first time.
func _init(
#    param_primary:AstroBody = primary,
#    param_a:float = a,
#    param_ecc:float = ecc,
#    param_mAnomaly:float = mAnomaly,
#    param_rotI:float = rotI,
#    param_rotW:float = rotW,
#    param_rotOmeg:float = rotOmeg,
)->void:
#    print(name)
#    primary = param_primary
#    a = param_a
#    ecc = param_ecc
#    mAnomaly = param_mAnomaly
#    rotI = param_rotI
#    rotW = param_rotW
#    rotOmeg = param_rotOmeg

     pass

func get_class()->String:
    return "Orbit"

##Part II: calculation functions:
    ##Source: http://www.bogan.ca/orbits/kepler/orbteqtn.html

# Calculate Periapsis
# * @function calcualtePeri
# * @returns {number} peri - (m) the periapsis of the orbit based on its orbital elements
# * @private
func _calculatePeri() -> float:
    if (ecc < 1): # circular or eliptical
        return (1.0 - ecc) * a
    elif (ecc == 1): # parabolic
        return a
    elif (ecc > 1): # hyperbolic
        return (1.0 - ecc) * a
    else: #Something has gone wrong
        var error_text:String = 'Error: _calculatePeri() for '
        error_text += get_parent().name
        error_text += ' did not have a valid ecc'
        print(error_text)
        return -INF


# Calculate Apoapsis
# * @function calculateApo
# * @returns {number} apo - (m) the apoapsis of the orbit based on its orbital elements
# * Infinity for parabolic and hyperbolic orbits
# * @private
func _calculateApo() -> float:
    if (ecc < 1): # circular or eliptical
        return (1 + ecc) * a
    elif (ecc >= 1): # parabolic or hyperbolic
        return INF
    else: #Something has gone wrong
        var error_text:String = 'Error: _calculateApo() for '
        error_text += get_parent().name
        error_text += ' did not have a valid ecc'
        print(error_text)
        return -INF

# Calculate Period (Sidereal Year)
# * @function calculateT
# * @returns {number} T - (s) The period (sidereal Year) for the orbit.
# * @private
func _calculateT() -> float:
    if (ecc < 1): # circular or eliptical
        return 2 * PI * pow( ( (a*a*a)/(mu) ) , 0.5)
    elif (ecc >= 1): # parabolic or hyperbolic
        return INF
    else: #Something has gone wrong
        var error_text:String = 'Error: _calculateT() for '
        error_text += get_parent().name
        error_text += ' did not have a valid ecc'
        print(error_text)
        return -INF

# Calculate Mean Motion
# * @function calculateMeanMotion
# * @returns {number} meanMotion - (rad/s) the Mean motion of the orbiting body
# * @private
func _calculateMeanMotion() -> float:
    if (ecc < 1): # circular or eliptical
        return pow( ( (a*a*a)/(mu) ) , -0.5)
    elif (ecc >= 1): #parabolic or hyperbolic
        return INF
    else: #Something has gone wrong
        var error_text:String = 'Error: _calculateMeanMotion() for '
        error_text += get_parent().name
        error_text += ' did not have a valid ecc'
        print(error_text)
        return -INF

# Calculate True Anomaly
# * @function calculateTAnomaly
# * @returns {number} tAnomaly - (rad) True anomaly, nu
# * @private
func _calculateTAnomaly() -> float:
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
        #@warning_ignore("integer_division")
        var A:float = (3.0/2.0) * mAnomaly;
        var B:float = pow(  A + sqrt( (A*A) + 1 )  , (1.0 / 3.0) );
        tAnomaly = 2 * atan( B - (1/B) );

        return tAnomaly
    elif (ecc >= 1):  #hyperbolic
        # cosh(F) = (ecc + cos(tAnomaly)) / (1+ecc*cos(tAnomaly))          //http://www.bogan.ca/orbits/kepler/orbteqtn.html
        # Using analogous method as elliptical solution
        E = _calculateE()
        var tanh_tAnomaly:float =  sqrt( (1+ecc) )* sin( E/2 ) / sqrt( (1-ecc) )* cos( E/2 )
        tAnomaly = 2 * atanh(tanh_tAnomaly)
        return tAnomaly
    else: #Something has gone wrong
        var error_text:String = 'Error: _calculateTAnomaly() for '
        error_text += get_parent().name
        error_text += ' did not have a valid ecc'
        print(error_text)
        return -INF

# Calculate Time of Periapsis
# * @function calculatePeriT
# * @returns {number} periT - (s) Time of periapsis
# * @private
func _calculatePeriT() -> float:
    if (ecc < 1): #circular or eliptical
        return pow( ( (a*a*a)/(mu) ) , 0.5) * mAnomaly
    elif (ecc == 1): #parabolic
        return pow( ( (2*a*a*a)/(mu) ) , 0.5) * mAnomaly
    elif (ecc > 1): #hyperbolic
        return pow( ( ((-a)*(-a)*(-a))/(mu) ) , 0.5) * mAnomaly;
    else: #Something has gone wrong
        var error_text:String = 'Error: _calculatePeriT() for '
        error_text += get_parent().name
        error_text += ' did not have a valid ecc'
        print(error_text)
        return -INF

# Calculate Eccentric anomaly
# * @function calculateE
# * @returns {number} E - (rad)   Eccentric anomaly
# * @private
func _calculateE() -> float:
    if (ecc == 0): #circular
        return mAnomaly
    elif (ecc < 1): # eliptical, per guidance from Markus
        var M:float = mAnomaly
        E = M
        var i:int = 0
        while (abs(E - (ecc * sin(E)) - M) > 0.0000000001):
            E -= (E - ecc * sin(E) - M) / (1 - ecc * cos(E))
            i += 1
            if (i>=1000):
                push_error('took too long to determine E for '+(str(get_instance_id())))
        return E
    elif (ecc == 1): #parabolic
        #D = tan(tAnomaly/2);
        tAnomaly = _calculateTAnomaly()
        var D:float = tan(tAnomaly/2)
        return D
    elif (ecc > 1): # hyperbolic
        var M:float = mAnomaly;
        #M = ecc * sinh(F) - F
        #0 = ecc * sinh(F) - F - M
        var F:float = M
        var i:int = 0
        while (abs((ecc * sinh(F)) - F - M) > 0.0000000001):
            F-= (ecc * sinh(F) - F - M) / (ecc * cosh(F) - 1)
            i += 1
            if (i>=1000):
                push_error('took too long to determine E for '+(str(get_instance_id())))
        return F
    else: #Something has gone wrong
        var error_text:String = 'Error: _calculateE() for '
        error_text += get_parent().name
        error_text += ' did not have a valid ecc'
        print(error_text)
        return -INF

# update Orbital Elements based on Cartesian Position and Velocity
# * @function keplerize
# * @param {Vector3} position - Vector of current position
# * @param {Vector3} velocity - Vector of current velocity
# * @see {@link http://microsat.sm.bmstu.ru/e-library/Ballistics/kepler.pdf}
# * @private
func _keplerize(new_mu:float,position:Vector3,velocity:Vector3) -> Dictionary:
    #var r = Vector3(position) # m
    #var v = Vector3(velocity) # m/s

    #This Function is meant to be used for determining orbits from points and velocities.
    #It is primarily used for applying velocity changes and calculating resulting new orbital elements
    var ang_momentum:Vector3 = position.cross(velocity) # m^2/s

    var rotOmeg2:float = atan( ang_momentum.x / (-ang_momentum.z) ) #radians
    var rotI2:float = atan( sqrt(ang_momentum.x*ang_momentum.x + ang_momentum.z*ang_momentum.z) / ang_momentum.y ) #radians

    #rotate position into orbital frame:
    var r:Vector3 = Vector3(position); # m

    var axis_omeg:Vector3 = Vector3(0,1,0) # Y+ = North
    r.rotated(axis_omeg, rotOmeg2)

    var axis_i:Vector3 = Vector3(1,0,0)
    r.rotated(axis_i,rotI2)

    #determine argument of latitude u, where tan(u) = tan(w+v) = p2/p1
    var arg_lat:float = atan(r.z/r.x)  #radians

    # a = (GM * r) / (2GM - r*velocity_scalar^2)
    # e = SQRT(1 - (h^2 / (GM*a))

    var mu2:float = new_mu;  #  m^3/s^2
    var r_scalar:float = position.length()  # m
    var v_scalar:float = velocity.length()  # m
    var h_scalar:float = ang_momentum.length(); # m^2/s
    var a2:float = (mu2 * r_scalar) / (2*mu2 - r_scalar*(v_scalar*v_scalar));  #  m^4/s^2  / (m^3/s^2 - m*m/s*m/s) = m^4/s^2  / m^3/s^2 = m
    var ecc2:float = sqrt(1 - (h_scalar*h_scalar / (mu2*a2)));  #  m^2/s * m^2/s  / (m^3/s^2 * m) = m^4/s^2 / m^4/s^2 = no units

    # radial velocity = position DOT velocity / position_scalar
    var rad_v:float = position.dot(velocity) / r_scalar;  #  m*m/s / m = m/s

    # cos(E) = (a-r)/(ae)
    # sin(E) = (r_scalar*rad_v)/(ecc*sqrt(mu*a))
    var sin_E:float = (a2-r_scalar)/(a2*ecc2);
    var cos_E:float = (r_scalar*rad_v)/(ecc*sqrt(mu2*a2));

    # tan(v) = ( sqrt(1-ecc*ecc)*sin_E )/( cos_E - ecc )
    var v:float = atan( (sqrt(1-ecc2*ecc2)*sin_E )/( cos_E - ecc2 ) );  #radians

    # u = w+v = arg_lat
    # w = arg_lat-v
    var rotW2:float = arg_lat - v; # radians

    # E - ecc*sin(E) = M

    var E2:float = asin(sin_E)  #radians
    var mAnomaly2:float = E2 - ecc2*sin_E;  #radians

    var elements:Dictionary = {	 'a'       :a2
                        ,'ecc'     :ecc2
                        ,'mAnomaly':mAnomaly2
                        ,'rotI'    :rotI2
                        ,'rotW'    :rotW2
                        ,'rotOmeg' :rotOmeg2
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
func _reverseRotations(vector:Vector3)-> Vector3:
    #NOTE: XZ plane is the (typical) plane of reference with X+ axis = reference direction and Y+ axis = "north"

    #Part I: Rotate orbit around y world-axis by angle -rotW so that periapsis lines up with reference direction
    var axisW:Vector3 = Vector3(0,1,0)
    vector.rotated(axisW,-rotW)

    #Part II: Rotate orbital plane around x world-axis by angle -rotI so that orbital plane lines up with reference plane
    var axisI:Vector3 = Vector3(1,0,0)
    vector.rotated(axisI,-rotI)

    #Part III: Rotate orbital plane around y world-axis by angle -rotOmeg so that ascending node lines up with reference direction
    var axisOmeg:Vector3 = Vector3(0,1,0)
    vector.rotated(axisOmeg,-rotOmeg);

    return vector;


# Get Cartesian position (x,y,z)
# * @function getPosition
# * @returns {KEPLER.Vector3} - Returns a KEPLER.Vector3 which defines the position in the orbit (INCORPORATES PRIMARY)
# * @see {@link http://microsat.sm.bmstu.ru/e-library/Ballistics/kepler.pdf}
# * @public
func getPosition() -> Vector3:

    #Part I: Update Orbital Elements
    #this.updateAllElements();
    updateElement("E")

    #Part II: Create initial elipse
    var position:Vector3 = Vector3(
         a*cos(E) -a*ecc
        ,a*sqrt(1-(ecc*ecc))*sin(E)
        ,0
    )

    #Part III: Conduct rotations (reversed):
    var positionFinal:Vector3 = _reverseRotations(position)

    #Part IV: Add position vector of primary:
    
    #if primary is NILL, then it's position is Vector3(0,0,0,0)
    if (primary == null):
        #positionFinal += Vector3(0,0,0)
        pass
    else:
        positionFinal += primary.getPosition()
    

    return positionFinal


# Get Cartesian velocity (x,y,z)
# * @function getVelocity
# * @returns {KEPLER.Vector3} - Returns a KEPLER.Vector3 which defines the position in the orbit (INCORPORATES PRIMARY)
# * @see {@link http://microsat.sm.bmstu.ru/e-library/Ballistics/kepler.pdf}
# * @public
func getVelocity() -> Vector3:

    #Part I: Update Orbital Elements
    #this.updateAllElements();
    updateElement("E")
    updateElement("meanMotion")

    #Part II: Create initial elipse
    var velocity:Vector3 = Vector3(
         ( (meanMotion*a)/( 1-(ecc*cos(E)) ) )*( -sin(E) )
        ,( (meanMotion*a)/( 1-(ecc*cos(E)) ) )*( sqrt(1-(ecc*ecc))*cos(E) )
        ,0
    )

    #Part III: Conduct rotations (reversed):
    var velocityFinal:Vector3 = _reverseRotations(velocity);

    #Part IV: Add position vector of primary:
    velocityFinal +=primary.getVelocity()

    return velocityFinal;


# remaining built-in virtual methods

# public methods

# Update single element
# * @member {object} updateElement
# * @example
# * //updates E (eccentric anomaly)
# * this.updateElement['E']()
# * @public
func updateElement(element:String) -> float:
    match element:
        "peri"        : peri = _calculatePeri();				return peri
        "apo"         : apo = _calculateApo();					return apo
        "T"           : T = _calculateT();						return T
        "meanMotion"  : meanMotion = _calculateMeanMotion(); 	return meanMotion
        "tAnomaly"    : tAnomaly = _calculateTAnomaly();		return tAnomaly
        "periT"       : periT = _calculatePeriT();				return periT
        "E"           : E = _calculateE(); 						return E
        _             : print("Error: updateElement invalid element"); return -INF

# Update all derivable elements
# * @function updateAllElements
# * @public
func updateAllElements() -> void:
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
func getElements() -> Dictionary:
    updateAllElements()
    var retObject:Dictionary = {
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
func clone() -> Orbit:
    #Part I: gather Orbital Elements
    updateAllElements()
    #var elements = getElements();

    #Part II: Create clone
    var new_clone:Orbit = Orbit.new()
    new_clone.primary = primary 
    new_clone.a = a
    new_clone.ecc = ecc
    new_clone.mAnomaly = mAnomaly
    new_clone.rotI = rotI
    new_clone.rotW = rotW
    new_clone.rotOmeg = rotOmeg


    return new_clone

 # Add Time: revolve object forward in time
# * @function addTime
# * @param {number} time - the time (in seconds) to adjust the object's movement
# * @returns {KEPLER.Orbit} - Returns this KEPLER.Orbit in it's new state after the transition
# * @public
func addTime(deltaTime:float) -> Orbit:
    #Adding Time can be completely accomplished with updating the mean Anomaly (mAnomaly) with a new value.
    #deltaMAnomaly = (deltaTime*meanMotion)%(2PI)
    updateElement("meanMotion"); # (rad/s)
    var deltaMAnomaly:float = fmod ( ( deltaTime * meanMotion ) , (2 * PI) )  # ( (s) * (rad/s) ) % (rad)
    mAnomaly = fmod (
            ( 
              fmod ( (mAnomaly+deltaMAnomaly),(2*PI)) 
            + (2*PI) 
            )
            ,(2*PI)
        )  # (rad), forces to always be between 0 and 2PI
    updateAllElements();
    return self;


# Subtract Time: revolve object backwards in time
# * @function subTime
# * @param {number} time - the time (in seconds) to adjust the object's movement
# * @returns {KEPLER.Orbit} - Returns this KEPLER.Orbit in it's new state after the transition
# * @public
func subTime(deltaTime:float)->Orbit:
    return addTime(-deltaTime)


# Add Cartesian velocity to this orbit to cause a change in orbital functions
# * @function addVelocity
# * @param {KEPLER.Vector3} deltaV - Vector to be added to the current velocity and adjust orbital elements
# * @returns {KEPLER.Orbit} - Returns a KEPLER.Vector3 which defines the position in the orbit (INCORPORATES PRIMARY)
# * @see {@link http://microsat.sm.bmstu.ru/e-library/Ballistics/kepler.pdf}
# * @public
func addVelocity(deltaV:Vector3)->Dictionary:
    #Part I: Get Cartesian elements
    var position:Vector3 = getPosition();
    var velocity:Vector3 = getVelocity();

    #Part II: Add deltaV to velocity;
    velocity += deltaV

    #Part III: Update orbital elements
    var result:Dictionary = _keplerize(mu,position,velocity);

    a           = result.a;
    ecc         = result.ecc;
    mAnomaly    = result.mAnomaly;
    rotI        = result.rotI;
    rotW        = result.rotW;
    rotOmeg     = result.rotOmeg;

    updateAllElements();

    return result;


# Check if this Orbit has an OrbitMesh object, if so, return the first one found
# * @function getOrbitMesh
# * @returns {KEPLER.OrbitMesh} - Returns a KEPLER.OrbitMesh which is a child of this Orbit
# * @public
func getOrbitMesh()->OrbitMesh:
    var children:Array = get_children()
    #var childQty:int = get_child_count()
    #var childTarget:OrbitMesh = null
    
    for child in children:
        if child.get_class() == "OrbitMesh":
            return child
        else:
            continue
    
    #No OrbitMesh found
    return null
    pass

# private methods

# subclasses
