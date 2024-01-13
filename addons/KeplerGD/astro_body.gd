@tool
class_name AstroBody
extends Kepler

# A function for creating AstroBody objects
# * @author Rotiahn / https://github.com/Rotiahn/
# * @class
# * @classdesc AstroBody is the root class for any object (such as Planets, Moons, Spacecraft)
# * @param {string} id - the ID of the AstroBody being created.
# * @param {number} mass - the mass (in kg) of the AstroBody being created.
# * @param {KEPLER.Orbit} orbit - the Orbit object for this AstroBody.
# * @example
# * //Gives an AstroBody for the Sun at the center of the solar System
# * EXAMPLE.Sol = new KEPLER.AstroBody(KEPLER.SOL_MASS,new KEPLER.NULL_ORBIT());
# * @example
# * //Gives an AstroBody for the Earth orbiting the sun
# * EXAMPLE.Earth = new KEPLER.AstroBody(
# *     5.97219e24                               //mass
# *     ,new KEPLER.Orbit(
# *         EXAMPLE.Sol                         //Primary
# *         ,1.000371833989169e+00*KEPLER.AU    //a
# *         ,1.704239716781501e-02              //ecc
# *         ,3.581891404220149e+02*KEPLER.DEGREE//mAnomaly
# *         ,2.669113820737183e-04*KEPLER.DEGREE//rotI
# *         ,2.977668064579176e+02*KEPLER.DEGREE//rotW
# *         ,1.639752443600624e+02*KEPLER.DEGREE//rotOmeg
# *     )
# * );
# * @module kepler

# signals
# enums
# constants
# @export variables

## PROPERTIES IN INSPECTOR ##

#@export var id: String = "AstroBody"		## ID to use for this AstroBody

#@export_exp_easing("positive_only") var mass2: float = 1.0e0
@export var mass: float = 1.0e0				     ##Mass of this AstroBody
@export var mass_units: MassEnum = MassEnum.KG   ##Units for Mass

@export var orbit: Orbit  					##The orbit object for this AstroBody

#@export_enum(
    #"KILOMETER",
    #"KM", 
    #"AU", 
    #"LIGHT_SECOND",
    #"LIGHT_MINUTE",
    #"LIGHT_HOUR", 
    #"LIGHT_DAY", 
    #"LIGHT_YEAR",
#) var units2: String
    #

#const LENGTHS_DICT: Dictionary = {
    #KILOMETER = 1000,
    #KM = LengthEnum.KILOMETER,
    #AU = 1.496e+11,
    #LIGHT_SECOND = 2.998e+8,
    #LIGHT_MINUTE = 1.799e+10,
    #LIGHT_HOUR = 1.079e+12,
    #LIGHT_DAY = 2.59e+13,
    #LIGHT_YEAR = 9.461e+15,    
#}


# public variables
var satellites:Array = []				##Array of satellites which have this AstroBody as primary
var primary:AstroBody

# private variables

# @onready variables

# optional built-in virtual _init method
func _init(
    #name: String = "AstroBody",
    #mass: float = mass,
    orbit: Orbit = orbit,
    
) -> void:
    print("AstroBody_init")
    print(name)
    print(mass)
    print(orbit)
    pass
    
# optional built-in virtual _enter_tree() method
# built-in virtual _ready method


func _ready()-> void:
    #If we don't have an orbit assigned, Give ourselves the NULL Orbit 
    #  This implies this AstroBody sits at center of system and doesn't move
    #print("AstroBody_Ready")
    print(name)
    #print(mass)
    #print(orbit)
    if orbit == null:
        orbit = Orbit.new()
    primary = orbit.primary

    #Convert Mass input from inspector to KG
    #print("mass")
    #print(mass)
    #print(mass_units)

    match mass_units:
        MassEnum.KG:
            #print("KG")
            pass
        MassEnum.TONNE:
            mass *= TONNE
            mass_units = MassEnum.KG
        MassEnum.EARTH_MASS:
            mass *= EARTH_MASS
            mass_units = MassEnum.KG
        MassEnum.SOL_MASS:
            mass *= SOL_MASS
            mass_units = MassEnum.KG
        _:
            var error_text:String = 'Error: invalid mass calculations for '
            error_text += get_parent().name
            error_text += ' did not have a valid mass_units'
            print(error_text)
            mass = -INF
 
    #print(mass)
    #Create various functions related to having an orbit
    #var primary = orbit.primary		##The Primary for this Astrobody is based on it orbit

# remaining built-in virtual methods
func get_class()->String:
    return "AstroBody"


# public methods

# @borrows KEPLER.Orbit.getElements as getElements 
func getElements()->Dictionary:
    return orbit.getElements()

# @borrows KEPLER.Orbit.getPosition as getPosition 
func getPosition() -> Vector3:
    return orbit.getPosition()

# @borrows KEPLER.Orbit.getVelocity as getVelocity 
func getVelocity() -> Vector3:
    return orbit.getVelocity()

# @borrows KEPLER.Orbit.addTime as addTime 
func addTime(seconds:float) -> Orbit:
    return orbit.addTime(seconds)

# @borrows KEPLER.Orbit.subTime as subTime 
func subTime(seconds:float) -> Orbit:
    return orbit.subTime(seconds)


# Create a new (clone) AstroBody with the same parameters at this one, and return it.
# * @function clone
# * @returns {KEPLER.AstroBody} - Returns a new KEPLER.AstroBody with the same parameters as this one.  The orbits will also be separate objects.
# * @example
# * //Returns AstroBodyB is a copy of AstroBodyA, but different objects
# * AstroBodyA = new KEPLER.AstroBody('1',50*KEPLER.TONNE,new KEPLER.Orbit({mass:KEPLER.SOL},100e3,0,0,0,0));
# * AstroBodyB = AstroBodyA.clone();
# * AstroBodyA === AstroBodyB; //false
# * //All True:
# * for (key in Object.keys(orbitA.getElements())) {console.log(key,':',(orbitA.getElements()[key]===orbitB.getElements()[key]));};
# * @public
func clone() -> AstroBody:
    #Part I: gather Orbital Elements
    #this.updateAllElements();
    var elements:Dictionary = getElements();
    #var id = id;
    #var mass = this.mass;
    #var satellites = this.satellites;

    var a2:float = elements["a"]
    var ecc2:float = elements["ecc"]
    var mAnomaly2:float = elements["mAnomaly"]
    var rotI2:float = elements["rotI"]
    var rotW2:float = elements["rotW"]
    var rotOmeg2:float = elements["rotOmeg"]

    #Part II: Create clone of Orbit
    var cloneOrbit:Orbit = orbit.clone()

    #Part III: Clone Astrobody
    #var cloneAstroBody:AstroBody = AstroBody.new(name,mass,cloneOrbit);
    var cloneAstroBody:AstroBody = AstroBody.new();
    cloneAstroBody.orbit=cloneOrbit
    for satellite:AstroBody in satellites:
        cloneAstroBody.addSatellite(satellite.clone());
    return cloneAstroBody;


# Add satellite
# * @function addSatellite
# * @param {KEPLER.AstroBody} satellite
# * @public
func addSatellite(satellite:AstroBody)->void:
    satellites.push_back(satellite)


# Remove satellite
# * Only satellites which are pointers to exactly the same object will be removed.
# * Satellites which are different, but identical values will not be removed
# * @function removeSatellite
# * @param {KEPLER.AstroBody} satellite
# * @public
func removeSatellite(satellite:AstroBody)->void:
    satellites = satellites.filter(func(x:AstroBody)->bool:return x != satellite)
   
# private methods
# subclasses
