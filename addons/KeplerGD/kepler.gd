class_name Kepler
extends Node


#  A library for handling orbital mechanics math and calculations
# @author Rotiahn / https://github.com/Rotiahn/
# @namespace kepler
# 


# signals

# enums

# constants
const VERSION = '0.1.0'
#const MAX_SAFE_LENGTH = Number.MAX_SAFE_INTEGER; # 9007199254740991 m ~= 0.95 light years

#Constants - Distances
const UNIT_LENGTH = 'meter'
const KILOMETER = 1000
const KM = KILOMETER
const AU = 1.496e+11
const LIGHT_SECOND = 2.998e+8
const LIGHT_MINUTE = 1.799e+10
const LIGHT_HOUR = 1.079e+12
const LIGHT_DAY = 2.59e+13
const LIGHT_YEAR = 9.461e+15
enum LengthEnum { 
    METER=1,
    KILOMETER=2,
    KM=2, 
    AU=3, 
    LIGHT_SECOND=4,
    LIGHT_MINUTE=5,
    LIGHT_HOUR=6, 
    LIGHT_DAY=7, 
    LIGHT_YEAR=8,
}

#Constants - Time
const UNIT_TIME = 'second'
#const MAX_SAFE_TIME = Number.MAX_SAFE_INTEGER; # 9007199254740991 s ~= 285 Million years
const MINUTE = 60
const HOUR = 3600
const DAY = 86400
const YEAR = 3.154e+7
enum TimeEnum {
    SECOND=1,
    MINUTE=2,
    HOUR=3,
    DAY=4,
    YEAR=5
}
#Constants - Rotation
const UNIT_ROTATION = 'Radians'
const DEGREE = 0.0174533
#const PI = PI
const DEGREES_PER_DAY = DEGREE / DAY

#Constants - Mass
const UNIT_mass = 'Kilograms'
const TONNE = 1000
const TON = TONNE
const EARTH_MASS = 5.974e24
const SOL_MASS = 1.9891e30
enum MassEnum {
    KG=1,
    TONNE=2,
    EARTH_MASS=3,
    SOL_MASS=4,
}
#Constants - Physics
const G = 6.674e-11 # Nm^2 / kg^2 = kg*m*(1/s^2)*m^2*(1/kg^2) = m^3/(kg*s^2)

# @export variables

# public variables

# private variables

# @onready variables

# Check if this Kepler Object has a child of a certain class, if so, return the first one found
# * @function getChildByClass
# * @param String - Name of the Class to search for.
# * @returns {KEPLER} - Returns a KEPLER object which is a child of this Kepler Object
# * @public
func getChildByClass(classTarget:String)->Node:
    var children:Array = get_children()
    #var childQty:int = get_child_count()
    #var childTarget:OrbitMesh = null
    
    for child in children:
        if child.get_class() == classTarget:
            return child
        else:
            continue
    
    #No OrbitMesh found
    return null
    pass

