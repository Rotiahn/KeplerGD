@tool
class_name OrbitMesh
extends MeshInstance3D


# Called when the node enters the scene tree for the first time.
func _ready() -> void:
    #draw_orbit()
    pass

# Called every frame. 'delta' is the elapsed time since the previous frame.
func _process(delta: float) -> void:
    pass
  
func get_class()->String:
    return "OrbitMesh"
  

func draw_orbit()->void:
    var orbit:Orbit = get_parent()
    var immediate_mesh:ImmediateMesh = ImmediateMesh.new()
    mesh = immediate_mesh
    var material:StandardMaterial3D = StandardMaterial3D.new()
    material.shading_mode = BaseMaterial3D.SHADING_MODE_UNSHADED
    var astroParent:AstroBody = orbit.get_parent()
    var astromesh:AstroMesh = astroParent.getChildByClass("AstroMesh")
    if (astromesh == null):  #We don't have an AstroMesh, add one
        astromesh = AstroMesh.new()
        astroParent.add_child(astromesh)
    var astromeshmesh:Mesh = astromesh.mesh
    if (astromeshmesh == null): #AstroMesh doesn't have a mesh, add one
        astromeshmesh = SphereMesh.new()
        astromesh.mesh = astromeshmesh
    var astromaterial:StandardMaterial3D = astromeshmesh.surface_get_material(0)
    if (astromaterial == null): #mesh doesn't have a material, add one
        astromaterial = StandardMaterial3D.new()
        astromeshmesh.surface_set_material(0,astromaterial)
        astromaterial.albedo_color = Color.DEEP_SKY_BLUE
    material.albedo_color = astromaterial.albedo_color
    
    

    
    #Begin draw
    immediate_mesh.surface_begin(Mesh.PRIMITIVE_LINE_STRIP,material)
    
    #Iterate through mAnomaly drawing a point every 1 degrees
    var calc_orbit:Orbit = orbit.clone()
    for i in range(361):
        calc_orbit.mAnomaly = i * Kepler.DEGREE
        var new_point:Vector3 = calc_orbit.getPosition()
        #print(self.name,new_point)
        immediate_mesh.surface_add_vertex(new_point)
    
    immediate_mesh.surface_end()
