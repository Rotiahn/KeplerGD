@tool
class_name AstroMesh
extends MeshInstance3D

# Called when the node enters the scene tree for the first time.
func _ready() -> void:
    var parent:AstroBody = get_parent()
    position = parent.getPosition() ##Vector3
    
    pass # Replace with function body.


# Called every frame. 'delta' is the elapsed time since the previous frame.
func _process(delta: float) -> void:
    var parent:AstroBody = get_parent()
    position = parent.getPosition() ##Vector3
    #print(name," ",parent.primary," ",position)
    
    pass

func get_class()->String:
    return "AstroMesh"
