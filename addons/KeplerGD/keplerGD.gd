@tool
extends EditorPlugin


func _enter_tree() -> void:
    # Initialization of the plugin goes here.
    #add_custom_type("Orbit Node", "Node", preload("orbit.gd"), preload("orbit_icon.png"))
    #add_custom_type("Orbit Mesh", "Node", preload("orbit.gd"), preload("orbit_mesh_icon.png"))
    #add_custom_type("AstroBody Node", "Node", preload("astro_body.gd"), preload("astro_body_icon.png"))
    #add_custom_type("Astro Mesh", "Node", preload("astro_mesh.gd"), preload("astro_mesh_icon.png"))
    
    add_custom_type("Orbit", "Kepler", preload("orbit.gd"), preload("orbit_icon.png"))
    add_custom_type("OrbitMesh", "MeshInstance3D", preload("orbit_mesh.gd"), preload("orbit_mesh_icon.png"))
    add_custom_type("AstroBody", "Kepler", preload("astro_body.gd"), preload("astro_body_icon.png"))
    add_custom_type("AstroMesh", "MeshInstance3D", preload("astro_mesh.gd"), preload("astro_mesh_icon.png"))

func _exit_tree() -> void:
    # Clean-up of the plugin goes here.
    remove_custom_type("Orbit Node")
    remove_custom_type("Orbit Mesh")
    remove_custom_type("AstroBody Node")
    remove_custom_type("Astro Mesh")
