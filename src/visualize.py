import trimesh
import trimesh.scene
from body import *
from pathlib import Path
import os

def create_point_geometry(point, radius=0.01, color=[255, 0, 0, 255]):
    """
    Create a trimesh.Trimesh object representing a 3D point as a small colored sphere.

    Parameters:
        point (array-like): 3D coordinates of the point (x, y, z).
        radius (float): Radius of the sphere representing the point.
        color (list): RGBA color for the sphere mesh.

    Returns:
        trimesh.Trimesh: A mesh representing the point.
    """
    point = np.asarray(point)

    # Create the sphere
    sphere = trimesh.creation.icosphere(radius=radius)
    sphere.apply_translation(point)

    # Apply color
    sphere.visual.face_colors = color

    return sphere

def create_arrow(origin, direction, length=1.0, shaft_radius=0.01, head_length=0.2, head_radius=0.02, color=[255, 0, 0, 255]):
    """
    Create an arrow mesh (cylinder + cone) representing a vector from origin in the given direction.
    """
    direction = direction / np.linalg.norm(direction)
    shaft_length = length * (1 - head_length)

    # Shaft (cylinder)
    shaft = trimesh.creation.cylinder(radius=shaft_radius, height=shaft_length, sections=16)
    shaft.apply_translation([0, 0, shaft_length / 2])

    # Head (cone)
    head = trimesh.creation.cone(radius=head_radius, height=length * head_length, sections=16)
    head.apply_translation([0, 0, shaft_length + (length * head_length / 2)])

    # Combine
    arrow = shaft.union(head)

    # Align to direction
    z_axis = np.array([0, 0, 1])
    rotation = trimesh.geometry.align_vectors(z_axis, direction)
    arrow.apply_transform(rotation)

    # Move to origin
    arrow.apply_translation(origin)

    # Color
    arrow.visual.face_colors = color

    return arrow

def show_volume(volume: Volume):
    scene = trimesh.Scene(geometry=volume.mesh)
    camera = scene.camera
    
    transformation_matrix = camera.look_at(volume.mesh.vertices,
                                            distance=-10)
    
    x_axis_vector = create_arrow(origin=(0,0,0),
                                 direction=(1,0,0),
                                 length=.5,
                                 shaft_radius=.005,
                                 head_length=.1,
                                 head_radius=.01,
                                 color=(255, 0, 0, 255))
    
    y_axis_vector = create_arrow(origin=(0,0,0),
                                 direction=(0,1,0),
                                 length=.5,
                                 shaft_radius=.005,
                                 head_length=.1,
                                 head_radius=.01,
                                 color=(0, 255, 0, 255))
    
    z_axis_vector = create_arrow(origin=(0,0,0),
                                 direction=(0,0,1),
                                 color=(0, 0, 255, 255))
    
    feet_color = (255, 0, 0, 255)
    left_foot_point, right_foot_point = volume._identify_feet(volume.mesh)
    left_foot = create_point_geometry(left_foot_point, color=feet_color)
    right_foot = create_point_geometry(right_foot_point, color=feet_color)

    crotch_color = (0, 255, 0, 255)
    crotch_point = volume._identify_crotch(volume.mesh)
    crotch = create_point_geometry(crotch_point, color=crotch_color)
    
    scene.add_geometry(left_foot)
    scene.add_geometry(right_foot)
    
    scene.add_geometry(crotch)
    
    scene.add_geometry(x_axis_vector)
    scene.add_geometry(y_axis_vector)
    scene.add_geometry(z_axis_vector)
    
    scene.apply_transform(transformation_matrix)        
    scene.show()
    
def show_all_meshes():
    mesh_dir = Path('test/mesh')
    filename_list = list(map(lambda x: Path(x), os.listdir(mesh_dir)))
    filename_list = list(filter(lambda x: x.suffix == '.obj', filename_list))
    fp_list = [mesh_dir / fp for fp in filename_list]
    
    for fp in fp_list:
        mesh = trimesh.load_mesh(fp)
        volume = Body(mesh)
        show_volume(volume)