import trimesh
import numpy as np

from visualize import *
from body import *

def visualize_test():
    mesh_file_path = 'test/mesh/penn-mesh-2.obj' 
    body_mesh = trimesh.load(mesh_file_path)
    volume = Body(body_mesh)
    
    show_volume(volume)
    
    
if __name__ == '__main__':
    show_all_meshes()