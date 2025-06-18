import trimesh
import numpy as np
from scipy.spatial import cKDTree
from typing import Tuple
from shapes import convexity_search
# Conventions from negative to positive
# x is left-right
# y is back-front
# z is down-up

class Volume:
    def __init__(self, mesh: trimesh.Trimesh):
        mesh: trimesh.Trimesh = self.clean_mesh(mesh)
        
        std = np.std(mesh.vertices)
        mean = np.mean(mesh.vertices)
        
        self._normalize = np.vectorize(lambda x: (x-mean) / std )
        self._denormalize = np.vectorize(lambda x: (x*std + mean))
        
        self.mesh = self.normalize(mesh)
        
    def clean_mesh(self, mesh: trimesh.Trimesh) -> trimesh.Trimesh:
        cleaned_mesh = mesh.copy()
        cleaned_mesh.update_faces(cleaned_mesh.unique_faces())
        cleaned_mesh.remove_unreferenced_vertices()
        cleaned_mesh.fill_holes()
        cleaned_mesh.update_faces(cleaned_mesh.nondegenerate_faces(height=1e-8))
        
        return cleaned_mesh
    
    def normalize(self, mesh: trimesh.Trimesh) -> trimesh.Trimesh:
        new_mesh = mesh.copy()
        new_mesh.vertices = self._normalize(new_mesh.vertices)
        return new_mesh 
    
    def denormalize(self, mesh: trimesh.Trimesh) -> trimesh.Trimesh:
        new_mesh = mesh.copy()
        new_mesh.vertices = self._denormalize(new_mesh.vertices)
        return new_mesh
    

class Body(Volume):
    def __init__(self, mesh: trimesh.Trimesh):
        super().__init__(mesh)
        self.mesh = self.orient(self.mesh)
        
        
    def orient(self, mesh: trimesh.Trimesh) -> trimesh.Trimesh:
        new_mesh = mesh.copy()
        
        # Center mesh
        midpoint = (np.max(new_mesh.vertices, axis=0) + np.min(new_mesh.vertices, axis=0)) / 2
        new_mesh.vertices -= midpoint
        
        # Find height
        major_axis = self.find_major_axis(new_mesh)
        if major_axis == 0:
            new_mesh.vertices[:, (0, 2)] = new_mesh.vertices[:, (2, 0)] 
        elif major_axis == 1:
            new_mesh.vertices[:, (1, 2)] = new_mesh.vertices[:, (2, 1)]
        
        # Determine left-right axis
        # This algorithm creates a vertical line at the minimum value at both the x
        # and y axes, then finds the nearest neighbors to those lines. Whichever
        # result leads to the lowest total distance between the nearest neighbors and
        # a vertical line at the midpoint is taken to be the y-axis (that is to say,
        # front to back)
                
        kdtree = cKDTree(new_mesh.vertices)
        
        top_of_head = new_mesh.vertices[np.argmax(new_mesh.vertices, axis=0)[2]]
        z_axis = np.linspace(0, top_of_head[2], 100)
        
        x_coord, y_coord = np.min(new_mesh.vertices, axis=0)[:2]
        x_vec = np.empty((len(z_axis), 3))
        y_vec = np.empty((len(z_axis), 3))
        
        x_vec[:, (0, 1)] = (x_coord, 0)
        x_vec[:, 2] = z_axis
        
        y_vec[:, (0, 1)] = (0, y_coord)
        y_vec[:, 2] = z_axis
        
        x_alignment = kdtree.query(x_vec)[1]
        y_alignment = kdtree.query(y_vec)[1]
        
        x_vertices = new_mesh.vertices[x_alignment]
        y_vertices = new_mesh.vertices[y_alignment]
        
        x_score = np.linalg.norm(x_vertices[:2], axis=0).sum()
        y_score = np.linalg.norm(y_vertices[:2], axis=0).sum()
       
        if y_score > x_score:
            print("Switching x and y axis...")
            new_mesh.vertices[:, (0, 1)] = new_mesh.vertices[:, (1, 0)]
                
        return new_mesh
        
    def find_major_axis(self, mesh: trimesh.Trimesh) -> int:
        mins = np.min(mesh.vertices, axis=0)
        maxs = np.max(mesh.vertices, axis=0)
        major_axis = np.argmax(maxs - mins)
        return major_axis
    
    def _identify_feet(self, mesh: trimesh.Trimesh) -> Tuple[np.ndarray, np.ndarray]:
        # This method assumes the mesh is oriented
        new_mesh = mesh.copy()
        
        kdtree = cKDTree(new_mesh.vertices)
        
        lower_half_vertices = new_mesh.vertices[new_mesh.vertices[:, 2] < 0, :]
        left_side = lower_half_vertices[lower_half_vertices[:, 0] < 0, :]
        left_foot = left_side[np.argmin(left_side, axis=0)[2], :]
        
        right_foot_estimate = left_foot.copy()
        right_foot_estimate[0] *= -1
        
        right_foot = new_mesh.vertices[kdtree.query(right_foot_estimate)[1], :]
        
        return left_foot, right_foot
        
    def _identify_crotch(self, 
                         mesh: trimesh.Trimesh) -> np.ndarray:
        new_mesh = mesh.copy()
        
        kdtree = cKDTree(new_mesh.vertices)
        
        minimum_z = new_mesh.vertices[np.argmin(new_mesh.vertices, axis=0)[2], 2]
        
        ray_origin = np.array([0, 0, minimum_z])
        ray_direction = np.array([0, 0, 1])
        
        intersects = new_mesh.ray.intersects_location(
            ray_origins=[ray_origin],
            ray_directions=[ray_direction]
        )[0]
        
        min_viable_point = intersects[np.argmin(intersects, axis=0)[2], :]
        viable_point_idx =  kdtree.query(min_viable_point)[1]
        viable_point = new_mesh.vertices[viable_point_idx]
        
        crotch_point_nearest = convexity_search(new_mesh, 
                                        rays=32,
                                        origin=viable_point)
        
        crotch_point_idx = kdtree.query(crotch_point_nearest)[1]
        crotch_point = new_mesh.vertices[crotch_point_idx]
        
        return crotch_point
        
if __name__ == '__main__':
    mesh_file_path = 'test/mesh/penn-mesh-2.obj' 
    body_mesh = trimesh.load_mesh(mesh_file_path)
    wrong_volume = Volume(body_mesh)
    