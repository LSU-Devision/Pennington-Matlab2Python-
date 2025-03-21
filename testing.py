import trimesh
import numpy as np


class Avatar:
    def __init__(self, mesh_path):
        self.mesh = trimesh.load(mesh_path)
        self.v = self.mesh.vertices  # Get the vertex coordinates
        self.mesh.remove_duplicate_faces()
        self.mesh.remove_unreferenced_vertices()
        self.mesh.fill_holes()
        self.mesh.remove_degenerate_faces()
        self.fix_orientation()

    def fix_orientation(self):
        # Step 1: Compute the mean of the vertices (centroid of the model)
        centroid = np.mean(self.v, axis=0)
        
        # Step 2: Center the model around the origin (subtract the centroid)
        centered_vertices = self.v - centroid
        
        # Step 3: Compute the covariance matrix
        cov_matrix = np.cov(centered_vertices.T)
        
        # Step 4: Perform eigenvalue decomposition to get principal components
        eigvals, eigvecs = np.linalg.eig(cov_matrix)
        
        # Step 5: The principal components are the eigenvectors corresponding to the eigenvalues
        # Sort the eigenvectors by the eigenvalues in descending order
        sorted_indices = np.argsort(eigvals)[::-1]
        sorted_eigvecs = eigvecs[:, sorted_indices]
        
        # Step 6: Get the rotation matrix (the first column of eigenvectors)
        rotation_matrix = sorted_eigvecs.T
        
        # Step 7: Apply the rotation to the mesh vertices
        rotated_vertices = np.dot(centered_vertices, rotation_matrix.T)

        # Apply the rotation to the mesh
        self.mesh.vertices = rotated_vertices + centroid

        # Step 8: Translate the model so that the lowest point (feet) is at the origin (or desired position)
        min_z = np.min(self.mesh.vertices[:, 2])  # Find the lowest Z value (typically the feet)
        translation_vector = np.array([0, 0, -min_z])  # Move the model so the feet are at the origin
        self.mesh.apply_translation(translation_vector)

        self.v = self.mesh.vertices

    def get_legs_min(self):
        v1, v2, v3 = self.v[:, 0], self.v[:, 1], self.v[:, 2]

        # Find median height (to separate upper & lower body)
        median_z = np.median(v3)

        # Filter out upper body (keep only lower half)
        lower_body_mask = v3 < median_z
        lower_v1 = v1[lower_body_mask]
        lower_v2 = v2[lower_body_mask]
        lower_v3 = v3[lower_body_mask]

        # Find lowest point on the left leg
        left_mask = lower_v1 < 0  # Left side vertices
        left_side = np.column_stack((lower_v1[left_mask], lower_v2[left_mask], lower_v3[left_mask]))
        left_leg = left_side[np.argmin(left_side[:, 2])] if len(left_side) > 0 else None

        # Find lowest point on the right leg
        right_mask = lower_v1 > 0  # Right side vertices
        right_side = np.column_stack((lower_v1[right_mask], lower_v2[right_mask], lower_v3[right_mask]))
        right_leg = right_side[np.argmin(right_side[:, 2])] if len(right_side) > 0 else None

        return left_leg, right_leg
    
    def visualize_legs(self):
        l_leg, r_leg = self.get_legs_min()

        # Convert to NumPy array to avoid issues
        points = np.array([l_leg, r_leg])

        # Define colors (Red)
        colors = np.array([[255, 0, 0], [255, 0, 0]])  # Red in RGB

        # Create a point cloud for visibility
        point_cloud = trimesh.points.PointCloud(points, colors=colors)

        # Create a scene and show
        scene = trimesh.Scene([self.mesh, point_cloud])
        scene.show()

# Example usage
mesh_path = "man.obj"  # Replace with your actual 3D model file
avatar = Avatar(mesh_path)
left_leg, right_leg = avatar.get_legs_min()

print("Left Leg:", left_leg)
print("Right Leg:", right_leg)


avatar.visualize_legs()