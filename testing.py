import trimesh
import numpy as np
from sklearn.cluster import DBSCAN

class Avatar:
    def __init__(self, mesh_path):
        self.mesh = trimesh.load(mesh_path)
        self.v = self.mesh.vertices
        self.scene = trimesh.Scene(self.mesh)
        self.clean_mesh()
        self.align_model()
        self.mark_extrema()
        self.mark_feet()
        self.mark_crotch()
        #self.draw_bounding_box()
        #self.draw_axes()
        self.scene.show()

    def clean_mesh(self):
        self.mesh.remove_duplicate_faces()
        self.mesh.remove_unreferenced_vertices()
        self.mesh.fill_holes()
        self.mesh.update_faces(self.mesh.nondegenerate_faces(height=1e-8))

    def align_model(self):
        bounds = self.mesh.bounds
        extents = bounds[1] - bounds[0]
        major_axis = np.argmax(extents)

        rotation_to_z = self._rotation_matrix_to_axis(major_axis, target_axis=2)
        self.mesh.apply_transform(rotation_to_z)

        center_z = np.mean(self.mesh.bounds[:, 2])
        z_dist = np.abs(self.mesh.vertices[:, 2] - center_z)
        closest_mask = z_dist < 0.01
        cross_section_vertices = self.mesh.vertices[closest_mask]

        lines = []
        errors = []
        for axis in [0, 1]:
            side_center = np.mean(self.mesh.bounds[:, axis])
            point_on_line = np.zeros(3)
            point_on_line[axis] = side_center
            point_on_line[2] = center_z

            direction = np.zeros(3)
            direction[2] = 1

            diffs = cross_section_vertices - point_on_line
            projections = np.outer(diffs @ direction, direction)
            reconstructed = point_on_line + projections
            mse = np.mean(np.linalg.norm(cross_section_vertices - reconstructed, axis=1) ** 2)
            lines.append((point_on_line, direction))
            errors.append(mse)

        x_axis = np.argmin(errors)
        rotation_to_x = self._rotation_matrix_to_axis(x_axis, target_axis=0)
        self.mesh.apply_transform(rotation_to_x)

        center = np.mean(self.mesh.vertices, axis=0)
        self.mesh.apply_translation(-center)

    def _rotation_matrix_to_axis(self, source_axis, target_axis):
        if source_axis == target_axis:
            return np.eye(4)
        axes = np.eye(3)
        v1 = axes[source_axis]
        v2 = axes[target_axis]
        cross = np.cross(v1, v2)
        dot = np.dot(v1, v2)
        skew = np.array([[0, -cross[2], cross[1]], [cross[2], 0, -cross[0]], [-cross[1], cross[0], 0]])
        R = np.eye(3) + skew + skew @ skew * ((1 - dot) / (np.linalg.norm(cross) ** 2))
        M = np.eye(4)
        M[:3, :3] = R
        return M

    def mark_extrema(self):
        bounds = self.mesh.bounds
        colors = [
            [255, 0, 0, 255],
            [0, 0, 255, 255],
            [0, 255, 0, 255]
        ]
        for i in range(3):
            for j, val in enumerate(bounds[:, i]):
                marker = trimesh.creation.icosphere(radius=0.01)
                pos = np.zeros(3)
                pos[i] = val
                marker.apply_translation(pos)
                marker.visual.face_colors = colors[i]
                self.scene.add_geometry(marker)



    def mark_feet(self):
        v = np.asarray(self.mesh.vertices)
        z = v[:, 2]
        lower_mask = z < np.percentile(z, 10)
        lower = v[lower_mask]
        left = lower[lower[:, 0] < 0]
        right = lower[lower[:, 0] > 0]
        left_foot = left[np.argmin(left[:, 2])] if len(left) > 0 else None
        right_foot = right[np.argmin(right[:, 2])] if len(right) > 0 else None
        for foot in [left_foot, right_foot]:
            if foot is not None:
                marker = trimesh.creation.icosphere(radius=0.01)
                marker.apply_translation(foot)
                marker.visual.face_colors = [128, 0, 128, 255]
                self.scene.add_geometry(marker)

    def mark_crotch(self):
        verts = self.mesh.vertices
        x_values = verts[:, 0]
        z_values = verts[:, 2]

        min_x, max_x = np.percentile(x_values, 10), np.percentile(x_values, 90)
        print(f"Min X (25th percentile): {min_x}")
        print(f"Max X (75th percentile): {max_x}")

        # Visualize 25th and 75th percentile x values as vertical lines
        for x_val, color in zip([min_x, max_x], [[255, 0, 0, 255], [0, 0, 255, 255]]):
            line = trimesh.load_path(np.array([
                [[x_val, 0, z_values.min()], [x_val, 0, z_values.max()]]
            ]))
            line.colors = np.tile(color, (len(line.entities), 1))
            self.scene.add_geometry(line)

        # Now try to find the crotch by sweeping through z slices and keeping only middle 50% in X
        z_range = np.linspace(z_values.min(), z_values.max(), 60)
        prev_n_clusters = None

        for z in z_range:
            slice_mask = (
                (np.abs(verts[:, 2] - z) < 0.05) &  # Thickness of the slice
                (verts[:, 0] >= min_x) &
                (verts[:, 0] <= max_x)
            )
            slice_verts = verts[slice_mask]
            if len(slice_verts) < 10:
                continue

            coords_2d = slice_verts[:, :2]
            clustering = DBSCAN(eps=0.03, min_samples=5).fit(coords_2d)
            labels = clustering.labels_
            n_clusters = len(set(labels)) - (1 if -1 in labels else 0)

            if prev_n_clusters == 2 and n_clusters == 1:
                crotch_point = np.mean(slice_verts, axis=0)
                print(f"Crotch found at: {crotch_point}")
                marker = trimesh.creation.icosphere(radius=0.02)
                marker.apply_translation(crotch_point)
                marker.visual.face_colors = [255, 165, 0, 255]
                self.scene.add_geometry(marker)
                break

            prev_n_clusters = n_clusters

    def draw_bounding_box(self):
        corners = self.mesh.bounding_box_oriented.vertices
        for v in corners:
            marker = trimesh.creation.icosphere(radius=0.2)
            marker.apply_translation(v)
            marker.visual.face_colors = [255, 255, 0, 255]
            self.scene.add_geometry(marker)

    def draw_axes(self):
        origin = np.array([0.0, 0.0, 0.0])
        length = 50
        axes = [
            ([length, 0, 0], [255, 0, 0, 255]),
            ([0, length, 0], [0, 0, 255, 255]),
            ([0, 0, length], [0, 255, 0, 255])
        ]
        for vec, color in axes:
            line = trimesh.load_path(np.array([[origin, origin + vec]]))
            line.colors = np.tile(color, (len(line.entities), 1))
            self.scene.add_geometry(line)

mesh_path = "man.obj"
avatar = Avatar(mesh_path)
