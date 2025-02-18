import open3d
import open3d.visualization

#
mesh = open3d.io.read_triangle_mesh("test/mesh/cow.ply")

# Draws mesh to screen
open3d.visualization.draw_geometries([mesh])

print(mesh)