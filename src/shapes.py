import numpy as np
import trimesh
from shapely.geometry import Polygon


def convexity_search(mesh: trimesh.Trimesh, 
                    rays: int,
                    origin=np.array([0, 0, 0]),
                    max_iter=20,
                    slice_width=.005,
                    convexity_threshold=.05):
    
    z = origin[2]
    
    thetas = np.linspace(0, 2*np.pi, rays)
    ray_directions = np.empty((rays, 3))
    ray_directions[:, 0] = np.cos(thetas)
    ray_directions[:, 1] = np.sin(thetas)
    
    ray_origins = np.empty((rays, 3))
    ray_origins[:, :] = origin
    
    convexity_scores = []
    for iteration in range(max_iter):
        ray_directions[:, 2] = iteration * slice_width + z
        ray_origins[:, 2] = iteration * slice_width + z
        
        intersections = mesh.ray.intersects_location(
            ray_origins=ray_origins,
            ray_directions=ray_directions
        )[0]
        
        center = np.mean(intersections, axis=0)
        assigned = assign_points_to_rays(intersections, ray_directions, center)
        outer_section = pick_outermost_points(assigned, center)
        loops = order_loop_from_points(outer_section)
        
        convexity_score = compute_convexity(loops)
        if convexity_score >= 1 - convexity_threshold:
            return center
        convexity_scores.append((convexity_score, center))
    
    convexity_scores.sort(key=lambda x: x, reverse=True)
    return convexity_scores[0][1]
        
def assign_points_to_rays(points, directions, centroid):
    """
    Assign each point to the ray direction it aligns with most (in angle).
    """
    assigned = [[] for _ in range(len(directions))]
    vectors = points - centroid  # shift to geometric center
    
    for p, v in zip(points, vectors):
        unit_v = v / np.linalg.norm(v)
        dots = directions @ unit_v
        idx = np.argmax(dots)  # most aligned ray
        assigned[idx].append(p)

    return assigned

def pick_outermost_points(assigned, centroid):
    """
    Pick the farthest point from centroid per ray direction.
    """
    result = []
    for hits in assigned:
        if not hits:
            continue
        hits = np.array(hits)
        dists = np.linalg.norm(hits - centroid, axis=1)
        result.append(hits[np.argmax(dists)])
    return np.array(result)

def order_loop_from_points(points):
    """
    Order points around their centroid by angle in XY plane.
    """
    centroid = points[:, :2].mean(axis=0)
    rel = points[:, :2] - centroid
    angles = np.arctan2(rel[:, 1], rel[:, 0])
    order = np.argsort(angles)
    return points[order]

def compute_convexity(loop):
    loop_2d = loop[:, :2]
    poly = Polygon(loop_2d)

    if not poly.is_valid or poly.area == 0:
        return 0, 0, 0

    area = poly.area
    hull = poly.convex_hull
    hull_area = hull.area
    convexity = area / hull_area

    return convexity