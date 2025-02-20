import open3d
import open3d.visualization
import os
from pathlib import Path

class Settings:
    pass

class Measurement:
    pass
    
class Head(Measurement):
    pass
    
class Torso(Measurement):
    pass 

class Leg(Measurement):
    pass

class Arm(Measurement):
    pass

class Avatar:    
    def __init__(self, fp: os.PathLike):
        self.path = Path(fp)
        self.mesh = None
        
        if self.path.suffix in ('.ply', '.obj'):
            self.mesh = open3d.io.read_triangle_mesh(self.path)
    
    def render(self):
        open3d.visualization.draw_geometries([self.mesh])

if __name__ == "__main__":
    obj = Avatar("test/mesh/cow.ply")
    obj.render()
        