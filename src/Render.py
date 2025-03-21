import open3d
import open3d.visualization
import os
from pathlib import Path
import pymeshlab

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

class Render:    
    def __init__(self, fp: os.PathLike):
        self.path = Path(fp)
        self.mesh = None
        
        if self.path.suffix == '.obj':
            obj = pymeshlab.MeshSet()
            obj.load_new_mesh(str(self.path))
            self.path = self.path.parent / Path(str(self.path.stem) + '.ply')
            obj.save_current_mesh(str(self.path))
            
        if self.path.suffix == '.ply':
            self.mesh = open3d.io.read_triangle_mesh(str(self.path))

        
        
    def render(self):
        open3d.visualization.draw_geometries([self.mesh])

if __name__ == "__main__":
    obj = Render("test/mesh/cow.ply")
    obj.render()
    