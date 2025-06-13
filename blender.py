import bpy,random,mathutils, os

# === Settings ===
xyz_path = "/Users/vsathyas/Downloads/p_1.xyz"
output_path = "/Users/vsathyas/Downloads/rendered_image.png"
scale_factor = 1.0
radii_scale = 3.0
material_cache = {}

# === UFF radii ===
UFF_RADII = {
    'H': 0.354, 'He': 0.849, 'Li': 1.336, 'Be': 1.074, 'B': 0.838, 'C': 0.757,
    'N': 0.700, 'O': 0.658, 'F': 0.668, 'Ne': 0.920, 'Na': 1.539, 'Mg': 1.421,
    'Al': 1.244, 'Si': 1.117, 'P': 1.117, 'S': 1.064, 'Cl': 1.044, 'Ar': 1.032,
    'K': 1.953, 'Ca': 1.761, 'Sc': 1.513, 'Ti': 1.412, 'V': 1.402, 'Cr': 1.345,
    'Mn': 1.382, 'Fe': 1.335, 'Co': 1.241, 'Ni': 1.164, 'Cu': 1.302, 'Zn': 1.193,
    'Ga': 1.260, 'Ge': 1.197, 'As': 1.211, 'Se': 1.190, 'Br': 1.192, 'Kr': 1.147,
    'Rb': 2.260, 'Sr': 2.052, 'Y': 1.698, 'Zr': 1.564, 'Nb': 1.473, 'Mo': 1.484,
    'Tc': 1.322, 'Ru': 1.478, 'Rh': 1.332, 'Pd': 1.338, 'Ag': 1.386, 'Cd': 1.403,
    'In': 1.459, 'Sn': 1.398, 'Sb': 1.407, 'Te': 1.386, 'I': 1.382, 'Xe': 1.267,
    'Cs': 2.570, 'Ba': 2.277, 'La': 1.943, 'Hf': 1.611, 'Ta': 1.511, 'W': 1.526,
    'Re': 1.372, 'Os': 1.372, 'Ir': 1.371, 'Pt': 1.364, 'Au': 1.262, 'Hg': 1.340,
    'Tl': 1.518, 'Pb': 1.459, 'Bi': 1.512, 'Po': 1.500, 'At': 1.545, 'Rn': 1.420,
    'default': 0.7
}

# === Define Element Colors (RGBA) ===
element_colors = {
    "H":  (1.00, 1.00, 1.00, 1),   # White
    "C":  (0.20, 0.20, 0.20, 1),   # Dark Gray
    "O":  (1.00, 0.00, 0.00, 1),   # Red
    "N":  (0.00, 0.00, 1.00, 1),   # Blue
    "F":  (0.00, 1.00, 0.00, 1),   # Green
    "Cl": (0.12, 0.94, 0.12, 1),   # Bright Green
    "Br": (0.65, 0.16, 0.16, 1),   # Dark Red
    "I":  (0.58, 0.00, 0.83, 1),   # Purple
    "S":  (1.00, 1.00, 0.18, 1),   # Yellow
    "P":  (1.00, 0.50, 0.00, 1),   # Orange
    "B":  (1.00, 0.70, 0.70, 1),   # Pale Pink
    "Li": (0.80, 0.50, 1.00, 1),   # Light Purple
    "Na": (0.67, 0.36, 0.95, 1),   # Violet
    "Mg": (0.54, 1.00, 0.00, 1),   # Light Green
    "Al": (0.75, 0.65, 0.65, 1),   # Light Gray
    "Si": (0.94, 0.78, 0.63, 1),   # Beige
    "K":  (0.56, 0.25, 0.83, 1),   # Purple
    "Ca": (0.24, 1.00, 0.00, 1),   # Greenish Yellow
    "Ti": (0.75, 0.76, 0.78, 1),   # Silver
    "Cr": (0.54, 0.60, 0.78, 1),   # Blue Gray
    "Mn": (0.61, 0.48, 0.78, 1),   # Violet Gray
    "Fe": (0.87, 0.40, 0.20, 1),   # Rust
    "Co": (0.94, 0.56, 0.63, 1),   # Pinkish Red
    "Ni": (0.31, 0.82, 0.31, 1),   # Mint
    "Cu": (1.00, 0.49, 0.00, 1),   # Copper Orange
    "Zn": (0.49, 0.50, 0.69, 1),   # Slate Blue
    "As": (0.50, 0.38, 0.37, 1),   # Brown
    "Se": (1.00, 0.63, 0.00, 1),   # Orange
    "Ag": (0.75, 0.75, 0.75, 1),   # Silver
    "Au": (1.00, 0.84, 0.00, 1),   # Gold
    "Hg": (0.72, 0.72, 0.82, 1),   # Pale Lavender
    "Pb": (0.34, 0.35, 0.36, 1),   # Dark Gray
}

def get_material(element, color,type='metal'):
    
    if type =='metal':
        if element not in material_cache:
            mat = bpy.data.materials.new(name=f"{element}_mat")
            mat.use_nodes = True
            bsdf = mat.node_tree.nodes["Principled BSDF"]
            bsdf.inputs["Base Color"].default_value = color
            bsdf.inputs["Roughness"].default_value = 1.0
            bsdf.inputs["Metallic"].default_value = 1.0
            
            material_cache[element] = mat
            return mat
        
    if type == 'glass':
        if element not in material_cache:
            mat = bpy.data.materials.new(name=f"{element}_mat")
            mat.use_nodes = True
            nodes = mat.node_tree.nodes
            links = mat.node_tree.links

            # Clear default nodes
            for n in nodes:
                nodes.remove(n)

            # Add Glass BSDF
            output = nodes.new(type="ShaderNodeOutputMaterial")
            glass = nodes.new(type="ShaderNodeBsdfGlass")

            glass.inputs["Color"].default_value = color  # RGBA
            glass.inputs["Roughness"].default_value = 0.05  # Optional tweak
            glass.inputs["IOR"].default_value = 1.45  # Index of refraction, e.g. ~1.45 for quartz

            links.new(glass.outputs["BSDF"], output.inputs["Surface"])

            mat.blend_method = 'BLEND'         # For Eevee: allow transparency
            mat.use_screen_refraction = True   # Eevee only
            mat.show_transparent_back = False
            mat.use_backface_culling = False
            mat.shadow_method = 'HASHED'       # Eevee shadows
            mat.refraction_depth = 0.1
            material_cache[element] = mat
            return mat
    return material_cache[element]

def read_xyz(xyz_path):
    atoms = dict()
    with open(xyz_path, 'r') as f:
        lines = f.readlines()[2:]  # Skip count and comment
        for count_line,line in enumerate(lines):
            tokens = line.strip().split()
            if len(tokens) < 4:
                continue
            
            atoms[count_line] = {'element':tokens[0],'position':tuple(float(x) for x in tokens[1:4]),'radii':UFF_RADII.get(tokens[0], UFF_RADII['default'])}
    
    return atoms

# === Clear Scene ===
bpy.ops.object.select_all(action='SELECT')
bpy.ops.object.delete(use_global=False)

# === Load XYZ ===

atoms = read_xyz(xyz_path)
atom_objs = []

for atom in atoms:
    
    e,x,y,z,r = atoms[atom]['element'],*atoms[atom]["position"],atoms[atom]["radii"]
    bpy.ops.mesh.primitive_uv_sphere_add(radius=r, location=(x * scale_factor, y * scale_factor, z * scale_factor))
    obj = bpy.context.object
    obj.name = e
    obj.data.materials.append(get_material(e,element_colors[e]))
    atom_objs.append(obj)

# === Camera Setup ===
cam_data = bpy.data.cameras.new("Camera")
cam = bpy.data.objects.new("Camera", cam_data)
bpy.context.collection.objects.link(cam)
bpy.context.scene.camera = cam

# === Frame Atoms ===
bbox_min = mathutils.Vector((min(v[i] for o in atom_objs for v in o.bound_box) for i in range(3)))
bbox_max = mathutils.Vector((max(v[i] for o in atom_objs for v in o.bound_box) for i in range(3)))
center = (bbox_min + bbox_max) / 2
size = (bbox_max - bbox_min).length

cam.location = center + mathutils.Vector((0, -2.0 * size, 0.8 * size))
cam.rotation_euler = (center - cam.location).to_track_quat('-Z', 'Y').to_euler()
cam.data.lens = 100

# === Render Settings ===
scene = bpy.context.scene
scene.render.engine = 'CYCLES'
scene.cycles.device = 'GPU'  # Or 'CPU'
scene.render.resolution_x = 1920
scene.render.resolution_y = 1080
scene.render.image_settings.file_format = 'PNG'
scene.render.filepath = output_path
os.makedirs(os.path.dirname(output_path), exist_ok=True)

# === Render ===
bpy.ops.render.render(write_still=True)
print(f"Image saved to: {output_path}")
