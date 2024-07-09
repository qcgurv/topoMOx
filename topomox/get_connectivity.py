import numpy as np

def euclidean_distance(a, b):
    return np.sqrt((a['x'] - b['x']) ** 2 + (a['y'] - b['y']) ** 2 + (a['z'] - b['z']) ** 2)

def angle_between_atoms(atom_a, atom_b, atom_c):
    vec_ab = np.array([atom_a['x'] - atom_b['x'], atom_a['y'] - atom_b['y'], atom_a['z'] - atom_b['z']])
    vec_cb = np.array([atom_c['x'] - atom_b['x'], atom_c['y'] - atom_b['y'], atom_c['z'] - atom_b['z']])
    dot_product = np.dot(vec_ab, vec_cb)
    norm_ab = np.linalg.norm(vec_ab)
    norm_cb = np.linalg.norm(vec_cb)
    angle_radians = np.arccos(dot_product / (norm_ab * norm_cb))
    angle_degrees = np.degrees(angle_radians)
    return angle_degrees
