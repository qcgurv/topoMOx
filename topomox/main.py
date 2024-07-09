import sys
import pandas as pd
from itertools import combinations
import os
import re
import numpy as np

from .constants import atomic_weights, thresholds, hetero
from .utilities import get_threshold, modify_element, get_atomic_weight, extract_atom_number
from .get_connectivity import euclidean_distance, angle_between_atoms
from .gro_generation import write_gro, format_row, format_bonds_row, format_angles_row
from .dependency_installer import install_and_import


def main():
    if len(sys.argv) != 3:
        print("Usage: topoMOx.py <input file in XYZ format> <ChelpG charges file>")
        sys.exit(1)
    
    infile = None
    charges = None
    for arg in sys.argv[1:]:
        if arg.endswith(".xyz"):
            infile = arg
        else:
            charges = arg

    if infile is None or charges is None:
        print("Error: Both an XYZ input file and a ChelpG charges file are required.")
        sys.exit(1)

    resname = infile.split('.')[0][:3].upper()
    xyz = pd.read_csv(infile, names=['atom', 'x', 'y', 'z'], skiprows=2, delim_whitespace=True, dtype={'x': float, 'y': float, 'z': float})
    xyz['element'] = xyz['atom'].apply(lambda x: re.match(r'^[A-Za-z]+', x).group(0))
    xyz['atom'] = xyz.apply(lambda row: f"{row['element']}{row.name + 1}", axis=1)
    xyz = xyz.round({'x': 5, 'y': 5, 'z': 5})
    
    identified_metals = set(xyz['element']).intersection(set(thresholds.keys()))

    for metal in set(xyz['element']):
        if metal not in identified_metals:
            threshold = get_threshold(metal, thresholds, hetero)
            thresholds[metal] = threshold

    bonds = []
    all_metals = set(xyz['element']).intersection(set(atomic_weights.keys())) - set(hetero) - set(['O'])
    
    for _, atom in xyz[xyz['element'].isin(all_metals | set(hetero))].iterrows():
        for _, oxygen_atom in xyz[xyz['element'] == 'O'].iterrows():
            distance = euclidean_distance(atom, oxygen_atom)
            if atom['element'] in identified_metals:
                min_distance, max_distance = thresholds[atom['element']]
                if min_distance <= distance <= max_distance:
                    bonds.append((atom['atom'], oxygen_atom['atom'], distance))
            elif atom['element'] in hetero:
                min_distance, max_distance = thresholds['XO']
                if min_distance <= distance <= max_distance:
                    bonds.append((atom['atom'], oxygen_atom['atom'], distance))
            elif atom['element'] in all_metals - identified_metals:
                min_distance, max_distance = (1.60, 2.99)
                if min_distance <= distance <= max_distance:
                    bonds.append((atom['atom'], oxygen_atom['atom'], distance))
    
    bonds_df = pd.DataFrame(bonds, columns=['Atom', 'Oxygen Atom', 'Distance'])
    bonded_atoms = {}
    for _, row in bonds_df.iterrows():
        if row['Atom'] not in bonded_atoms:
            bonded_atoms[row['Atom']] = [row['Oxygen Atom']]
        else:
            bonded_atoms[row['Atom']].append(row['Oxygen Atom'])
        if row['Oxygen Atom'] not in bonded_atoms:
            bonded_atoms[row['Oxygen Atom']] = [row['Atom']]
        else:
            bonded_atoms[row['Oxygen Atom']].append(row['Atom'])

    angles = []
    for atom2 in bonded_atoms.keys():
        bonded_to_atom2 = bonded_atoms[atom2]
        for atom1, atom3 in combinations(bonded_to_atom2, 2):
            angle = angle_between_atoms(xyz.loc[xyz['atom'] == atom1].iloc[0], 
                                        xyz.loc[xyz['atom'] == atom2].iloc[0], 
                                        xyz.loc[xyz['atom'] == atom3].iloc[0])
            angles.append((atom1, atom2, atom3, angle))
    
    angles_df = pd.DataFrame(angles, columns=['Atom A', 'Atom B', 'Atom C', 'Angle'])
    charge_df = pd.read_csv(charges, usecols=[2], names=['charge'], delim_whitespace=True)
    
    xyz['attype'] = xyz.apply(lambda row: modify_element(row, identified_metals, hetero), axis=1)
    xyz['mass'] = xyz.apply(get_atomic_weight, axis=1)
    
    new_columns = {
        'nr': np.arange(1, len(xyz) + 1),
        'resid': np.ones(len(xyz), dtype=int),
        'resname': [resname] * len(xyz),
        'cg': np.arange(1, len(xyz) + 1)
    }
    
    molec_df = pd.concat([pd.DataFrame(new_columns), xyz['attype'], xyz['atom'], xyz['mass'], charge_df], axis=1)
    molec_df = molec_df[['nr', 'attype', 'resid', 'resname', 'atom', 'cg', 'charge', 'mass']]
    
    bond_params_df = pd.DataFrame()
    bond_params_df['i'] = bonds_df['Atom'].apply(extract_atom_number)
    bond_params_df['j'] = bonds_df['Oxygen Atom'].apply(extract_atom_number)
    bond_params_df['n'] = 1
    bond_params_df['distance'] = bonds_df['Distance'] / 10
    bond_params_df['kb'] = 900000.0
    
    angle_params_df = pd.DataFrame()
    angle_params_df['i'] = angles_df['Atom A'].apply(extract_atom_number)
    angle_params_df['j'] = angles_df['Atom B'].apply(extract_atom_number)
    angle_params_df['k'] = angles_df['Atom C'].apply(extract_atom_number)
    angle_params_df['n'] = 1
    angle_params_df['angle'] = angles_df['Angle']
    angle_params_df['k_theta'] = 900.0
    
    header = ";  nr   attype  resid    resname  atom     cg      charge         mass"
    formatted_rows = molec_df.apply(lambda row: format_row(row, identified_metals, hetero), axis=1)
    formatted_rows = formatted_rows.dropna()
    formatted_output = "\n".join(["[ atoms ]", header] + list(formatted_rows))
    
    with open("molec", "w") as f:
        f.write(formatted_output)
    
    bonds_header = ";    i     j   n   distance      kb"
    formatted_bonds_rows = bond_params_df.apply(lambda row: format_bonds_row(row, identified_metals, hetero), axis=1)
    formatted_bonds_output = "\n".join(["[ bonds ]", bonds_header] + list(formatted_bonds_rows))
    
    with open("bonds", "w") as f:
        f.write(formatted_bonds_output)
    
    angles_header = ";   i      j       k    n   angle     k_theta"
    formatted_angles_rows = angle_params_df.apply(lambda row: format_angles_row(row, identified_metals, hetero), axis=1)
    formatted_angles_output = "\n".join(["[ angles ]", angles_header] + list(formatted_angles_rows))
    
    with open("angles", "w") as f:
        f.write(formatted_angles_output)
    
    header = f"[ moleculetype ]\n;  name    nrexcl\n{resname:14}3"
    with open("ffbonded.itp", "w") as ffbonded:
        ffbonded.write(header)
        ffbonded.write("\n\n")
        with open("molec", "r") as molec_file:
            molec_content = molec_file.read()
            ffbonded.write(molec_content)
        ffbonded.write("\n\n")
        with open("bonds", "r") as bonds_file:
            bonds_content = bonds_file.read()
            ffbonded.write(bonds_content)
        ffbonded.write("\n\n")
        with open("angles", "r") as angles_file:
            angles_content = angles_file.read()
            ffbonded.write(angles_content)
    
    for file in ['angles', 'bonds', 'molec']:
        try:
            os.remove(file)
        except FileNotFoundError:
            print(f"File '{file}' not found.")
    
    gro_output_filename = f"{resname}.gro"
    atoms = [(row['atom'], row['x'], row['y'], row['z']) for _, row in xyz.iterrows()]
    write_gro(atoms, resname, gro_output_filename)
    
    ending = f"""
        ┌─────────────────────────────────────┐
        │         Thank you for using         │
        │              topoMOx!               │
        │                                     │
        │  We appreciate your support and hope│
        │  that our tool has been helpful in  │
        │   your research on metal-oxide      │
        │         clusters.                   │
        │                                     │
        │   Please acknowledge the author and │
        │    tool in any resulting work.      │
        │                                     │
        │    Masip-Sánchez, A. topoMOx:       │
        │    Topologies for Metal Oxides,     │
        │              1.0.0, 2024            │
        │                                     │
        │       Albert Masip-Sánchez          │
        │        albert.masip@urv.cat         │
        │   Universitat Rovira i Virgili      │
        └─────────────────────────────────────┘
    """
    
    print(ending)

if __name__ == "__main__":
    dependencies = {
        'numpy': 'np',
        'pandas': 'pd',
        'itertools': 'itertools',
        'os': 'os',
        're': 're',
    }
    for package, import_name in dependencies.items():
        install_and_import(package, import_name)
        
    main()
