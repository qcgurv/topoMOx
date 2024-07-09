import re
from .constants import atomic_weights, thresholds, hetero

def get_threshold(element, thresholds, hetero):
    if element in hetero:
        return thresholds['XO']
    if element == 'O':
        return None
    if element in thresholds:
        return thresholds[element]
    else:
        return 1.60, 2.99

def modify_element(row, identified_metals, hetero):
    if row['element'] == 'H':
        return 'Hp'
    elif row['element'] in identified_metals:
        return f'{row["element"]}P' if len(row['element']) == 1 else row['element']
    elif row['element'] in hetero:
        return row['element']
    else:
        return f'{row["element"]}P' if len(row['element']) == 1 else row['element']

def get_atomic_weight(row):
    element = row['atom'].strip()
    if len(element) > 1 and element[1].islower():
        element = element[:2]
    else:
        element = element[0]
    return atomic_weights[element]

def extract_atom_number(atom_name):
    return int(re.findall(r'\d+', atom_name)[0])
