<div style="text-align: center;">
  <img src="images/topomox.png" width="350"/>
</div>

# topoMOx: Topologies for Metal Oxides

- [Introduction](#Introduction)
- [Installation](#Installation)
- [Usage](#Usage)
- [How to Cite](#Howtocite)
- [License and Copyright Information](#licenseandcopyrightinformation) 
- [Support and Contact](#supportandcontact)
- [References](#references)

## Introduction
TopoMOx is a tool specifically designed to construct semi-rigid force fields for classical molecular dynamics simulations, currently tailored for the GROMACS simulation program.

The "semi-rigid" model has been widely employed in classical simulations of metal oxides in solution, particularly polyoxometalates (POMs).<sup>1-3</sup>. This model assigns high values to bond and angle constants, which has proven to be a valid approach since most interactions within metal oxides are through van der Waals or electrostatic forces, as defined at the quantum mechanical (QM) level. This validity has been demonstrated by Bonet et al.,<sup>4</sup> Leroy et al.,<sup>5</sup> and Chaumont and Wipff.<sup>6</sup>

**Author**: Albert Masip-Sánchez<sup>a</sup>

<small><sup>a</sup> Quantum Chemistry Group, Physical and Inorganic Chemistry Department, Universitat Rovira i Virgili - Tarragona - Spain</small>

**Developer(s)**: Albert Masip-Sánchez

## Installation
TopoMOx can be installed via git and pip by executing the following commands:
```
git clone https://github.com/qcgurv/topoMOx.git
cd topomox
git checkout v1.0.0
pip install -e .
```

## Usage
Once installed successfully, TopoMOx takes an optimized xyz structure and a file containing partial atomic charges to generate the GRO file and the ffbonded.itp file. These files are ready for use as the topology of the metal oxide in a classical molecular dynamics simulation. Execute it as follows:

```
python -m topomox.main <xyz file> <charges file>
```

In the [Examples](examples) section, you can find a demonstration for a [PW<sub>12</sub>O<sub>40</sub>]<sup>3-</sup> system.

## How to cite
We appreciate your support and hope that our tool has been helpful in your research on metal-oxide clusters.

> Masip-Sánchez, A. topoMOx: Topologies for Metal Oxides, *1.0.0*, **2024**.

## Support and Contact
In case you should encounter problems or bugs, please write a short message to albert.masip@urv.cat.

## References
<sup>1</sup> A. Solé-Daura, V. Goovaerts, K. Stroobants, G. Absillis, P. Jiménez-Lozano, J. M. Poblet, J. D. Hirst, T. N. Parac-Vogt, J. J. Carbó, *Chem. Eur. J.* **2016**, 22, 15280.

<sup>2</sup> A. Solé-Daura, A. Notario-Estévez, J. J. Carbó, J. M. Poblet, C. de Graaf, K. Y. Monakhov, and X. López, *Inorg. Chem.* **2019**, 58, 6, 3881-3894

<sup>3</sup> A. Tzaguy, A. Masip-Sánchez, L. Avram, A. Solé-Daura, X. López, J. M. Poblet, and R. Neumann, *J. Am. Chem. Soc.* **2023**. 145, 36, 19912-19924

<sup>4</sup> X. López, C. Nieto-Draghi, C. Bo, J. Bonet-Avalos, and J. M. Poblet, *J. Phys. Chem. A* **2005**, 109, 6, 1216-1222

<sup>5</sup> F. Leroy, P. Miró, J. M. Poblet, C. Bo, and J. Bonet-Ávalos, *J. Phys. Chem. B.* **2008**, 112, 29, 8591-8599

<sup>6</sup> A. Chaumont, and G. Wipff, *Phys. Chem. Chem. Phys.* **2008**, 10, 6940-6953
