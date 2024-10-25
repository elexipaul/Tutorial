<!--StartFragment-->

**MOLECULAR DYNAMICS SIMULATION USING GROMACS**

**AIM**: Molecular Dynamics Simulation using Gromacs

**MATERIAL/TOOLS:** Linux, GROMACS, GRACE, PyMol

**Introduction:**\
\
GROMACS is one of the most widely used open-source and free software codes in chemistry, used primarily for dynamical simulations of biomolecules. It provides a rich set of calculation types, preparation and analysis tools. Several advanced techniques for free-energy calculations are supported.

Molecular structures, which were solved by X-ray crystallography, nuclear magnetic resonance, or in any other way are stored in the Protein Data Bank. To make those data easily accessible to researchers, the PDB file format had been developed. This type of file contains information about atom names, their positions in the Cartesian coordinate system, bindings to other atoms, and some other auxiliary information such as locations of the secondary structure elements. The PDB formatted files are plain text, so they can be easily read and manipulated by most programs dedicated to molecular structure analysis. They can also be directly read by any text viewing program. 

Molecular dynamics (MD) simulations are powerful computational tools used to explore the structural, dynamic, and thermodynamic properties of biomolecules at an atomic level. In this study, we performed a comprehensive MD simulation of a protein system using GROMACS, a widely used software for simulating the molecular mechanics of proteins, lipids, and other biomolecules. The simulation workflow included preprocessing steps such as cleaning the protein structure, solvating the protein in a water box, ion addition to neutralize the system, energy minimization, and equilibration phases. We then ran a production MD simulation, capturing the dynamic behavior of the protein over time. The results were analyzed to obtain essential properties such as Root Mean Square Deviation (RMSD), Radius of Gyration, and system density over time.

The objectives of this MD simulation were:

1. **To analyze structural stability** of the protein throughout the simulation by calculating RMSD and Radius of Gyration.

2. **To assess system stability** during the equilibration phases (NVT and NPT) by tracking temperature, pressure, and density.

\
\


**PROCEDURE:**

**Downloading protein structure**:

Go to the RCSB website and download the PDB text for the crystal structure. Once the structure has been downloaded, it can be visualized using a viewing program such as VMD, Chimera, PyMOL, etc.

**Deleting water molecules:**

use grep to delete these lines very easily

**Preparing input file for gromacs using pdb2gmx**: Now that the crystal waters are gone and we have verified that all the necessary atoms are present, the PDB file should contain only protein atoms, and is ready to be input into the first GROMACS module, pdb2gmx. The purpose of pdb2gmx is to generate three files:

- The topology for the molecule. 

- A position restraint file.

- A post-processed structure file.

The topology (topol.top by default) contains all the information necessary to define the molecule within a simulation. This information includes nonbonded parameters (atom types and charges) as well as bonded parameters (bonds, angles, and dihedrals). Execute pdb2gmx 

**Generated three new files:**

1. **1AKI processed.gro**: 

1AKI\_processed.gro is a GROMACS-formatted structure file that contains all the atoms defined within the force field (Le, H atoms have been added to the amino acids in the protein)

1. **Topol.top:** The topol.top file is the system topology.

2) **Posre.itp**: The posre itp file contains information used to restrain the positions atoms.

**Defining box for solvation:**

It is possible to simulate proteins and other molecules in different solvents, provided that good parameters are available for all species involved. There are two steps to defining the box and filling it with solvent: Define the box dimensions using the editconf module. Fill the box with water using the solvate module.

**Solvating the box:**

Now that we have defined a box, we can fill it with solvent (water). Solvation is accomplished using solvate 

gmx solvate -cp 1AKI\_newbox gro -cs spc216.gro -o 1AKI\_solv.gro -p topol.top The configuration of the protein (-cp) is contained in the output of the previous editconf step, and the configuration of the solvent (-es) is part of the standard GROMACS installation. We are using spc216.gro, which is a generic equilibrated 3-point solvent model. You can use spe216.gro as the solvent configuration for SPC, SPC/E, or TIP3P water, since they are all three-point water models. The output is called 1AKI solv.gro, and we tell solvate the name of the topology file (topol.top) so it can be modified.\
**Adding ions:**\
Assemble tpr file

Now we have an atomic-level description of our system in the binary file ions.tpr.

**Energy Minimization**:

The solvated, electroneutral system is now assembled. Before we can begin dynamics, we must ensure that the system has no steric clashes or inappropriate geometry. The structure is relaxed through a process called energy minimization (EM). 


### **GROMACS Molecular Dynamics Simulation Protocol**<a id="h.sbjaan185epd"></a>

#### **1. Remove Water Molecules from the PDB File**<a id="h.at5gjjbcsfr6"></a>

grep -v HOH 1aki.pdb > output\_clean.pdb

- 1aki.pdb original PDB file name.

- output\_clean.pdb will be the cleaned PDB file with water molecules removed.

***


#### **2. Convert PDB to GROMACS Format**<a id="h.x28kgau5jezb"></a>

gmx pdb2gmx -f output\_clean.pdb -o processed.gro -water spce

- output\_clean.pdb is the input cleaned PDB file.

- processed.gro is the output file in GROMACS format.

***


#### **3. Define the Simulation Box**<a id="h.paehdm8xmnpw"></a>

gmx editconf -f processed.gro -o newbox.gro -c -d 1.0 -bt cubic

- Adjust -d 1.0 to set the box size (distance in nm between solute and box edge).

- newbox.gro is the output GROMACS file for the box.

***


#### **4. Solvate the System**<a id="h.d50tymmn3g88"></a>

gmx solvate -cp newbox.gro -cs spc216.gro -o solvated.gro -p topol.top

- solvated.gro is the output with solvent molecules added.

***


#### **5. Add Ions**<a id="h.lgqog74xc9oo"></a>

**Prepare Input for Ion Addition:**

gmx grompp -f ions.mdp -c solvated.gro -p topol.top -o ions.tpr

**Neutralize the System by Adding Ions:**

gmx genion -s ions.tpr -o solvated\_ions.gro -p topol.top -pname NA -nname CL -neutral

- solvated\_ions.gro is the output with ions added.

***


#### **6. Energy Minimization**<a id="h.9vt3u5dkocri"></a>

**Prepare for Minimization:**

gmx grompp -f minim.mdp -c solvated\_ions.gro -p topol.top -o em.tpr

**Run Minimization:**

gmx mdrun -v -deffnm em

**Analyze Potential Energy:**

gmx energy -f em.edr -o potential.xvg

***


#### **7. Equilibration**<a id="h.43dw0d34r9gm"></a>

##### **7.1 NVT Equilibration (Constant Volume, Temperature)**<a id="h.gkerd6p5tk4c"></a>

**Prepare NVT Equilibration:**

gmx grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr

**Run NVT Equilibration:**

gmx mdrun -deffnm nvt

**Analyze Temperature Profile:**

gmx energy -f nvt.edr -o temperature.xvg

_(Select option 16 0 when prompted)_

***


##### **7.2 NPT Equilibration (Constant Pressure, Temperature)**<a id="h.undcvc5p30r3"></a>

**Prepare NPT Equilibration:**

gmx grompp -f npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -o npt.tpr

**Run NPT Equilibration:**

gmx mdrun -deffnm npt

**Analyze Pressure and Density Profiles:**

gmx energy -f npt.edr -o pressure.xvg

gmx energy -f npt.edr -o density.xvg

_(Select options 18 0 for pressure and 24 0 for density)_

***


#### **8. Production Molecular Dynamics Run**<a id="h.9m4b0mcsk83j"></a>

**Prepare Production MD Run:**

gmx grompp -f md.mdp -c npt.gro -t npt.cpt -p topol.top -o md\_0\_1.tpr

**Run Production MD (GPU Enabled):**

gmx mdrun -deffnm md\_0\_1 -nb gpu

***


#### **9. Post-Processing and Analysis**<a id="h.wt0f43exgqon"></a>

**Remove Periodic Boundary Conditions:**

gmx trjconv -s md\_0\_1.tpr -f md\_0\_1.xtc -o md\_0\_1\_noPBC.xtc -pbc mol -center

**Calculate RMSD of Trajectory:**

gmx rms -s md\_0\_1.tpr -f md\_0\_1\_noPBC.xtc -o rmsd.xvg -tu ns

gmx rms -s em.tpr -f md\_0\_1\_noPBC.xtc -o rmsd\_xtal.xvg -tu ns

**Calculate Radius of Gyration:**

gmx gyrate -s md\_0\_1.tpr -f md\_0\_1\_noPBC.xtc -o gyrate.xvg

\
\
\
\
\
\
\
\
\
\
\
\
\
\
\


**Results Summary**

1\. Solvation Box Setup

\
\
\
\
![](https://lh7-rt.googleusercontent.com/docsz/AD_4nXcEsaFJsiEzjlfp-lnCLOrlCeZrFwI1J0SNrDMxwn4Ir5ZKB1hCNG9RCggdSgrO_4yT0BK_6d9yv_pFdMAdeA7aU8J7abLoBC8ggOXCYjUa6P7h-Gy00rHi0Kjm6z1N4TlOxiDMXFG034-JqmVrd7riyRA?key=ClqgZUQ9fz0lZaeWbU0Flg)\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\


Figure: Solvation Box with Protein and Water Molecules

The protein was placed in a cubic simulation box with a 1 nm buffer distance from the edges, as shown in Figure X (Solvation Box). This buffer space provides enough room for solvent interactions around the protein without excessive computational overhead. After defining the box, we added water molecules using the SPC/E water model, resulting in a fully solvated system ready for ionization and minimization steps.

\
\
\
\
\
\
\
\
\
\


2\. Energy Minimization

![](https://lh7-rt.googleusercontent.com/docsz/AD_4nXcZhdFkFpKbB095a0BDTxSyrcz3OiyC0NYGCEKAujWPvSTUINSXiyOvq8o4SayGwjodCpASPEZ2t46pT1vHLYKPBz-dTFrWB0SWnrIUbhYQTALxzFEYr6lS0xRqxX9Cje3gdgrdFr8P2HJaeJXOpZV1Y7Q?key=ClqgZUQ9fz0lZaeWbU0Flg)

Figure: Potential Energy of the System During Minimization

The energy minimization step is essential to remove steric clashes and unfavorable atomic overlaps in the system. As shown in Figure Y (Potential Energy), the potential energy decreased significantly over the minimization steps, converging to a stable value. This reduction confirms that the system has achieved a local energy minimum, indicating readiness for the equilibration steps.

3\. Equilibration Phase

3.1 NVT Equilibration (Constant Volume, Temperature)

\
\
\
\
\
\


![](https://lh7-rt.googleusercontent.com/docsz/AD_4nXdudOQHVvq-Ngg2ca33SiUYRGVcCcM0f-OKm12TCm1IEEPfLNiyClVYNmQL1VriWfgO6oFFCFFcncmAauwMQvUq0X1ixB05FHjTMtr4PEuNBkBi698S4iyD0o3T6f7gKOxSufY6REimwdZK1WXxMeYiteA?key=ClqgZUQ9fz0lZaeWbU0Flg)

Figure: Temperature Profile During NVT Equilibration

During the NVT equilibration phase, the system’s temperature was stabilized to the target of 300 K, as shown in Figure Z (Temperature Profile). The plot demonstrates a stable temperature achieved after an initial fluctuation period, confirming that the system reached thermal equilibrium under constant volume conditions.

3.2 NPT Equilibration (Constant Pressure, Temperature)

\
\
\
\
\
\
\
\
\
\
\


![](https://lh7-rt.googleusercontent.com/docsz/AD_4nXcOyy4oIsPwy5h7eGCoekIMli52akxX1VyiZfo6XLZYmaUCsPeZIf9m3iO6deLaR0vWdbKVV3XgNRsAIyNgXjrn_6s1dM2qanJ1vuVdUNQ_OGbcmTJIHaWA-MQXL2exYX0c1FS5b45m-zHOgzQLYjzAbA?key=ClqgZUQ9fz0lZaeWbU0Flg)![](https://lh7-rt.googleusercontent.com/docsz/AD_4nXcqfDxZFI1TJfJbP1RJJQDLgD5-jTq5wvgR3lOQGxbDMdpX0OU9uuoHvngF_qW8sAbaWMRAO9kZOCHlx8YGumWB3UETC4WUJotD_h3jehi5H3n7eJ4qATW_p0hLmi21X-pA8Z07R4NqLLnAGYkqxzckVA?key=ClqgZUQ9fz0lZaeWbU0Flg)

Figure: Pressure and Density Profiles During NPT Equilibration

In the NPT phase, the system was equilibrated under constant pressure to reach a density close to experimental water density. Figure A (Pressure Profile) and Figure B (Density Profile). The density profile confirms successful equilibration, ensuring system stability before the production MD run.

4\. Production MD Simulation

4.1 Structural Stability: RMSD Analysis

\


![](https://lh7-rt.googleusercontent.com/docsz/AD_4nXeb72hzO9PRZd1kwizvWxXhBQ1iGtXzvAJeSgiRmTukmC0wiOsrwotwsxZLOqiNA1yp6RlBNvBpGUvRMljWq9s9b-x55IfvIE5X8CnpuS-u-csAvGRyCu6uKQPzpw5gd3BIfKZJMkRfaM_PA0j8XlNKfA?key=ClqgZUQ9fz0lZaeWbU0Flg)

Figure: RMSD of Protein Backbone Over Time

To assess the structural stability of the protein, we calculated the RMSD of the backbone atoms throughout the MD simulation, as shown in Figure C (RMSD Profile). The RMSD plot indicates initial structural adaptation followed by stabilization, suggesting that the protein maintains a stable conformation during the simulation, with fluctuations within an acceptable range.

4.2 Compactness: Radius of Gyration

![](https://lh7-rt.googleusercontent.com/docsz/AD_4nXciyGvdWfQfXRWAD9Y1C3MgHDvEixK_HZQ78esZwOUJl8Wn0vVQeJHxTMkfITu9YE-_QIq5jKZSSi6D5vzNMEDPOpcfBcK9Yc97T4eljb01EMW__R1ArhZX0D4wmQPz8cqdKBh20hKRtNJh0qyXVjTrMA?key=ClqgZUQ9fz0lZaeWbU0Flg)

Figure: Radius of Gyration Over Time

The Radius of Gyration (Rg) provides insight into the compactness and folding stability of the protein


### **4.3RMSF Analysis**<a id="h.7o4br1p5tcnu"></a>

Root Mean Square Fluctuation (RMSF) helps in analyzing the flexibility of residues in a protein by showing their average deviation from a reference position. RMSF values are useful for identifying regions of the protein that are more flexible or stable throughout the simulation.

![](https://lh7-rt.googleusercontent.com/docsz/AD_4nXeTfxutzvpqHuAYPgO2xvQUuQTg194DDXgKHHNsIQ2cwDtRd7vB32hVZFfICA3Yob03DMBPL-ht6qNzeBfmSW7R1OL0w3BziNP7luSloQ9T_In32mtLrNcaEIzcLnP-Dtjt5J4zXZhwjQZ0tRMvtUdPVzc?key=ClqgZUQ9fz0lZaeWbU0Flg)


### **4.4. Hydrogen Bond Analysis**<a id="h.tautk19am771"></a>

Hydrogen bonds are key in stabilizing secondary and tertiary structures in proteins. This analysis helps in tracking the number and strength of hydrogen bonds over time in your system.

gmx hbond computes and analyzes hydrogen bonds. Hydrogen bonds are determined based on cutoffs for the angle Hydrogen - Donor- Acceptor (zero is extended) and the distance Donor- Acceptor (or Hydrogen - Acceptor using -noda) OH and NH groups are regarded as donors, O is an acceptor always, N is an acceptor by default, but this can be switched using nitace. Dummy hydrogen atoms are assumed to be connected to the first preceding non-hydrogen atom.

\
\
\
\
\
\
\
\
\


![](https://lh7-rt.googleusercontent.com/docsz/AD_4nXfDxQeNU0XaDsyD55zelLhRrt8QGX1VrRK_gg0nM6D_geBIZd7yo7xDJwE_4-tqJ7YJnijfKWmfqolC-ucFa9Q38_lxXVNr3DK-Wm0fL4TXMk_z6TyXOmWjgC_kHLJDY5XJVm6bAIYSIwOdfAsSt2PgfiM?key=ClqgZUQ9fz0lZaeWbU0Flg)

\
\
\
\
\
\
\
\
\
\
\
\
\
\
\


5\. Solvent and Ion Distribution

![](https://lh7-rt.googleusercontent.com/docsz/AD_4nXe1QiU-8IRwdKWlPk9qZHNig6mI0JoglAsZZByuohKASJvLQaVTAbGFfsbtTHVlOoKA1QjYjgiyK6ACNBOgpyNJC_I1GOJwyPUdK4_XgtGV9-vJJSHbJUsqcIT3p9C7jGbWFmC_KBoXPox5JbATemTAz1c?key=ClqgZUQ9fz0lZaeWbU0Flg)

. Figure E (Ion Distribution in Solvation Box) shows a uniform distribution of ions around the protein, ensuring the system's neutrality and stability. This step was essential for preventing artifacts in electrostatic interactions during MD simulation.

\
\
\
\
\
\
\
\


6\. Conclusion

The MD simulation successfully captured the dynamic behavior and stability of the protein within a solvated and neutralized environment. The stable temperature, pressure, and density profiles confirm the adequate equilibration of the system. At the same time, RMSD and Radius of Gyration analyses indicate a stable and compact protein structure throughout the simulation. Future studies may delve deeper into the conformational transitions observed in this system to further understand their relevance to the protein’s biological functions.

\
\


<!--EndFragment-->
