from simtk.openmm import app
import simtk.openmm as mm
from simtk import unit
from sys import stdout

from pdbfixer import PDBFixer

from mdtraj.reporters import HDF5Reporter
import mdtraj as mdt
import ast
import os
import re
import numpy as np
import pandas as pd


# general parameters
boxtype = 'rectangular'
box_padding = 1.0# nanometers
production_steps = 500000# with timestep 2fs

# GPU parameters
gpu_index = '0'

# force fields
sim_force_field = 'amber14-all.xml'
sim_water_model = 'amber14/tip3p.xml'

prot = 'HG3_H2O.pdb'

protonation_dict = {}


sim_ph = 7

sim_temperature = 308.15

trajectory_out_atoms = 'protein'
restrained_eq_atoms = 'protein and name CA'

df = pd.read_csv("kemp_shuffle_variants_pdb_mut.csv", sep=";")
traj_folder = 'minimized_shuffle_structures'
os.system('mkdir {0}'.format(traj_folder))
for number, row in df.iterrows():
    variant_id = row["Variant"]
    result=None
    variant = row['Variant']
    mutations = row['PDB_mutations'].replace('[', '').replace(']', '')
    pattern = re.compile(r'(?<!\'A\'),')
    mutation_list = re.split(pattern, mutations)
    mutation_list = [ast.literal_eval(s) for s in mutation_list]    
    xml_list = [sim_force_field, sim_water_model]
    forcefield = app.ForceField(*xml_list)
    template = PDBFixer(filename=prot)
    for mutation in mutation_list:
        replacement = mutation[1]
        chain = mutation[0]
        template.applyMutations([replacement,], chain)
    template.findNonstandardResidues()
    template.findMissingResidues()
    template.findMissingAtoms()
    template.addMissingAtoms()    
    app.PDBFile.writeFile(template.topology, template.positions, open('variant.pdb', 'w'), keepIds=True)
    protein_pdb = app.PDBFile('variant.pdb')
    protein_pdb.topology.createStandardBonds()
    protein_mod = app.Modeller(protein_pdb.topology, protein_pdb.positions)
    protonation_list = []
    if len(protonation_dict.keys()) > 0:
    	for chain in protein_mod.topology.chains():
            chain_id = chain.id
            
            protonations_in_chain_dict = {}
            for protonation_tuple in protonation_dict:
            	if chain_id == protonation_tuple[0]:
            		residue_number = protonation_tuple[1]
            		protonations_in_chain_dict[residue_number] = protonation_dict[protonation_tuple]
            
            for residue in chain.residues():
            	residue_id = residue.id
            	if residue_id in protonations_in_chain_dict.keys():
            		protonation_list.append(protonations_in_chain_dict[residue_id])
            	else:
            		protonation_list.append(None)   
    protein_mod.addHydrogens(forcefield, pH=sim_ph)


    x_list = []
    y_list = []
    z_list = []

    for index in range(len(protein_mod.positions)):
    	x_list.append(protein_mod.positions[index][0]._value)
    	y_list.append(protein_mod.positions[index][1]._value)
    	z_list.append(protein_mod.positions[index][2]._value)

    x_span = (max(x_list) - min(x_list))
    y_span = (max(y_list) - min(y_list))
    z_span = (max(z_list) - min(z_list))


    d =  max(x_span, y_span, z_span) + (2 * box_padding)

    d_x = x_span + (2 * box_padding)
    d_y = y_span + (2 * box_padding)
    d_z = z_span + (2 * box_padding)

    prot_x_mid = min(x_list) + (0.5 * x_span)
    prot_y_mid = min(y_list) + (0.5 * y_span)
    prot_z_mid = min(z_list) + (0.5 * z_span)

    box_x_mid = d_x * 0.5
    box_y_mid = d_y * 0.5
    box_z_mid = d_z * 0.5

    shift_x = box_x_mid - prot_x_mid
    shift_y = box_y_mid - prot_y_mid
    shift_z = box_z_mid - prot_z_mid

    solvated_model = app.Modeller(protein_mod.topology, protein_mod.positions)

    for index in range(len(solvated_model.positions)):
    	solvated_model.positions[index] = (solvated_model.positions[index][0]._value + shift_x, solvated_model.positions[index][1]._value + shift_y, solvated_model.positions[index][2]._value + shift_z)*unit.nanometers
    while result is None:
        try:
            if boxtype == 'cubic':
            	solvated_model.addSolvent(forcefield, model='tip3p', neutralize=True, ionicStrength=0*unit.molar, boxVectors=(mm.Vec3(d, 0., 0.), mm.Vec3(0., d, 0.), mm.Vec3(0, 0, d)))
            elif boxtype == 'rectangular':
            	solvated_model.addSolvent(forcefield, model='tip3p', neutralize=True, ionicStrength=0*unit.molar, boxVectors=(mm.Vec3(d_x, 0., 0.), mm.Vec3(0., d_y, 0.), mm.Vec3(0, 0, d_z)))
            result=True
        except:
            pass
    selection_reference_topology = mdt.Topology().from_openmm(solvated_model.topology)
    trajectory_out_indices = selection_reference_topology.select(trajectory_out_atoms)
    restrained_eq_indices = selection_reference_topology.select(restrained_eq_atoms)
    system = forcefield.createSystem(solvated_model.topology, nonbondedMethod=app.PME, nonbondedCutoff=1.0*unit.nanometers, constraints=app.HBonds, ewaldErrorTolerance=0.0005, rigidWater=True)

    integrator = mm.LangevinIntegrator(sim_temperature*unit.kelvin, 1/unit.picosecond, 0.002*unit.picoseconds)

    platform = mm.Platform.getPlatformByName('CUDA')
    properties = {'Precision': 'single', 'DeviceIndex': gpu_index}

    simulation = app.Simulation(solvated_model.topology, system, integrator, platform, properties)
    simulation.context.setPositions(solvated_model.positions)
    print(f'Minimizing {variant_id}')
    simulation.minimizeEnergy(tolerance=0.1 * unit.kilojoules_per_mole / unit.nanometer, maxIterations=0)
    min_pos = simulation.context.getState(getPositions=True).getPositions()
    app.PDBFile.writeFile(simulation.topology, min_pos, open(traj_folder + '/' + f'MIN_{variant_id}.pdb', 'w'))
    del(simulation)	