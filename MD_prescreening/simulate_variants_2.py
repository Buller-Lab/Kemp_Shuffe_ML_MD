from simtk.openmm import app
import simtk.openmm as mm
from simtk import unit
from sys import stdout

import parmed as pmd

from pdbfixer import PDBFixer

from mdtraj.reporters import HDF5Reporter
import mdtraj as mdt
import ast
import os
import re
import numpy as np
import pandas as pd



###########################################################
###########################################################
# general parameters
boxtype = 'rectangular'
box_padding = 1.0# nanometers
production_steps = 500000# with timestep 2fs

# GPU parameters
gpu_index = '1'

# force fields
sim_force_field = 'amber14-all.xml'
sim_water_model = 'amber14/tip3p.xml'
sim_gaff = 'gaff.xml'

# protein structure
prot = 'HG3_H2O.pdb'

# ligand parameters
ligand_names = ["H5J"]
ligand_xml_files = ['H5J.xml']
ligand_pdb_files = ["H5J.pdb"]

# additional residue definitions (for ligands or novel residues)
additional_residue_definitions_file = "add_residue_def.xml"

protonation_dict = {}#('A', 226): 'ASH', ('A', 195): 'GLH'}

# pH
sim_ph = 7

# temperature
sim_temperature = 308.15

# trajectory
trajectory_out_atoms = 'protein or resname H5J or chainid 1'# MDTraj selection syntax
trajectory_out_interval = 500# write to trajectory every ... steps

# restrained atoms during equilibration
restrained_eq_atoms = 'protein and name CA or resname H5J or chainid 1'# MDTraj selection syntax

###########################################################
###########################################################

df = pd.read_csv("kemp_shuffle_variants_pdb_mut_2.csv", sep=";")

import subprocess

def execute_bash_command(command, output_file):
    # Run the bash command and capture the output
    result = subprocess.run(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

    # Write the output to a text file
    with open(output_file, 'w') as file:
        file.write(result.stdout)

    # Print the last row of the file
    with open(output_file, 'r') as file:
        lines = file.readlines()
        if lines:
            print("Free Energy of Unfolding:")
            print(lines[-3].strip().split("=")[1])
            row["deltaG_stability"]=lines[-3].strip().split("=")[1]
        else:
            print("The output file is empty.")

for index, row in df.iterrows():
    result=None
    variant = row['Variant']
    mutations = row['PDB_mutations'].replace('[', '').replace(']', '')
    pattern = re.compile(r'(?<!\'A\'),')
    mutation_list = re.split(pattern, mutations)
    mutation_list = [ast.literal_eval(s) for s in mutation_list]
    traj_folder = f'MD_{variant}'
    os.system('mkdir {0}'.format(traj_folder))
    xml_list = [sim_force_field, sim_water_model, sim_gaff]
    for lig_xml_file in ligand_xml_files:
    	xml_list.append(lig_xml_file)
    forcefield = app.ForceField(*xml_list)########### hier unpacking der Liste xml_list mit dem Aterisken (*), weil Funktion ForceField() keine Liste akzeptiert
    template = PDBFixer(filename=prot)
    if len(mutation_list) > 0:
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
    if additional_residue_definitions_file != None:
    	protein_pdb.topology.loadBondDefinitions(additional_residue_definitions_file)
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
    protein_mod.addHydrogens(forcefield, pH=sim_ph)#, variants = protonation_list)

    # add ligand structures to the model
    for lig_pdb_file in ligand_pdb_files:
    	ligand_pdb = app.PDBFile(lig_pdb_file)
    	protein_mod.add(ligand_pdb.topology, ligand_pdb.positions)

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
    for i in range(20):
        system = forcefield.createSystem(solvated_model.topology, nonbondedMethod=app.PME, nonbondedCutoff=1.0*unit.nanometers, constraints=app.HBonds, ewaldErrorTolerance=0.0005, rigidWater=True)

        integrator = mm.LangevinIntegrator(sim_temperature*unit.kelvin, 1/unit.picosecond, 0.002*unit.picoseconds)

        platform = mm.Platform.getPlatformByName('CUDA')
        properties = {'Precision': 'single', 'DeviceIndex': gpu_index}

        simulation = app.Simulation(solvated_model.topology, system, integrator, platform, properties)
        simulation.context.setPositions(solvated_model.positions)
        print('Minimizing')
        simulation.minimizeEnergy()
        min_pos = simulation.context.getState(getPositions=True).getPositions()
        app.PDBFile.writeFile(simulation.topology, min_pos, open(traj_folder + '/' + f'MIN_{i}.pdb', 'w'))
        del(simulation)	
        force = mm.CustomExternalForce("(k/2)*periodicdistance(x, y, z, x0, y0, z0)^2")
        force.addGlobalParameter("k", 100.0*unit.kilojoules_per_mole/unit.angstroms**2)
        force.addPerParticleParameter("x0")
        force.addPerParticleParameter("y0")
        force.addPerParticleParameter("z0")
        for res_atom_index in restrained_eq_indices:
        	force.addParticle(int(res_atom_index), min_pos[int(res_atom_index)].value_in_unit(unit.nanometers))
        system.addForce(force)
        integrator = mm.LangevinIntegrator(sim_temperature*unit.kelvin, 1/unit.picosecond, 0.002*unit.picoseconds)
        simulation = app.Simulation(solvated_model.topology, system, integrator, platform, properties)
        simulation.context.setPositions(min_pos)
        simulation.context.setVelocitiesToTemperature(sim_temperature*unit.kelvin)
        simulation.reporters.append(app.StateDataReporter(stdout, 1000, step=True, potentialEnergy=True, temperature=True, progress=True, remainingTime=True, speed=True, totalSteps=40000, separator='\t'))
        simulation.reporters.append(HDF5Reporter(traj_folder + '/' + f'EQ_NVT_{i}.h5', 10000, atomSubset=trajectory_out_indices))
        print('restrained NVT equilibration (ligand and backbone)')
        simulation.step(40000)
        state_nvt_EQ = simulation.context.getState(getPositions=True, getVelocities=True)
        positions = state_nvt_EQ.getPositions()
        app.PDBFile.writeFile(simulation.topology, positions, open(traj_folder + '/' + f'NVT_EQ_{i}.pdb', 'w'), keepIds=True)
        del(simulation)  
        system.addForce(mm.MonteCarloBarostat(1*unit.atmospheres, sim_temperature*unit.kelvin, 25))
        integrator = mm.LangevinIntegrator(sim_temperature*unit.kelvin, 1/unit.picosecond, 0.002*unit.picoseconds)
        simulation = app.Simulation(solvated_model.topology, system, integrator, platform, properties)
        simulation.context.setState(state_nvt_EQ)
        simulation.reporters.append(app.StateDataReporter(stdout, 1000, step=True, potentialEnergy=True, temperature=True, progress=True, remainingTime=True, speed=True, totalSteps=production_steps, separator='\t'))
        simulation.reporters.append(HDF5Reporter(traj_folder + '/' + f'EQ_NPT_{i}.h5', 10000, atomSubset=trajectory_out_indices))
        print('restrained NPT equilibration (ligand and backbone)')
        simulation.step(production_steps)
        state_npt_EQ = simulation.context.getState(getPositions=True, getVelocities=True)
        positions = state_npt_EQ.getPositions()
        app.PDBFile.writeFile(simulation.topology, positions, open(traj_folder + '/' + f'NPT_EQ_{i}.pdb', 'w'), keepIds=True)
        del(simulation)
        n_forces = len(system.getForces())
        system.removeForce(n_forces-2)
        integrator = mm.LangevinIntegrator(sim_temperature*unit.kelvin, 1/unit.picosecond, 0.002*unit.picoseconds)
        simulation = app.Simulation(solvated_model.topology, system, integrator, platform, properties)
        simulation.context.setState(state_npt_EQ)
        simulation.reporters.append(app.StateDataReporter(stdout, 1000, step=True, potentialEnergy=True, temperature=True, progress=True, remainingTime=True, speed=True, totalSteps=production_steps, separator='\t'))
        simulation.reporters.append(HDF5Reporter(traj_folder + '/' + f'EQ_NPT_free_{i}.h5', 10000, atomSubset=trajectory_out_indices))
        print('free NPT equilibration')
        simulation.step(production_steps)
        state_free_EQ = simulation.context.getState(getPositions=True, getVelocities=True)
        positions = state_free_EQ.getPositions()
        app.PDBFile.writeFile(simulation.topology, positions, open(traj_folder + '/' + f'free_NPT_EQ_{i}.pdb', 'w'), keepIds=True)
        del(simulation)
        integrator = mm.LangevinIntegrator(sim_temperature*unit.kelvin, 1/unit.picosecond, 0.002*unit.picoseconds)
        simulation = app.Simulation(solvated_model.topology, system, integrator, platform, properties)
        simulation.context.setState(state_free_EQ)
        simulation.reporters.append(app.StateDataReporter(stdout, 10000, step=True, potentialEnergy=True, temperature=True, progress=True, remainingTime=True, speed=True, totalSteps=production_steps, separator='\t'))
        simulation.reporters.append(HDF5Reporter(traj_folder + '/' + f'production_{i}.h5', trajectory_out_interval, atomSubset=trajectory_out_indices))
        simulation.step(production_steps)
        state_production = simulation.context.getState(getPositions=True, getVelocities=True)
        positions = state_production.getPositions()			
        app.PDBFile.writeFile(simulation.topology, positions, open(traj_folder + '/' + f'production_{i}.pdb', 'w'), keepIds=True)
        del(simulation)		
        print('production run')
    