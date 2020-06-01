import os, sys
import pandas as pd
import glob

from modules import chimeraplip
from modules import config
from modules import detection
from modules import plipremote
from modules import mp
from modules import plipxml
from modules.preparation import *
from modules import report
from modules import supplemental
from modules import webservices
from modules import visualize
from modules import pymolplip

import gc

import unittest
import shutil

try: # for openbabel < 3.0.0
    import pybel
    
except ImportError: # for openbabel >= 3.0.0
    from openbabel import pybel

pd.set_option('display.max_columns', 20)
pd.set_option('display.max_row', 20)

tmpmol = PDBComplex()

pdbs = glob.glob('D:\\task\\adamNet2.0\\code_an\\PLIP\\Plibpdb\\PDB2\\*.pdb')
names = [os.path.splitext(filename)[0] for filename in os.listdir('D:\\task\\adamNet2.0\\code_an\\PLIP\\Plibpdb\\PDB2')]
output_path = 'D:\\task\\adamNet2.0\\code_an\\PLIP\\Plibpdb\\PDB_output2'
for fname, ids in zip(pdbs, names):
    try:

        tmpmol.load_pdb(fname)

        bsids = []

        for i in tmpmol.__str__()[37:].split('\n'):
            bsids.append(i)

        bsids = sorted(bsids)

        filepath = os.path.join(output_path, ids)

        if bsids != ['']:

            if not os.path.isdir(filepath):
                os.makedirs(filepath)
                pass

            for i in range(len(bsids)):
                filepath_b = os.path.join(output_path, ids, bsids[i].replace(':', ','))

                if not os.path.isdir(filepath_b):
                    os.makedirs(filepath_b)
                    pass

            for ligand in tmpmol.ligands:
                if ':'.join([ligand.hetid, ligand.chain, str(ligand.position)]) in bsids:
                    tmpmol.characterize_complex(ligand)

            for i in range(len(bsids)):

                print('pdb ', fname[-8:-4])
                print('residue_n ', bsids[i])

                s = tmpmol.interaction_sets[bsids[i]]

                # hydrophobic interaction
                # hydrophobics = hydroph = {hydroph.resnr for hydroph in s.hydrophobic_contacts}

                # hydrophobic interaction
                hydrophobics = {hydrophobic.resnr for hydrophobic in s.all_hydrophobic_contacts} #

                # Hydrogen bonds
                bonds_l = {str(hbond.resnr)+hbond.reschain for hbond in s.hbonds_ldon} #
                bonds_p = {str(hbond.resnr)+hbond.reschain for hbond in s.hbonds_pdon} #

                # Salt bridge
                saltb_l = {saltbridge.resnr for saltbridge in s.saltbridge_lneg} #
                saltb_p = {saltbridge.resnr for saltbridge in s.saltbridge_pneg} #

                # Aromatic stacking
                pistackres = {pistack.resnr for pistack in s.pistacking} #

                # Halogen Bonding
                halogens = {halogen.resnr for halogen in s.halogen_bonds} #

                # Water bridges
                waterbridges = {str(wb.resnr)+wb.reschain for wb in s.water_bridges} #

                # Pi-Cation interaction
                picat_l = {pication.resnr for pication in s.pication_laro} #
                picat_p = {pication.resnr for pication in s.pication_paro} #

                # Coordination
                metalres = {mres.restype for mres in s.metal_complexes} #

                """
                print('hydrophobics : ', hydrophobics) #
                print('bonds_l : ', bonds_l) #
                print('bonds_p : ', bonds_p) #  
                print('saltb_l : ', saltb_l) #
                print('saltb_p : ', saltb_p) #
                print('pistackres : ', pistackres) #
                print('halogens : ', halogens) #
                print('waterbridges : ', waterbridges) #
                print('picat_l : ', picat_l) #
                print('picat_p : ', picat_p) #
                print('metalres : ', metalres) #
                print(sep='/n')
                """

                # hydrophobic interactions
                data_hi = pd.DataFrame()

                if len(s.hydrophobic_contacts) != 0:

                    for k in range(len(s.hydrophobic_contacts)):

                        data_hi_temp = []

                        data_hi_temp.append(s.hydrophobic_contacts[k][6])
                        data_hi_temp.append(s.hydrophobic_contacts[k][1])
                        data_hi_temp.append(s.hydrophobic_contacts[k][3])
                        data_hi_temp.append(s.hydrophobic_contacts[k][4])    
                        data_hi_temp.append(s.hydrophobic_contacts[k][5])
                        data_hi_temp.append(s.hydrophobic_contacts[k][7])
                        data_hi_temp.append(s.hydrophobic_contacts[k][8])
                        data_hi_temp.append(s.hydrophobic_contacts[k][9])
                        data_hi_temp.append(s.hydrophobic_contacts[k][10])

                        data_hi_temp = pd.DataFrame(data_hi_temp).T
                        data_hi = data_hi.append(data_hi_temp)

                    data_hi.columns = ['resnr','bsatom_orig_idx','ligatom_orig_idx','distance','restype','reschain','restype_l','resnr_l','reschain_l']
                    data_hi = data_hi.sort_values(by='resnr',ascending=True)

                    if len(data_hi) != 0:
                        data_hi.to_csv(output_path + '\\' + ids + '\\' + bsids[i].replace(':', ',') + '\\' + 'hydrophobic interactions.csv', index=False)
                        """
                        print('hydrophobic interactions')
                        print(data_hi)
                        print(sep='/n')
                        """
                # hydrogen bonds ligand
                data_hb_l = pd.DataFrame()

                if len(s.hbonds_ldon) != 0:

                    for k in range(len(s.hbonds_ldon)):

                        data_hb_l_temp = []

                        data_hb_l_temp.append(s.hbonds_ldon[k][10])
                        data_hb_l_temp.append(s.hbonds_ldon[k][1])
                        data_hb_l_temp.append(s.hbonds_ldon[k][3])
                        data_hb_l_temp.append(s.hbonds_ldon[k][5])
                        data_hb_l_temp.append(s.hbonds_ldon[k][6])
                        data_hb_l_temp.append(s.hbonds_ldon[k][7])
                        data_hb_l_temp.append(s.hbonds_ldon[k][8])
                        data_hb_l_temp.append(s.hbonds_ldon[k][9])
                        data_hb_l_temp.append(s.hbonds_ldon[k][11])
                        data_hb_l_temp.append(s.hbonds_ldon[k][12])
                        data_hb_l_temp.append(s.hbonds_ldon[k][13])
                        data_hb_l_temp.append(s.hbonds_ldon[k][14])
                        data_hb_l_temp.append(s.hbonds_ldon[k][15])
                        data_hb_l_temp.append(s.hbonds_ldon[k][16])
                        data_hb_l_temp.append(s.hbonds_ldon[k][17])
                        data_hb_l_temp.append(s.hbonds_ldon[k][18])

                        data_hb_l_temp = pd.DataFrame(data_hb_l_temp).T
                        data_hb_l = data_hb_l.append(data_hb_l_temp)

                    data_hb_l.columns = ['resnr','a_orig_idx','d_orig_idx','distance_ah','distance_ad','angle','type','protisdon',
                                         'restype','reschain','resnr_l','restype_l','reschain_l','sidechain','atype','dtype']
                    data_hb_l = data_hb_l.sort_values(by='resnr',ascending=True)

                    if len(data_hb_l) != 0:
                        data_hb_l.to_csv(output_path + '\\' + ids + '\\' + bsids[i].replace(':', ',') + '\\' + 'hydrogen bonds ligand.csv', index=False)
                        """
                        print('hydrogen bonds ligand')
                        print(data_hb_l)
                        print(sep='/n')
                        """
                # hydrogen bonds protein
                data_hb_p = pd.DataFrame()

                if len(s.hbonds_pdon) != 0:

                    for k in range(len(s.hbonds_pdon)):

                        data_hb_p_temp = []
                        data_hb_p_temp.append(s.hbonds_pdon[k][10])
                        data_hb_p_temp.append(s.hbonds_pdon[k][1])
                        data_hb_p_temp.append(s.hbonds_pdon[k][3])
                        data_hb_p_temp.append(s.hbonds_pdon[k][5])
                        data_hb_p_temp.append(s.hbonds_pdon[k][6])
                        data_hb_p_temp.append(s.hbonds_pdon[k][7])
                        data_hb_p_temp.append(s.hbonds_pdon[k][8])
                        data_hb_p_temp.append(s.hbonds_pdon[k][9])
                        data_hb_p_temp.append(s.hbonds_pdon[k][11])
                        data_hb_p_temp.append(s.hbonds_pdon[k][12])
                        data_hb_p_temp.append(s.hbonds_pdon[k][13])
                        data_hb_p_temp.append(s.hbonds_pdon[k][14])
                        data_hb_p_temp.append(s.hbonds_pdon[k][15])
                        data_hb_p_temp.append(s.hbonds_pdon[k][16])
                        data_hb_p_temp.append(s.hbonds_pdon[k][17])
                        data_hb_p_temp.append(s.hbonds_pdon[k][18])

                        data_hb_p_temp = pd.DataFrame(data_hb_p_temp).T
                        data_hb_p = data_hb_p.append(data_hb_p_temp)

                    data_hb_p.columns = ['resnr','a_orig_idx','d_orig_idx','distance_ah','distance_ad','angle','type','protisdon',
                                         'restype','reschain','resnr_l','restype_l','reschain_l','sidechain','atype','dtype']
                    data_hb_p = data_hb_p.sort_values(by='resnr',ascending=True)

                    if len(data_hb_p) != 0:
                        data_hb_p.to_csv(output_path + '\\' + ids + '\\' + bsids[i].replace(':', ',') + '\\' + 'hydrogen bonds protein.csv', index=False)
                        """
                        print('hydrogen bonds protein')
                        print(data_hb_p)
                        print(sep='/n')
                        """
                # saltbridge ligand
                data_sb_l = pd.DataFrame()

                if len(s.saltbridge_lneg) != 0:

                    for k in range(len(s.saltbridge_lneg)):

                        data_sb_l_temp = []

                        data_sb_l_temp.append(s.saltbridge_lneg[k][4])
                        data_sb_l_temp.append(s.saltbridge_lneg[k][2])
                        data_sb_l_temp.append(s.saltbridge_lneg[k][3])
                        data_sb_l_temp.append(s.saltbridge_lneg[k][5])
                        data_sb_l_temp.append(s.saltbridge_lneg[k][6])
                        data_sb_l_temp.append(s.saltbridge_lneg[k][7])
                        data_sb_l_temp.append(s.saltbridge_lneg[k][8])
                        data_sb_l_temp.append(s.saltbridge_lneg[k][9])
                        data_sb_l_temp.append(s.saltbridge_lneg[k][0][1])
                        data_sb_l_temp.append(s.saltbridge_lneg[k][0][2])
                        data_sb_l_temp.append(s.saltbridge_lneg[k][0][3])
                        data_sb_l_temp.append(s.saltbridge_lneg[k][0][4])
                        data_sb_l_temp.append(s.saltbridge_lneg[k][0][5])
                        data_sb_l_temp.append(s.saltbridge_lneg[k][0][6])
                        data_sb_l_temp.append(s.saltbridge_lneg[k][1][2])
                        data_sb_l_temp.append(s.saltbridge_lneg[k][1][3])
                        data_sb_l_temp.append(s.saltbridge_lneg[k][1][4])
                        data_sb_l_temp.append(s.saltbridge_lneg[k][1][5])

                        data_sb_l_temp = pd.DataFrame(data_sb_l_temp).T
                        data_sb_l = data_sb_l.append(data_sb_l_temp)

                    data_sb_l.columns = ['resnr','distance','protispos','restype','reschain','resnr_l','restype_l','reschain_l','atoms_orig_idx','type',
                                         'center','restype','resnr1','reschain','atoms_orig_idx','type','center','fgroup']
                    data_sb_l = data_sb_l.sort_values(by='resnr',ascending=True)

                    if len(data_sb_l) != 0:
                        data_sb_l.to_csv(output_path + '\\' + ids + '\\' + bsids[i].replace(':', ',') + '\\' + 'saltbridge ligand.csv', index=False)
                        """
                        print('saltbridge ligand')
                        print(data_sb_l)
                        print(sep='/n')
                        """
                # saltbridge protein
                data_sb_p = pd.DataFrame()

                if len(s.saltbridge_pneg) != 0:

                    for k in range(len(s.saltbridge_pneg)):

                        data_sb_p_temp = []

                        data_sb_p_temp.append(s.saltbridge_pneg[k][4])
                        data_sb_p_temp.append(s.saltbridge_pneg[k][2])
                        data_sb_p_temp.append(s.saltbridge_pneg[k][3])
                        data_sb_p_temp.append(s.saltbridge_pneg[k][5])
                        data_sb_p_temp.append(s.saltbridge_pneg[k][6])
                        data_sb_p_temp.append(s.saltbridge_pneg[k][7])
                        data_sb_p_temp.append(s.saltbridge_pneg[k][8])
                        data_sb_p_temp.append(s.saltbridge_pneg[k][9])
                        data_sb_p_temp.append(s.saltbridge_pneg[k][0][2])
                        data_sb_p_temp.append(s.saltbridge_pneg[k][0][3])
                        data_sb_p_temp.append(s.saltbridge_pneg[k][0][4])
                        data_sb_p_temp.append(s.saltbridge_pneg[k][0][5])
                        data_sb_p_temp.append(s.saltbridge_pneg[k][1][1])
                        data_sb_p_temp.append(s.saltbridge_pneg[k][1][2])
                        data_sb_p_temp.append(s.saltbridge_pneg[k][1][3])
                        data_sb_p_temp.append(s.saltbridge_pneg[k][1][4])
                        data_sb_p_temp.append(s.saltbridge_pneg[k][1][5])
                        data_sb_p_temp.append(s.saltbridge_pneg[k][1][6])

                        data_sb_p_temp = pd.DataFrame(data_sb_p_temp).T
                        data_sb_p = data_sb_p.append(data_sb_p_temp)

                    data_sb_p.columns = ['resnr','distance','protispos','restype','reschain','resnr_l','restype_l','reschain_l','atoms_orig_idx',
                                         'type','center','fgroup','atoms_orig_idx','type','center','restype','resnr1','reschain']
                    data_sb_p = data_sb_p.sort_values(by='resnr',ascending=True)

                    if len(data_sb_p) != 0:
                        data_sb_p.to_csv(output_path + '\\' + ids + '\\' + bsids[i].replace(':', ',') + '\\' + 'saltbridge protein.csv', index=False)
                        """
                        print('saltbridge protein')
                        print(data_sb_p)
                        print(sep='/n')
                        """
                # metal complex
                data_metc = pd.DataFrame()

                if len(s.metal_complexes) != 0:

                    for k in range(len(s.metal_complexes)):

                        data_metc_temp = []

                        data_metc_temp.append(s.metal_complexes[k][8])
                        data_metc_temp.append(s.metal_complexes[k][1])
                        data_metc_temp.append(s.metal_complexes[k][2])
                        data_metc_temp.append(s.metal_complexes[k][3][1])
                        data_metc_temp.append(s.metal_complexes[k][3][2])
                        data_metc_temp.append(s.metal_complexes[k][3][3])
                        data_metc_temp.append(s.metal_complexes[k][3][4])
                        data_metc_temp.append(s.metal_complexes[k][3][5])
                        data_metc_temp.append(s.metal_complexes[k][3][6])
                        data_metc_temp.append(s.metal_complexes[k][4])
                        data_metc_temp.append(s.metal_complexes[k][5])
                        data_metc_temp.append(s.metal_complexes[k][6])
                        data_metc_temp.append(s.metal_complexes[k][7])
                        data_metc_temp.append(s.metal_complexes[k][9])
                        data_metc_temp.append(s.metal_complexes[k][10])
                        data_metc_temp.append(s.metal_complexes[k][11])
                        data_metc_temp.append(s.metal_complexes[k][12])
                        data_metc_temp.append(s.metal_complexes[k][13])
                        data_metc_temp.append(s.metal_complexes[k][14])
                        data_metc_temp.append(s.metal_complexes[k][15])
                        data_metc_temp.append(s.metal_complexes[k][16])
                        data_metc_temp.append(s.metal_complexes[k][17])
                        data_metc_temp.append(s.metal_complexes[k][18])    
                        if len(s.metal_complexes[k][3]) == 9:
                            data_metc_temp[3] = s.metal_complexes[k][3][2]
                            data_metc_temp[4] = s.metal_complexes[k][3][3]
                            data_metc_temp[5] = s.metal_complexes[k][3][5]
                            data_metc_temp[6] = s.metal_complexes[k][3][6]
                            data_metc_temp[7] = s.metal_complexes[k][3][7]
                            data_metc_temp[8] = s.metal_complexes[k][3][8]

                        data_metc_temp = pd.DataFrame(data_metc_temp).T
                        data_metc = data_metc.append(data_metc_temp)

                    data_metc.columns = ['resnr','metal_orig_idx','metal_type','atom_orig_idx','type','restype','resnr1','reschain','location',
                                         'target_orig_idx','target_type','coordination_num','distance','restype','reschain','restype_l','reschain_l',
                                         'resnr_l','location','rms','geometry','num_partners','complexnum']
                    data_metc = data_metc.sort_values(by='resnr',ascending=True)

                    if len(data_metc) != 0:
                        data_metc.to_csv(output_path + '\\' + ids + '\\' + bsids[i].replace(':', ',') + '\\' + 'metal complex.csv', index=False)
                        """
                        print('metal complex')
                        print(data_metc)
                        print(sep='/n')
                        """
                # Pication interactions ligand
                data_pic_l = pd.DataFrame()

                if len(s.pication_laro) != 0:

                    for k in range(len(s.pication_laro)):

                        data_pic_l_temp = []

                        data_pic_l_temp.append(s.pication_laro[k][6])
                        data_pic_l_temp.append(s.pication_laro[k][2])
                        data_pic_l_temp.append(s.pication_laro[k][3])
                        data_pic_l_temp.append(s.pication_laro[k][4])
                        data_pic_l_temp.append(s.pication_laro[k][5])
                        data_pic_l_temp.append(s.pication_laro[k][7])
                        data_pic_l_temp.append(s.pication_laro[k][8])
                        data_pic_l_temp.append(s.pication_laro[k][9])
                        data_pic_l_temp.append(s.pication_laro[k][10])
                        data_pic_l_temp.append(s.pication_laro[k][11])
                        data_pic_l_temp.append(s.pication_laro[k][1][2])
                        data_pic_l_temp.append(s.pication_laro[k][1][3])
                        data_pic_l_temp.append(s.pication_laro[k][1][4])
                        data_pic_l_temp.append(s.pication_laro[k][1][5])

                        data_pic_l_temp = pd.DataFrame(data_pic_l_temp).T
                        data_pic_l = data_pic_l.append(data_pic_l_temp)

                    data_pic_l.columns = ['resnr','distance','offset','type','restype','reschain','restype_l','resnr_l','reschain_l','protcharged',
                                          'atoms_orig_idx','type','center','fgroup']
                    data_pic_l = data_pic_l.sort_values(by='resnr',ascending=True)

                    if len(data_pic_l) != 0:
                        data_pic_l.to_csv(output_path + '\\' + ids + '\\' + bsids[i].replace(':', ',') + '\\' + 'Pication interactions ligand.csv', index=False)
                        """
                        print('Pication interactions ligand')
                        print(data_pic_l)
                        print(sep='/n')
                        """
                # Pication interactions protein
                data_pic_p = pd.DataFrame()

                if len(s.pication_paro) != 0:

                    for k in range(len(s.pication_paro)):

                        data_pic_p_temp = []

                        data_pic_p_temp.append(s.pication_paro[k][6])
                        data_pic_p_temp.append(s.pication_paro[k][2])
                        data_pic_p_temp.append(s.pication_paro[k][3])
                        data_pic_p_temp.append(s.pication_paro[k][4])
                        data_pic_p_temp.append(s.pication_paro[k][5])
                        data_pic_p_temp.append(s.pication_paro[k][7])
                        data_pic_p_temp.append(s.pication_paro[k][8])
                        data_pic_p_temp.append(s.pication_paro[k][9])
                        data_pic_p_temp.append(s.pication_paro[k][10])
                        data_pic_p_temp.append(s.pication_paro[k][11])
                        data_pic_p_temp.append(s.pication_paro[k][1][2])
                        data_pic_p_temp.append(s.pication_paro[k][1][3])
                        data_pic_p_temp.append(s.pication_paro[k][1][4])
                        data_pic_p_temp.append(s.pication_paro[k][1][5])

                        data_pic_p_temp = pd.DataFrame(data_pic_p_temp).T
                        data_pic_p = data_pic_p.append(data_pic_p_temp)

                    data_pic_p.columns = ['resnr','distance','offset','type','restype','reschain','restype_l','resnr_l','reschain_l','protcharged',
                                          'atoms_orig_idx','type','center','fgroup']
                    data_pic_p = data_pic_p.sort_values(by='resnr',ascending=True)

                    if len(data_pic_p) != 0:
                        data_pic_p.to_csv(output_path + '\\' + ids + '\\' + bsids[i].replace(':', ',') + '\\' + 'Pication interactions protein.csv', index=False)
                        """
                        print('Pication interactions proteind')
                        print(data_pic_p)
                        print(sep='/n')
                        """
                # waterbridges
                data_wb = pd.DataFrame()

                if len(s.water_bridges) != 0:

                    for k in range(len(s.water_bridges)):

                        data_wb_temp = []
                        data_wb_temp.append(s.water_bridges[k][14])
                        data_wb_temp.append(s.water_bridges[k][1])
                        data_wb_temp.append(s.water_bridges[k][2])
                        data_wb_temp.append(s.water_bridges[k][4])
                        data_wb_temp.append(s.water_bridges[k][5])
                        data_wb_temp.append(s.water_bridges[k][8])
                        data_wb_temp.append(s.water_bridges[k][9])
                        data_wb_temp.append(s.water_bridges[k][10])
                        data_wb_temp.append(s.water_bridges[k][11])
                        data_wb_temp.append(s.water_bridges[k][12])
                        data_wb_temp.append(s.water_bridges[k][13])
                        data_wb_temp.append(s.water_bridges[k][15])
                        data_wb_temp.append(s.water_bridges[k][16])
                        data_wb_temp.append(s.water_bridges[k][17])
                        data_wb_temp.append(s.water_bridges[k][18])
                        data_wb_temp.append(s.water_bridges[k][19])
                        data_wb_temp.append(s.water_bridges[k][20])

                        data_wb_temp = pd.DataFrame(data_wb_temp).T
                        data_wb = data_wb.append(data_wb_temp)

                    data_wb.columns = ['resnr','a_orig_idx','atype','AAd_orig_idx','dtype','water_orig_idx','distance_aw',
                                       'distance_dw','d_angle','w_angle','type',
                                       'restype','reschain','resnr_l','restype_l','reschain_l','protisdon']
                    data_wb = data_wb.sort_values(by='resnr',ascending=True)

                    if len(data_wb) != 0:
                        data_wb.to_csv(output_path + '\\' + ids + '\\' + bsids[i].replace(':', ',') + '\\' + 'waterbridges.csv', index=False)
                        """
                        print('waterbridgesd')
                        print(data_wb)
                        print(sep='/n')
                        """
                # halogen bonds
                data_halb = pd.DataFrame()

                if len(s.halogen_bonds) != 0:

                    for k in range(len(s.halogen_bonds)):

                        data_halb_temp = []

                        data_halb_temp.append(s.halogen_bonds[k][8])
                        data_halb_temp.append(s.halogen_bonds[k][4])
                        data_halb_temp.append(s.halogen_bonds[k][5])
                        data_halb_temp.append(s.halogen_bonds[k][6])
                        data_halb_temp.append(s.halogen_bonds[k][7])
                        data_halb_temp.append(s.halogen_bonds[k][9])
                        data_halb_temp.append(s.halogen_bonds[k][10])
                        data_halb_temp.append(s.halogen_bonds[k][11])
                        data_halb_temp.append(s.halogen_bonds[k][12])
                        data_halb_temp.append(s.halogen_bonds[k][13])
                        data_halb_temp.append(s.halogen_bonds[k][14])
                        data_halb_temp.append(s.halogen_bonds[k][15])
                        data_halb_temp.append(s.halogen_bonds[k][0][1])
                        data_halb_temp.append(s.halogen_bonds[k][0][3])
                        data_halb_temp.append(s.halogen_bonds[k][1])
                        data_halb_temp.append(s.halogen_bonds[k][2][2])
                        data_halb_temp.append(s.halogen_bonds[k][2][4])
                        data_halb_temp.append(s.halogen_bonds[k][3])

                        data_halb_temp = pd.DataFrame(data_halb_temp).T
                        data_halb = data_halb.append(data_halb_temp)

                    data_halb.columns = ['resnr','distance','don_angle','acc_angle','restype','reschain','restype_l','resnr_l',
                                         'reschain_l','donortype','acctype','sidechain','o_orig_idx','y_orig_idx','acc_orig_idx','x_orig_idx',
                                         'c_orig_idx','don_orig_idx']
                    data_halb = data_halb.sort_values(by='resnr',ascending=True)

                    if len(data_halb) != 0:
                        data_halb.to_csv(output_path + '\\' + ids + '\\' + bsids[i].replace(':', ',') + '\\' + 'halogen bonds.csv', index=False)
                        """
                        print('halogen bonds')
                        print(data_halb)
                        print(sep='/n')
                        """
                # pistack
                data_pi = pd.DataFrame()

                if len(s.pistacking) != 0:

                    for k in range(len(s.pistacking)):


                        data_pi_temp = []
                        data_pi_temp.append(s.pistacking[k][7])
                        data_pi_temp.append(s.pistacking[k][2])
                        data_pi_temp.append(s.pistacking[k][3])
                        data_pi_temp.append(s.pistacking[k][4])
                        data_pi_temp.append(s.pistacking[k][5])
                        data_pi_temp.append(s.pistacking[k][6])
                        data_pi_temp.append(s.pistacking[k][8])
                        data_pi_temp.append(s.pistacking[k][9])
                        data_pi_temp.append(s.pistacking[k][10])
                        data_pi_temp.append(s.pistacking[k][11])
                        data_pi_temp.append(s.pistacking[k][1][2])
                        data_pi_temp.append(s.pistacking[k][1][3])
                        data_pi_temp.append(s.pistacking[k][1][5])
                        data_pi_temp.append(s.pistacking[k][1][6])

                        data_pi_temp = pd.DataFrame(data_pi_temp).T
                        data_pi = data_pi.append(data_pi_temp)

                    data_pi.columns = ['resnr','distance','angle','offset','type','restype',
                                       'reschain','restype_l','resnr_l','reschain_l','atoms_orig_idx','normal','center','type']
                    data_pi = data_pi.sort_values(by='resnr',ascending=True)

                    if len(data_hi) != 0:
                        data_pi.to_csv(output_path + '\\' + ids + '\\' + bsids[i].replace(':', ',') + '\\' + 'pistack.csv', index=False)
                        """
                        print('pistack')
                        print(data_pi)
                        print(sep='/n')
                        """
        src = 'D:\\task\\adamNet2.0\\code_an\\PLIP\\Plibpdb\\PDB2\\'
        dir = 'D:\\task\\adamNet2.0\\code_an\\PLIP\\Plibpdb\\PDB2\\finished\\'
        shutil.move(src + ids + '.pdb', dir + ids + '.pdb')

        gc.collect()
    except:
        src = 'D:\\task\\adamNet2.0\\code_an\\PLIP\\Plibpdb\\PDB2\\'
        dir = 'D:\\task\\adamNet2.0\\code_an\\PLIP\\Plibpdb\\PDB2\\pdberror\\'
        shutil.move(src + ids + '.pdb', dir + ids + '.pdb')