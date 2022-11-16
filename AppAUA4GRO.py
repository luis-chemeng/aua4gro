#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 18 11:13:37 2022

Script to generate simulation files for Gromacs

@author: luis
"""

import numpy as np
from itertools import product as cart_prod
from openbabel import pybel
import AUA_parameters_2 as params
from utilities import AUA4_mol
from utilities import find_height_pyramid
from utilities import compare_list
from utilities import root_path
import os
import shutil


try:
    """
    Receive and validate a input SIMLES code
    """
    
    smiles = input('Type SMILES code of a valid molecule: ')
    # smiles = 'C(CO)N'  # MEA
    
    # Check supported funtional groups
    smi_chars = set(smiles)
    issupp = len(smi_chars - params.supp_smi_chars)==0
    
    if not issupp:
        print('The SMILES code read contain not supported functional groups.\n')
        print('It may cause erros.')
        ans = input('Would you like to continue? [YES/no]')
        if ans.lower()=='no':
            raise InterruptedError('Thanks.')
    
    isopt = False
    if smiles in params.opt_tops.keys():
        print('For this alkanolamine exist an optimized set of charges.')
        print('Aditional files with optimized set of charges will be genrated.')
        isopt = True
    
    mol = pybel.readstring('smi', smiles.upper())
    mol.addh()
    mol.make3D()
    
    dict_aua_mol = dict()
    for atom in mol.atoms:
        dict_aua_mol.update({atom.idx : []})
    
    """
    Identify Functional Groups
    """
    
    dict_chgs_found = dict()
    for code, chg_params in params.dict_chgs.items():
        smarts = pybel.Smarts(chg_params[0])
        chgs_found = smarts.findall(mol)
        if len(chgs_found) > 0:
            dict_chgs_found.update({chg_params[1] : chgs_found})
            print('It was identified {} atoms of type {}'.format(len(chgs_found), chg_params[1]))
        for chg_found in chgs_found:
            dict_aua_mol.update({chg_found[0] : [code] + chg_params})
       
    dict_fgs_found = dict()
    for fg, smart_code in params.dict_fgs.items():
        smarts = pybel.Smarts(smart_code)
        fgs_found = smarts.findall(mol)
        if len(fgs_found) > 0:
            dict_fgs_found.update({fg : fgs_found})
            print('It was identified {} functional group(s) of type {}'.format(len(fgs_found), fg))
    
    # If there is any atom of type aua0, raise an alert and finish 
    
    # Count the number of functional groups
    idxs_fgs = []
    for fg, list_idxs in dict_fgs_found.items():
        for i in range(len(list_idxs)):
            idxs_fgs.append(list_idxs[i])
    
    if len(idxs_fgs)>1:
        ismultifunc = True
    else:
        ismultifunc = False
    
    """
    Assign electric charge to alpha carbons.
    """
    
    auamol = AUA4_mol(dict_aua_mol, mol)
    
    for fg, list_idxs in dict_fgs_found.items():
        if fg == 'Alcohol':
            for idxs in list_idxs:
                idxoi = idxs[0] 
                auamol.ModCharge(idxoi, 0.265)
        elif fg == 'Ether':
            for idxs in list_idxs:
                idxoi1 = idxs[0]
                idxoi2 = idxs[1]
                idxoi3 = idxs[2]
                if smiles=='COC' or ismultifunc:
                    auamol.ModCharge(idxoi1, 0.223)
                    auamol.ModCharge(idxoi2, -0.446)
                    auamol.ModCharge(idxoi3, 0.223)
                else:
                    auamol.ModCharge(idxoi1, 0.185)
                    auamol.ModCharge(idxoi2, -0.37)
                    auamol.ModCharge(idxoi3, 0.185)
        elif fg == 'Aldehyde':
            for idxs in list_idxs:
                idxoi = idxs[0]
                auamol.ModCharge(idxoi, 0.46)
        elif fg == 'Ketone':
            for idxs in list_idxs:
                idxoi = idxs[0]
                auamol.ModCharge(idxoi, 0.49)
        elif fg == 'CarboxAcid':
            for idxs in list_idxs:
                idxoi1 = idxs[0]
                idxoi2 = idxs[1]
                auamol.ModCharge(idxoi1, -0.120)
                auamol.ModCharge(idxoi2, 0.658)
        elif fg == 'Ester':
            for idxs in list_idxs:
                idxoi1 = idxs[0]
                idxoi2 = idxs[1]
                idxoi3 = idxs[-1]
                auamol.ModCharge(idxoi1, -0.089)
                auamol.ModCharge(idxoi2, 0.484)
                auamol.ModCharge(idxoi3, 0.319)
        elif fg == 'Amine1':
            for idxs in list_idxs:
                idxoi = idxs[0]
                auamol.ModCharge(idxoi, 0.180)
        elif fg == 'Amine2':
            for idxs in list_idxs:
                idxoi1 = idxs[0]
                idxoi2 = idxs[-1]
                auamol.ModCharge(idxoi1, 0.176)
                auamol.ModCharge(idxoi2, 0.176)
        elif fg == 'Amine3':
            for idxs in list_idxs:
                idxoi1 = idxs[0]
                idxoi2 = idxs[-1]
                idxoi3 = idxs[-2]
                auamol.ModCharge(idxoi1, 0.230)
                auamol.ModCharge(idxoi2, 0.230)
                auamol.ModCharge(idxoi3, 0.230)
                
    
    """
    Assign Bonded interactions
    """
    
    # Bonds
    for smart_code, b_params in params.bondtypes.items():
        smarts = pybel.Smarts(smart_code)
        b_found = smarts.findall(mol)
        for idx in b_found:
            auamol.AddBond(list(idx)+b_params)
        
    # Angles
    for smart_code, a_params in params.angletypes.items():
        smarts = pybel.Smarts(smart_code)
        a_found = smarts.findall(mol)
        for idx in a_found:
            auamol.AddAngle(list(idx)+a_params)
    
    # Dihedral angles'CN(CCO)CCO'
    for smart_code, d_params in params.dihtypes.items():
        smarts = pybel.Smarts(smart_code)
        d_found = smarts.findall(mol)
        for idx in d_found:
            auamol.AddDihedral(list(idx)+d_params)
    
    # Verify if all the parameters were assigned
    
    smarts = pybel.Smarts('[*;!$([#1]-[#6])]~[*;!$([#1]-[#6])]')
    b_found = smarts.findall(mol)
    int12 = [(idx[0], idx[-1]) for idx in b_found]
    smarts = pybel.Smarts('[*;!$([#1]-[#6])]~[*]~[*;!$([#1]-[#6])]')
    a_found = smarts.findall(mol)
    int13 = [(idx[0], idx[-1]) for idx in a_found]
    smarts = pybel.Smarts('[*;!$([#1]-[#6])]~[*]~[*]~[*;!$([#1]-[#6])]')
    d_found = smarts.findall(mol)
    int14 = [(idx[0], idx[-1]) for idx in d_found]
    int14_for_exc = [(idx[0], idx[-1]) for idx in d_found]
    
    # Interaction 1-4 necessary for MEA, ethylenediamne and ethylene glycol.
    smarts = pybel.Smarts('[N,O][*][*][N,O]')
    d_found = smarts.findall(mol)
    int14_extra = [(idx[0], idx[-1]) for idx in d_found]
    
    for int14i in int14_extra:    
        int14_for_exc.remove(int14i)
    
    ints = set(int12 + int13 + int14)
    
    auamol_bonds = [tuple(bond[:2]) for bond in auamol.bonds]
    auamol_angles = [tuple(angle[0:3]) for angle in auamol.angles]
    auamol_dihedrals = [tuple(dihedral[0:4]) for dihedral in auamol.dihedrals]
    
    int12_na, int12_wa = compare_list(b_found, auamol_bonds)
    int13_na, int13_wa = compare_list(a_found, auamol_angles)
    int14_na, int14_wa = compare_list(d_found, auamol_dihedrals)
    
    
    # na stands for non-assignated and wa stands for wrong-assignated
    if len(int12_na) != 0:
        print('There are not parameters for the following pairs of bonded atoms:')
        for i in list(int12_na):
            print('Bond {} - {}'.format(i[0], i[1]))
            
    if len(int13_na) != 0:
        print('There are not parameters for the following triads of bonded atoms:')
        for i in list(int13_na):
            print('Angle {} - {} - {}'.format(i[0], i[1], i[2]))
            
    if len(int14_na) != 0:
        print('There are not parameters for the following quads of bonded atoms:')
        for i in list(int14_na):
            print('Dihedral {} - {} - {} - {}'.format(i[0], i[1], i[2], i[3]))
    
    if len(int12_wa) != 0 or len(int12_wa) != 0 or len(int12_wa) != 0:
        raise ValueError('There was a wrong assignation of bonded interactions.')
    
    
    """
    Find a vector to set molecule in the box
    """
    xs = []
    ys = []
    zs = []
    
    for atom in auamol.atoms:
        x, y, z = atom.patom.coords
        xs.append(x)
        ys.append(y)
        zs.append(z)
        
    xs = np.array(xs)/10
    ys = np.array(ys)/10
    zs = np.array(zs)/10
    
    dx = xs.max()-xs.min()
    dy = ys.max()-ys.min()
    dz = zs.max()-zs.min()
    
    dimax = round(np.array([dx, dy, dz]).max()+0.05, 2)
    
    x_ = dimax/2-xs.mean()
    y_ = dimax/2-ys.mean()
    z_ = dimax/2-zs.mean()
    
    box = np.array([dimax, dimax, dimax])
    vec = np.array([x_, y_, z_])
    
    """
    Create and Assign Virtual sites
    """
    
    auamol.GenNewIdxsDict()        
    auamol.GenNumberVS()
    
    for atom in auamol.atoms:
        if atom.hasvs:
            nh_atoms = atom.FindNA(mol)
            vs_idx = auamol.atomidx_vs[atom.idx]
            vsc = params.vsclasf[atom.code]
            if vsc == '2a':
                b_d = auamol.FindBondDistance(atom.idx, nh_atoms[0])
                a = round(-1*atom.delta/b_d, 4)
                vs_par = [vs_idx, atom.idx, nh_atoms[0], 1, a]
                auamol.AddVS(vs_par, '2')
                
                p0 = np.array(atom.patom.coords)/10 + vec
                p1 = np.array(auamol.atoms[nh_atoms[0]-1].patom.coords)/10 + vec
                v1 = p1-p0
                pvs = a*v1+p0
                
            elif vsc == '2b':
                for nhi in nh_atoms:
                    if auamol.atoms[nhi-1].chg == 'Alcohol H' or auamol.atoms[nhi-1].chg == 'CarboxAcid H' or auamol.atoms[nhi-1].chg == 'Amine2 H':
                        nhi_sel = nhi
                b_d = auamol.FindBondDistance(atom.idx, nhi_sel)
                a = round(atom.delta/b_d, 4)
                vs_par = [vs_idx, atom.idx, nhi_sel, 1, a]
                auamol.AddVS(vs_par, '2')
                
                p0 = np.array(atom.patom.coords)/10 + vec
                p1 = np.array(auamol.atoms[nhi_sel-1].patom.coords)/10 + vec
                v1 = p1-p0
                pvs = a*v1+p0
                
            elif vsc == '3a':
                a = round(-1*atom.delta, 4)
                vs_par = [vs_idx, atom.idx, *nh_atoms, 2, 0.5, a]
                auamol.AddVS(vs_par, '3')
                
                p0 = np.array(atom.patom.coords)/10 + vec
                p1 = np.array(auamol.atoms[nh_atoms[0]-1].patom.coords)/10 + vec
                p2 = np.array(auamol.atoms[nh_atoms[1]-1].patom.coords)/10 + vec
                v1 = (p1-p0)/np.linalg.norm(p1-p0)
                v2 = (p2-p0)/np.linalg.norm(p2-p0)
                v3 = (v1+v2)/np.linalg.norm(v1+v2)
                pvs = a*v3+p0
                
            elif vsc == '3b':
                nhi_sel = []
                for nhi in nh_atoms:
                    if auamol.atoms[nhi-1].chg == 'Amine1 H':
                        nhi_sel.append(nhi)
                a = round(atom.delta, 4)
                vs_par = [vs_idx, atom.idx, *nhi_sel, 2, 0.5, a]
                auamol.AddVS(vs_par, '3')
                
                p0 = np.array(atom.patom.coords)/10 + vec
                p1 = np.array(auamol.atoms[nhi_sel[0]-1].patom.coords)/10 + vec
                p2 = np.array(auamol.atoms[nhi_sel[1]-1].patom.coords)/10 + vec
                v1 = (p1-p0)/np.linalg.norm(p1-p0)
                v2 = (p2-p0)/np.linalg.norm(p2-p0)
                v3 = (v1+v2)/np.linalg.norm(v1+v2)
                pvs = a*v3+p0
                
            elif vsc == 'n&2':
                raise TypeError('This molecule include a CH unitted atom. The construction of such virtual site currently is not supported.')
                vs_aux_idx = auamol.GenAuxiliarVS(atom.idx)
                vs_par = [vs_aux_idx, 1, *nh_atoms]
                auamol.AddVS(vs_par, 'n')
                b_ds = [auamol.FindBondDistance(atom.idx, atomi) for atomi in nh_atoms]
                a_ms = []
                a_ms.append(auamol.FindAngleMag(nh_atoms[0], atom.idx, nh_atoms[1]))
                a_ms.append(auamol.FindAngleMag(nh_atoms[1], atom.idx, nh_atoms[2]))
                a_ms.append(auamol.FindAngleMag(nh_atoms[0], atom.idx, nh_atoms[2]))
                h = find_height_pyramid(*b_ds, a_ms[0], a_ms[1], a_ms[2])
                a = round(-1*atom.delta/h, 4)
                vs_par = [vs_idx, atom.idx, vs_aux_idx, 1, a]
                auamol.AddVS(vs_par, '2_')
                
                p0 = np.array(atom.patom.coords)/10 + vec
                p1 = np.array(auamol.atoms[nh_atoms[0]-1].patom.coords)/10 + vec
                p2 = np.array(auamol.atoms[nh_atoms[1]-1].patom.coords)/10 + vec
                p3 = np.array(auamol.atoms[nh_atoms[2]-1].patom.coords)/10 + vec
                pvs1 = (p1+p2+p3)/3
                v1 = pvs1-p0
                pvs = a*v1+p0
                auamol.vs_pos.update({vs_aux_idx : pvs1})
            
            auamol.vs_pos.update({vs_idx : pvs})
    
    """
    Treat non-bonded intramolecular interactions: exclusions, pairs.
    """
    
    auamol.GenExclusions(int12, int13, int14_for_exc)
    
    extra_pairs = []
    for idxs in idxs_fgs:
        idxs_fgs_ = idxs_fgs.copy()
        idxs_fgs_.remove(idxs)
        for co_idxs in idxs_fgs_:
            pairs = [tuple(sorted(pair)) for pair in list(cart_prod(idxs, co_idxs))]
            extra_pairs = extra_pairs + pairs
    
    test1 = extra_pairs
            
    extra_pairs = set(extra_pairs).intersection(ints)
    
    """
    Write topology file
    """
    
    if isopt:
        molname = params.opt_tops[smiles]
        auamol.short_name = molname
        src = root_path+'/opt/'+molname+'_opt'+'.gro'
        dst = os.getcwd()+'/'+molname+'_opt'+'.gro'
        shutil.copy2(src, dst)
        src = root_path+'/opt/'+molname+'_opt'+'.top'
        dst = os.getcwd()+'/'+molname+'_opt'+'.top'
        shutil.copy2(src, dst)
        
    
    
    f = open(auamol.short_name+'.top', 'w')
    
    f.write(auamol.WriteDefaults(ismultifunc))
    f.write(auamol.WriteSec1())
    f.write(auamol.WriteMoleculetype())
    f.write(auamol.WriteSec2())
    f.write(auamol.WriteSec3(int_na=list(int12_na)))
    f.write(auamol.WriteSec4(int_na=list(int13_na)))
    f.write(auamol.WriteSec5(int_na=list(int14_na)))
    f.write(auamol.WriteVS())
    f.write(auamol.WriteExclusions())
    f.write(auamol.WritePairs(extra_pairs))
    f.write(auamol.WriteSystem())
    
    f.close()
    
    
    """
    Write coordinate file
    """
    
    f = open(auamol.short_name+'.gro', 'w')
    f.write(auamol.WriteGro(vec, box))
    f.close()

except Exception as ex:
    print(ex)

    


# smarts = pybel.Smarts('[CX4H2]')
# test1 = smarts.findall(mol)
# test1

# for a in auamol.atoms:
#     print(a)
