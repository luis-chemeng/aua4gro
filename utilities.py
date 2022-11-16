#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 18 21:41:46 2022

@author: luis
"""

import numpy as np
from openbabel import pybel
import copy

import os

root_path = os.path.dirname(os.path.abspath(__file__))

class AUA4_atom(object):
    def __init__(self, idx, pybelatom, aua_code='aua0', smart_code='', chg='', mass=0, charge=0, sigma=0, epsilon=0, delta=0):
        # Mass, Charge, Sigma [nm], Epsilon [kJ/mol], delta [nm] 
        self.idx = idx
        self.patom = pybelatom
        self.code = aua_code
        self.smart_code = smart_code
        self.chg = chg
        if chg != '':
            self.ch_ele = chg.split(' ')[-1]
            self.fg = chg.split(' ')[0]
        else:
            self.ch_ele = '_'
            self.fg = '_'
        self.mass = mass
        self.charge = charge
        self.sig = sigma
        self.eps = epsilon
        self.delta = delta
        if delta != 0:
            self.hasvs = True
        else:
            self.hasvs = False
            
        self.name_gro = ''
        self.code_unique = ''

    def __str__(self):
        descriptor = 'Atom {} of {} in {} ({}).\n'.format(self.idx, self.ch_ele, self.fg, self.code)
        return descriptor
    
    def SetNameGro(self, auaidx):
        if auaidx < 10:
            self.name_gro = self.ch_ele[0] + '0' + str(auaidx)
            self.code_unique = self.code + '_0' + str(auaidx)
        else:
            self.name_gro = self.ch_ele[0] + str(auaidx)
            self.code_unique = self.code + '_' + str(auaidx)
    
    def FindNA(self, mol):
        couples = pybel.Smarts(self.smart_code+'~[*;!$([#1]-[#6])]').findall(mol)
        n_atoms = []
        for couple in couples:
            if couple[0] == self.idx:
                n_atoms.append(couple[1])
            elif couple[1] == self.idx:
                n_atoms.append(couple[0])
        
        return n_atoms
            
class AUA4_mol(object):
    def __init__(self, dict_pre_atoms, pybelmol):
        list_atoms = []
        for idx, props in dict_pre_atoms.items():
            list_atoms.append(AUA4_atom(idx, pybelmol.atoms[idx-1], *props))
        
        self.atoms = list_atoms
        self.pmol = pybelmol
        self.bonds = []
        self.angles = []
        self.dihedrals = []
        self.aua_idxs = dict()
        self.vs_defs = dict()
        self.vs_props = dict()
        self.atomidx_vs = dict()
        self.aux_vs = dict()
        self.vs2 = []
        self.vs2_ = []
        self.vs3 = []
        self.vsn = []
        self.exc = dict()
        self.vs_pos = dict()
        
        self.short_name = 'UNK'
    
    def __str__(self):
        
        descriptor = ''
        for atom in self.atoms:
            descriptor += atom.__str__()
        return descriptor
    
    def ModCharge(self, idx, new_value):
        for atom in self.atoms:
            if atom.idx == idx:
                if atom.charge == 0:
                    atom.charge = new_value
                else:
                    raise ValueError('The atom already has charge.')
                break
            
    def GenNewIdxsDict(self):
        aua_idx = 0
        for atom in self.atoms:
            if atom.code != 'aua0':
                aua_idx += 1
                self.aua_idxs.update({atom.idx : aua_idx})                
    
    def AddBond(self, bond_par):
        assert len(bond_par) == 5
        self.bonds.append(bond_par)
        
    def AddAngle(self, angle_par):
        assert len(angle_par) == 7
        self.angles.append(angle_par)
        
    def AddDihedral(self, dihedral_par):
        assert len(dihedral_par) == 12
        self.dihedrals.append(dihedral_par)
    
    def AddVS(self, vs_par, vstype):
        if vstype == '2':
            self.vs2.append(vs_par)
        elif vstype == '2_':
            self.vs2_.append(vs_par)
        elif vstype == '3':
            self.vs3.append(vs_par)
        elif vstype == 'n':
            self.vsn.append(vs_par)
    
    def FindBondDistance(self, atom1, atom2):
        l1 = [atom1, atom2]
        l2 = [atom2, atom1]
        for bond_par in self.bonds:
            if bond_par[0:2] == l1 or bond_par[0:2] == l2:
                value_found = bond_par[3]
                break
            value_found = None
        
        return value_found
    
    def FindAngleMag(self, atom1, atom2, atom3):
        a1 = [atom1, atom2, atom3]
        a2 = [atom3, atom2, atom1]
        for angle_par in self.angles:
            if angle_par[0:3] == a1 or angle_par[0:3] == a2:
                value_found = angle_par[4]
                break
            value_found = None
        
        return value_found
        
    def GenNumberVS(self):
        aua_idx = max(self.aua_idxs.values())
        for atom in self.atoms:
            if atom.hasvs:
                aua_idx += 1
                self.atomidx_vs.update({atom.idx : aua_idx})
                vs_name = 'V'+str(self.aua_idxs[atom.idx])
                params = [str(aua_idx), vs_name, '1', self.short_name, vs_name, '1', '0.000', '0.000']
                self.vs_props.update({aua_idx : params})
                params = [vs_name, vs_name, '0.000', '0.000', 'D', stround(atom.sig, 4), stround(atom.eps, 4)]
                self.vs_defs.update({aua_idx : params})
                
    def GenAuxiliarVS(self, atom_idx):
        aua_idx = max(self.aua_idxs.values())+len(self.atomidx_vs.values())+len(self.aux_vs.items())
        aua_idx += 1
        self.aux_vs.update({atom_idx : aua_idx})
        vs_name = 'V'+str(self.aua_idxs[atom_idx])+'a'
        params = [str(aua_idx), vs_name, '1', self.short_name, vs_name, '1', '0.000', '0.000']
        self.vs_props.update({aua_idx : params})
        params = [vs_name, vs_name, '0.000', '0.000', 'D', '0.000', '0.000']
        self.vs_defs.update({aua_idx : params})
        return aua_idx
    
    def WriteDefaults(self, ismultifunc):
        cad = '\n[  defaults  ]\n'
        cad += fill(';nbfunc', 10) + fill('comb-rule', 10) + fill('gen-pairs', 12)	+ fill('fudgeLJ', 10) + fill('fudgQQ', 10) + '\n'
        if ismultifunc:
            cad += fill('1', 10) + fill('2', 10) + fill('yes', 12) + fill('0', 10) + fill('1', 10) + '\n'
        else:
            cad += fill('1', 10) + fill('2', 10) + fill('no', 12) + '\n'
        
        return cad
    
    def WriteSec1(self):
        cad = '\n[  atomtypes  ]\n'
        for atom in self.atoms:
            if atom.code != 'aua0':
                auaidx = self.aua_idxs[atom.idx]
                atom.SetNameGro(auaidx)
                cad += fill(atom.code_unique, 10) + fill(atom.name_gro, 8)
                cad += fill(str(atom.mass), 8) + fill(str(atom.charge), 8)
                if atom.hasvs:
                    cad += fill('A', 8) + fill('0.000', 8) + fill('0.000', 8) + '\n'
                else:
                    cad += fill('A', 8) + fill(stround(atom.sig, 4), 8) + fill(stround(atom.eps, 4), 8) + '\n'
                
        for vsidx, vsterms in self.vs_defs.items():
            cad += fill(vsterms[0], 10)
            for vsterm in vsterms[1:]:
                cad += fill(vsterm, 8)
            cad += '\n'
                
        return cad
    
    def WriteMoleculetype(self):
        cad = '\n[  moleculetype  ]\n'
        cad += fill(';name', 8) + fill('nrexcl', 8) + '\n'
        cad += fill(self.short_name, 8) + fill('3', 8) + '\n'
        return cad
    
    def WriteSec2(self):
        cad = '\n[  atoms  ]\n'
        for atom in self.atoms:
            if atom.code != 'aua0':
                auaidx = self.aua_idxs[atom.idx]
                cad += fill(str(auaidx), 8) + fill(atom.code_unique, 10)
                cad += fill('1', 8) + fill(self.short_name, 8)
                cad += fill(atom.name_gro, 8) + fill('1', 8)
                cad += fill(stround(atom.charge, 4), 8) + fill(stround(atom.mass, 4), 8) + '\n'
        
        for vsidx, vsterms in self.vs_defs.items():
            cad += fill(str(vsidx), 8) + fill(vsterms[0], 10)
            cad += fill('1', 8) + fill(self.short_name, 8)
            cad += fill(vsterms[0], 8) + fill('1', 8)
            cad += fill(vsterms[2], 8) + fill(vsterms[3], 8)
            cad += '\n'
                
        return cad
    
    def WriteSec3(self, int_na=[]):
        cad = '\n[  constraints  ]\n'
        for bond in self.bonds:
            cad += fill(str(self.aua_idxs[bond[0]]), 4)
            cad += fill(str(self.aua_idxs[bond[1]]), 4)
            cad += fill(str(bond[2]), 4)
            cad += fill(stround(bond[3], 4), 9)
            cad += fill('    ;'+bond[4], 30, rev=True) + '\n'
        
        for bond in int_na:
            cad += ';'
            cad += fill(str(self.aua_idxs[bond[0]]), 4)
            cad += fill(str(self.aua_idxs[bond[1]]), 4)
            cad += '; Missing bond\n'
            
        return cad    
        
    def WriteSec4(self, int_na=[]):
        cad = '\n[  angles  ]\n'
        for angle in self.angles:
            cad += fill(str(self.aua_idxs[angle[0]]), 4)
            cad += fill(str(self.aua_idxs[angle[1]]), 4)
            cad += fill(str(self.aua_idxs[angle[2]]), 4)
            cad += fill(str(angle[3]), 4)
            cad += fill(stround(angle[4], 2), 9)
            cad += fill(stround(angle[5], 3), 9)
            cad += fill('    ;'+angle[6], 35, rev=True) + '\n'
        
        for angle in int_na:
            cad += ';'
            cad += fill(str(self.aua_idxs[angle[0]]), 4)
            cad += fill(str(self.aua_idxs[angle[1]]), 4)
            cad += fill(str(self.aua_idxs[angle[2]]), 4)
            cad += '; Missing angle\n'
        
        return cad
    
    def WriteSec5(self, int_na=[]):
        cad = '\n[  dihedrals  ]\n'
        for dih in self.dihedrals:
            cad += fill(str(self.aua_idxs[dih[0]]), 4)
            cad += fill(str(self.aua_idxs[dih[1]]), 4)
            cad += fill(str(self.aua_idxs[dih[2]]), 4)
            cad += fill(str(self.aua_idxs[dih[3]]), 4)
            cad += fill(str(dih[4]), 4)
            cad += fill(stround(dih[5], 4), 9)
            cad += fill(stround(dih[6], 4), 9)
            cad += fill(stround(dih[7], 4), 9)
            cad += fill(stround(dih[8], 4), 9)
            cad += fill(stround(dih[9], 4), 9)
            cad += fill(stround(dih[10], 4), 9)
            cad += fill('    ;'+dih[11], 40, rev=True) + '\n'
        
        for dih in int_na:
            cad += ';'
            cad += fill(str(self.aua_idxs[dih[0]]), 4)
            cad += fill(str(self.aua_idxs[dih[1]]), 4)
            cad += fill(str(self.aua_idxs[dih[2]]), 4)
            cad += fill(str(self.aua_idxs[dih[3]]), 4)
            cad += '; Missing dihedral\n'
        
        return cad
    
    def WriteVS(self):
        cad = '\n[  virtual_sites2  ]\n'
        for vs in self.vs2:
            cad += fill(str(vs[0]), 4)
            cad += fill(str(self.aua_idxs[vs[1]]), 4)
            cad += fill(str(self.aua_idxs[vs[2]]), 4)
            cad += fill(str(vs[3]), 4)
            cad += fill(str(vs[4]), 9)
            cad += '\n'
            
        for vs in self.vs2_:
            cad += fill(str(vs[0]), 4)
            cad += fill(str(self.aua_idxs[vs[1]]), 4)
            cad += fill(str(vs[2]), 4)
            cad += fill(str(vs[3]), 4)
            cad += fill(str(vs[4]), 9)
            cad += '\n'
            
        cad += '\n[  virtual_sites3  ]\n'
        for vs in self.vs3:
            cad += fill(str(vs[0]), 4)
            cad += fill(str(self.aua_idxs[vs[1]]), 4)
            cad += fill(str(self.aua_idxs[vs[2]]), 4)
            cad += fill(str(self.aua_idxs[vs[3]]), 4)
            cad += fill(str(vs[4]), 4)
            cad += fill(str(vs[5]), 9)
            cad += fill(str(vs[6]), 9)
            cad += '\n'
            
        cad += '\n[  virtual_sitesn  ]\n'
        for vs in self.vsn:
            for i in vs:
                cad += fill(str(i), 4)
            cad += '\n'
            
        return cad
    
    def GenExclusions(self, int12s, int13s, int14s):
        for atom in self.atoms:
        
            l_exc = []
            
            if atom.sig != 0:
                if not atom.hasvs:
                    auaidx = self.aua_idxs[atom.idx]
                elif atom.hasvs:        
                    auaidx = self.atomidx_vs[atom.idx]
                    
                    
                for int12 in int12s:
                    if atom.idx == int12[0]:
                        pair_idx = int12[1]
                        coatom = self.atoms[pair_idx-1]
                        
                        if coatom.sig != 0 and not coatom.hasvs:
                            auacoidx = self.aua_idxs[coatom.idx]
                        elif coatom.hasvs:
                            auacoidx = self.atomidx_vs[coatom.idx]
                        
                        l_exc.append(auacoidx)
                        
                    elif atom.idx == int12[1]:
                        pair_idx = int12[0]
                        coatom = self.atoms[pair_idx-1]
                        
                        if coatom.sig != 0 and not coatom.hasvs:
                            auacoidx = self.aua_idxs[coatom.idx]
                        elif coatom.hasvs:
                            auacoidx = self.atomidx_vs[coatom.idx]
                            
                        l_exc.append(auacoidx)
                        
                
                for int13 in int13s:
                    if atom.idx == int13[0]:
                        pair_idx = int13[1]
                        coatom = self.atoms[pair_idx-1]
                        
                        if coatom.sig != 0 and not coatom.hasvs:
                            auacoidx = self.aua_idxs[coatom.idx]
                        elif coatom.hasvs:
                            auacoidx = self.atomidx_vs[coatom.idx]
                        
                        l_exc.append(auacoidx)
                        
                    elif atom.idx == int13[1]:
                        pair_idx = int13[0]
                        coatom = self.atoms[pair_idx-1]
                        
                        if coatom.sig != 0 and not coatom.hasvs:
                            auacoidx = self.aua_idxs[coatom.idx]
                        elif coatom.hasvs:
                            auacoidx = self.atomidx_vs[coatom.idx]
                    
                        l_exc.append(auacoidx)
                       
                
                for int14 in int14s:
                    if atom.idx == int14[0]:
                        pair_idx = int14[1]
                        coatom = self.atoms[pair_idx-1]
                        
                        if coatom.sig != 0 and not coatom.hasvs:
                            auacoidx = self.aua_idxs[coatom.idx]
                        elif coatom.hasvs:
                            auacoidx = self.atomidx_vs[coatom.idx]
                        
                        l_exc.append(auacoidx)
                        
                    elif atom.idx == int14[1]:
                        pair_idx = int14[0]
                        coatom = self.atoms[pair_idx-1]
                        
                        if coatom.sig != 0 and not coatom.hasvs:
                            auacoidx = self.aua_idxs[coatom.idx]
                        elif coatom.hasvs:
                            auacoidx = self.atomidx_vs[coatom.idx]
                    
                        l_exc.append(auacoidx)
        
                self.exc.update({auaidx : list(set(l_exc))})        
        
        l_aux_exc = list(self.aux_vs.values())
        l_norm_exc = list(self.exc.keys())
        
        for k_exc, v_exc in self.exc.items():
            self.exc.update({k_exc : v_exc + l_aux_exc})
            
        for auaidx in self.aux_vs.values():
            l_exc = l_norm_exc + l_aux_exc 
            l_exc.remove(auaidx)
            self.exc.update({auaidx : l_exc})
        
    
    def WriteExclusions(self):
        cad = '\n[  exclusions  ]\n'
        
        for idx1, idxs in self.exc.items():
            cad += fill(str(idx1), 4)
            cad += 2*' '
            for idx in idxs:
                cad += fill(str(idx), 4)
            cad += '\n'
        
        return cad
    
    def WritePairs(self, pairs):
        cad = '\n[  pairs  ]\n'
        
        for pair in list(pairs):
            cad += fill(str(self.aua_idxs[pair[0]]), 4)
            cad += fill(str(self.aua_idxs[pair[1]]), 4)
            cad += fill('1', 4) + '\n'
    
        return cad
    
    def WriteSystem(self):
        cad = '\n[  system  ]\n'
        cad += self.short_name + '\n'
        cad += '\n[  molecules  ]\n'
        cad += self.short_name + fill('_NMOL_', 8) + '\n'
        return cad
        
    
    def WriteGro(self, vec, box):
        cad = '\n'
        n_sites = str(len(self.aua_idxs) + len(self.vs_defs))
        cad += n_sites + '\n'
        
        for idx, aua_idx in self.aua_idxs.items():
            atom = self.atoms[idx-1]
            cad += fill('1'+self.short_name, 8)
            cad += fill(atom.name_gro, 7)
            cad += fill(str(aua_idx), 5)
            x = stround(atom.patom.coords[0]/10 + vec[0], 3)
            y = stround(atom.patom.coords[1]/10 + vec[1], 3)
            z = stround(atom.patom.coords[2]/10 + vec[2], 3)
            cad += fill(x, 8) + fill(y, 8) + fill(z, 8) + '\n'
            
        for aua_idx, vs_def in self.vs_defs.items():
            cad += fill('1'+self.short_name, 8)
            cad += fill(vs_def[0], 7)
            cad += fill(str(aua_idx), 5)
            x = stround(self.vs_pos[aua_idx][0], 3)
            y = stround(self.vs_pos[aua_idx][1], 3)
            z = stround(self.vs_pos[aua_idx][2], 3)
            cad += fill(x, 8) + fill(y, 8) + fill(z, 8) + '\n'
        
        cad += fill(stround(box[0], 5), 11)
        cad += fill(stround(box[1], 5), 11)
        cad += fill(stround(box[2], 5), 11) + '\n'
        # cad += '   10.00000   10.00000   10.00000\n'
        
        return cad
            
    

def find_height_pyramid(u, v, w, anga, angb, angc):
    
    v_ = np.sqrt(u**2 + w**2 -2*u*w*np.cos(angc*np.pi/180))
    u_ = np.sqrt(v**2 + w**2 -2*v*w*np.cos(angb*np.pi/180))
    w_ = np.sqrt(u**2 + v**2 -2*u*v*np.cos(anga*np.pi/180))       
    
    x = (u_-v+w)*(v-w+u_)
    y = (v_-w+u)*(w-u+v_)
    z = (w_-u+v)*(u-v+w_)
    x_ = (w-u_+v)*(u_+v+w)
    y_ = (u-v_+w)*(v_+u+w)
    z_ = (v-w_+u)*(w_+u+v)
    
    a = np.sqrt(x*y_*z_)
    b = np.sqrt(x_*y*z_)
    c = np.sqrt(x_*y_*z)
    d = np.sqrt(x*y*z)
    
    term = (-a+b+c+d)
    term *= (a-b+c+d)
    term *= (a+b-c+d)
    term *= (a+b+c-d)
    
    vol = np.sqrt(term)/(192*u*v*w)
    
    sp = (u_+v_+w_)/3
    term2 = sp
    term2 *= (sp-u)
    term2 *= (sp-v)
    term2 *= (sp-w)
    
    area = np.sqrt(term2)
    
    h = 3*vol/area
    
    return h


def fill(stroi, n_spaces, rev=False):
    len_str = len(stroi)
    assert len(stroi) <= n_spaces
    if rev:
        return stroi + (n_spaces-len_str)*' '
    else:
        return (n_spaces-len_str)*' ' + stroi    
        
    
def stround(floatoi, decimals):
    if type(floatoi) == type(0):
        num = float(floatoi)
    elif type(floatoi) == type(0.0):
        num = floatoi
    else:
        try:
            num = float(floatoi)
        except:
            raise ValueError
        
    res = str(round(num, decimals))
    res_decimals = len(res.split('.')[-1])
    while res_decimals < decimals:
        res += '0'
        res_decimals = len(res.split('.')[-1])
        
    return res

def compare_list(li, lj):
    lis = set(li)  
    ljs = set(lj)
    
    dif_lis = lis-ljs
    dif_ljs = ljs-lis

    dif_li = list(dif_lis)
    dif_lj = list(dif_ljs)    
    
    # Reverse the tuple and compare if it must be removed
    for i in dif_li:
        ri = i[::-1]
        if ri in ljs:
            dif_lis.remove(i)
    
    for j in dif_lj:
        rj = j[::-1]
        if rj in lis:
            dif_ljs.remove(j)
            
    return dif_lis, dif_ljs
            
        
               