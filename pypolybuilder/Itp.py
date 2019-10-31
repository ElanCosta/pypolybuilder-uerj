'''
Created on Dec 22, 2013

@author: root
'''

from Angle    import Angle
from Dihedral import Dihedral
from Exclusions import Exclusions
import multiprocessing
import numpy as np
import math
import random
from Utils import unit_vector
from Utils import vectorFromAtoms

class Itp(object):
    atom_list      = []
    bond_list      = []
    angle_list     = []
    dihedral_list  = []
    branch_list    = []
    exclusion_list = []
    pair_list = []
    nonbonded_pairs = [] # Used as a filter to eliminate bonded pairs
    cutoff = 1.2
    cutoff2 = cutoff*cutoff

    def __init__(self, path):                    
        pass

    def create_initial_coords(self):
        atom = self.get_atom_list()[0]
        if atom.get_x() is None:
            atom.setposition(1, 1, 1)
        neighbors = []
        neighbors.append(atom)
        previousAtoms = {}
        previousAtoms[atom.get_nr()] = atom
        while len(neighbors) > 0:
            atom = neighbors.pop()
            insertedNeighbors = []
            count = -1
            for neighbor in self.getAtomNeighbors(atom):
                if neighbor != previousAtoms[atom.get_nr()]: insertedNeighbors.append(neighbor)
                if neighbor.get_x() is not None:
                    continue
                count += 1
                previousAtoms[neighbor.get_nr()] = atom
                self.setPositionNeighbor(neighbor, atom, count, insertedNeighbors, previousAtoms[atom.get_nr()])
                neighbors.append(neighbor)


    def getAtomNeighbors(self, atom):
        neighbors = []
        for bond in self.get_bond_list():
            if (bond.contains_atom(atom)):
                neighbors.append(bond.get_other_atom(atom))
        return neighbors




    def setPositionNeighbor(self, neighbor, atom, count, insertedNeighbors, previousAtom):
        bond = 0.10 if ((neighbor.get_atomtype() == 'H') or (neighbor.get_atomtype() == 'HC')) else 0.15
        if previousAtom == atom:
            if count == 0:
                #print ("COUNT 0 - level 0")
                neighbor.setposition(atom.get_x() + bond, atom.get_y(), atom.get_z())
                return
            if count == 1:
                #print ("COUNT 1 - level 0!!!")
                z = 1
                cos = -0.5
                sin = 0.866
                x = atom.get_x() + bond * cos
                y = atom.get_y() + bond * sin

                neighbor.setposition(x, y, z)
                return
            if count == 2:
                #print("COUNT 2 - level 0!!!")
                v1 = unit_vector(vectorFromAtoms(insertedNeighbors[0], atom))
                v2 = unit_vector(vectorFromAtoms(insertedNeighbors[1], atom))
                v3 = unit_vector(np.array(v1) + np.array(v2))
                position = np.array([atom.get_x(), atom.get_y(), atom.get_z()]) - bond * np.array(v3)
                neighbor.setposition(position[0], position[1], position[2])
                return
            if count == 3:
                neighbor.setposition(atom.get_x(), atom.get_y(), atom.get_z() - 0.1 + random.uniform(0, 0.05))
                return

        if count == 0:
            dihedral = self.find_dihedral(previousAtom, atom)
            if dihedral != None:
                #print("COUNT 0a - level 1!!!")
                otherAtom = None
                if dihedral.get_a_1() == neighbor: otherAtom = dihedral.get_a_4()
                if dihedral.get_a_4() == neighbor: otherAtom = dihedral.get_a_1()
                v1 = unit_vector(vectorFromAtoms(otherAtom, previousAtom))
                v2 = unit_vector(vectorFromAtoms(previousAtom, atom))
                v4 = unit_vector(np.cross(v1, v2))
                if dihedral.get_conformation() is not None:
                    v5 = unit_vector(np.cross(v2, v4))
                    if dihedral.get_conformation() == 'cis':
                        v6 = bond * -0.5 * v2 - bond * 0.866 * -v5
                    else:
                        v6 = bond * -0.5 * v2 - bond * 0.866 * v5
                    position = np.array([atom.get_x(), atom.get_y(), atom.get_z()]) + v6
                    neighbor.setposition(position[0], position[1], position[2])
                    return

                v5 = unit_vector(np.cross(v2, v4)) + (0, 0, random.uniform(-0.5, 0.5))
                v6 = bond * -0.5 * v2 - bond * 0.866 * unit_vector(v5)
                position = np.array([atom.get_x(), atom.get_y(), atom.get_z()]) + v6
                neighbor.setposition(position[0], position[1], position[2])
                return

            if dihedral == None:
                #print("COUNT 0b - level 1!!!")
                cos = -0.5
                sin = 0.866
                x = atom.get_x() + bond * cos
                y = atom.get_y() + bond * sin
                z = atom.get_z()
                neighbor.setposition(x, y, z)
                return





        if count == 1:
            #print("COUNT 1 - level 1!!!")
            sp3angle = math.radians(109.5)
            v2 = unit_vector(vectorFromAtoms(previousAtom,atom))
            v3 = unit_vector(vectorFromAtoms(insertedNeighbors[0],atom))
            v4 = unit_vector(-(v2+v3))
            v5 = np.cross(v2,v3)

            position = np.array([atom.get_x(), atom.get_y(), atom.get_z()]) + bond * math.sin(0.5 * sp3angle) * v5 + bond * math.cos(0.5*sp3angle) * v4
            neighbor.setposition(position[0], position[1], position[2])
            return



        if count == 2:
            #print("COUNT 2 - level 1!!!")
            v1 = unit_vector(vectorFromAtoms(previousAtom, atom))
            v2 = unit_vector(vectorFromAtoms(insertedNeighbors[0], atom))
            v3 = unit_vector(vectorFromAtoms(insertedNeighbors[1], atom))
            v4 = v1 + v2 + v3
            position = np.array([atom.get_x(), atom.get_y(), atom.get_z()]) - bond * np.array(v4)
            neighbor.setposition(position[0], position[1], position[2])
            return

    def create_all_angles(self, bond_list):
        angle_list =  self.angle_list                                          #new angle_list
        for x in range(0, len(bond_list)):
            for y in range(x+1, len(bond_list)):
                a = Angle.makeFromBonds(bond_list[x],bond_list[y])    #create angle
                if (a != None) and (a not in angle_list):           #check if this angle already exists on the angle_list
                    angle_list.append(a)
        return angle_list
    
    def create_all_dihedrals(self, bond_list):
        dihedral_list = self.dihedral_list                                                         #new angle_list
        for atom in self.get_atom_list():
            d = Dihedral.newMake(atom)
            if(len(d) >= 1):
                dihedral_list.extend(d)
            #if(d != None):
            #    dihedral_list.append(d)
        return dihedral_list

    def clear_repeated_dihedrals(self):
        d_list = self.get_dihedral_list()
        remove_list = []
        for x in range(0, len(d_list)-1):
            d = d_list[x]
            for y in range(x+1, len(d_list)-1):
                d2 = d_list[y]
                
                if(d2.containsAtom(d.get_a_1()) and d2.containsAtom(d.get_a_2()) and d2.containsAtom(d.get_a_3()) and d2.containsAtom(d.get_a_4())):
                    if(d2.isDeterminated == False):
                        remove_list.append(d2)

        for d in remove_list:
            if d in d_list: d_list.remove(d)

    def clear_repeated_angles(self):
        a_list = self.get_angle_list()
        remove_list = []
        for x in range(0, len(a_list) - 1):
            a = a_list[x]
            for y in range(x + 1, len(a_list) - 1):
                a2 = a_list[y]

                if (a2.contains_atom(a.get_a_1()) and a2.contains_atom(a.get_a_2()) and a2.contains_atom(
                        a.get_a_3())):
                    if (a2.isDeterminated == False):
                        remove_list.append(a2)

        for a in remove_list:
            if a in a_list: a_list.remove(a)

    def generate_exclusions(self):
        for atom in self.get_atom_list():
            atom.__exclusion_list = list()

        for atom in self.get_atom_list():
            exclusions = list()
            for neighbor in atom.get_atom_list():
                exclusions.append(int(neighbor.get_nr()))
                for second_neighbor in neighbor.get_atom_list():
                    exclusions.append(int(second_neighbor.get_nr()))

            exclusions = exclusions + atom.get_exclusion_list()
            exclusion_list = Exclusions(exclusions)
            atom.set_exclusion_list(exclusion_list)

                
    
    def createDihedral(self,bond_list,x,y,z,dihedral_list):
        d = Dihedral.makeFromBonds(self,bond_list[x],bond_list[y],bond_list[z])    #create angle
        if (d != None) and (d not in dihedral_list):                        #check if this angle already exists on the angle_list
            dihedral_list.append(d)
      #  print(a1.get_nr() + a2.get_nr() + a3.get_nr() + a4.get_nr() + a5.get_nr() + a6.get_nr())
    #NAO SEI SE FUNCIONA
    def get_common_atoms(self,bond_list):
        common_atoms = []
        i = 0
        len = 0
        for bond in bond_list:
            atom = bond[i]
            for bond2  in bond_list:
                if(bond2.contains_atom(atom) == False):
                    break
                elif(len == bond_list.lenght):
                    common_atoms.append(atom)
            len = len+1
            i   = i+1
                    
        return common_atoms



    def createPairList(self):
        self.pair_list = list()
        for atom in self.get_atom_list():
            for neighbor in atom.get_neighbor_list():
                if (atom != neighbor and not self.is_bond(atom,neighbor) and not self.pair_belongs_to_angle(atom,neighbor)):
                    pair = [atom,neighbor]
                    self.pair_list.append(pair)

    def calculateAtomsNeighbors(self):
        for atom in self.get_atom_list():
            atom.neighbor_list = list()

        for i in range(0, len(self.get_atom_list())):
            atom1  = self.get_atom_list()[i]
            atom1.neighbor_list = list()
            for j in range(i, len(self.get_atom_list())):
                atom2 = self.get_atom_list()[j]
                atom1.calcNeighbor(atom2,self.cutoff2)
     
    def createNonBondedPairs(self):
        self.nonbonded_pairs = list()
        for atom1 in self.get_atom_list():
            for atom2 in atom1.get_neighbor_list():
                if (atom1 != atom2 and not (atom2 in (atom1.connected_list))):
                    pair = [atom1, atom2]
                    self.nonbonded_pairs.append(pair)

    def standardPairList(self):
        self.pair_list = list()
        for pair in self.get_nonbonded_pairs():
            atom1 = pair[0]
            atom2 = pair[1]
            if(atom1.testAtomIsClose(atom2,self.cutoff2)):
                self.pair_list.append(pair)
            
            
        

    #Determina se dois atomos formam um bond
    #@classmethod
    def is_bond(self,atom_1,atom_2):
        for bond in self.get_bond_list():
            if(atom_1.get_nr() != atom_2.get_nr() and bond.contains_atom(atom_1) and bond.contains_atom(atom_2) ):
                return True
        return False

    def get_bond_from_atoms(self, atom_1, atom_2):
        for bond in self.get_bond_list():
            if(atom_1.get_nr() != atom_2.get_nr() and bond.contains_atom(atom_1) and bond.contains_atom(atom_2) ):
                return bond



    def pair_belongs_to_angle(self,atom_1,atom_2):
        for angle in self.get_angle_list():
            if (atom_1.get_nr() != atom_2.get_nr() and angle.contains_atom(atom_1) and angle.contains_atom(atom_2)):
                return True
        return False

    def pair_belongs_to_dihedral(self,atom_1,atom_2):
        for dihedral in self.get_dihedral_list():
            if(atom_1.get_nr() != atom_2.get_nr() and dihedral.containsAtom(atom_1) and dihedral.containsAtom(atom_2) ):
                return True
        return False

    def quartet_belongs_to_dihedral(self, atom_1, atom_2, atom_3, atom_4):
        for dihedral in self.get_dihedral_list():
            if(atom_1.get_nr() != atom_2.get_nr() and atom_1.get_nr() != atom_3.get_nr() and atom_1.get_nr() != atom_4.get_nr() and
                       atom_2.get_nr() != atom_3.get_nr() and atom_2.get_nr() != atom_4.get_nr() and atom_3.get_nr() != atom_4.get_nr()
               and dihedral.containsAtom(atom_1) and dihedral.containsAtom(atom_2) and dihedral.containsAtom(atom_3) and dihedral.containsAtom(atom_4) ):
                return True
        return False

    def triplet_belongs_to_angle(self, atom_1, atom_2, atom_3):
        for angle in self.get_angle_list():
            if(atom_1.get_nr() != atom_2.get_nr() and atom_1.get_nr() != atom_3.get_nr() and atom_2.get_nr() != atom_3.get_nr()
               and angle.contains_atom(atom_1) and angle.contains_atom(atom_2) and angle.contains_atom(atom_3)):
                return True
        return False


    def pair_defines_rotation_axis(self,atom_1,atom_2):
        for dihedral in self.get_dihedral_list():
            if(dihedral.containsAtom(atom_1) and dihedral.containsAtom(atom_2) and atom_1.get_nr() != atom_2.get_nr()):
                if(dihedral.get_func() == 1 ):
                    if((dihedral.get_a_2().get_nr() == atom_1.get_nr() and dihedral.get_a_3().get_nr() == atom_2.get_nr()) or (dihedral.get_a_2().get_nr() == atom_2.get_nr() and dihedral.get_a_3().get_nr() == atom_1.get_nr())):
                        return True
        return False

    def returnFourthAtom(self, a1, a2, a3):
        for dihedral in self.get_dihedral_list():
            if (dihedral.containsAtom(a1) and dihedral.containsAtom(a2) and dihedral.containsAtom(a3)):
                if (dihedral.get_a_1() !=  a1 and dihedral.get_a_1() !=  a2 and dihedral.get_a_1() !=  a3):
                    return dihedral.get_a_1()
                if (dihedral.get_a_2() !=  a1 and dihedral.get_a_2() !=  a2 and dihedral.get_a_2() !=  a3):
                    return dihedral.get_a_2()
                if (dihedral.get_a_3() !=  a1 and dihedral.get_a_3() !=  a2 and dihedral.get_a_3() !=  a3):
                    return dihedral.get_a_3()


    
    @classmethod
    def get_repeated_atoms(self,atom_list):
        repeated_atoms = []
        
        
        for x in range(0, len(atom_list)):
            for y in range(x+1, len(atom_list)):
                if(atom_list[x].get_nr() == atom_list[y].get_nr()):
                    repeated_atoms.append(atom_list[x])
                    

            
        if(len(repeated_atoms) == 2):
            return repeated_atoms
        
        return None
    
    def find_angle(self,atom1,atom2,atom3):
        for angle in self.get_angle_list():
            if angle.contains_atom(atom1) and angle.contains_atom(atom2) and angle.contains_atom(atom3):
                return angle

    def find_dihedral(self, atom2, atom3):
        for dihedral in self.get_dihedral_list():
            if dihedral.get_a_2() == atom2 and dihedral.get_a_3() == atom3:
                return dihedral
            if dihedral.get_a_3() == atom2 and dihedral.get_a_2() == atom3:
                return dihedral
        return None

    def get_dihedral_from_atoms(self, atom_1, atom_2, atom_3, atom_4):
        for dihedral in self.get_dihedral_list():
            if(atom_1.get_nr() != atom_2.get_nr() and atom_2.get_nr() != atom_3.get_nr() and atom_3.get_nr() != atom_4.get_nr()
               and atom_1.get_nr() != atom_3.get_nr() and dihedral.containsAtom(atom_1) and dihedral.containsAtom(atom_2) and dihedral.containsAtom(atom_3) and dihedral.containsAtom(atom_4)):
                return dihedral


    def print_atoms(self):
        for a in self.get_atom_list():
            print(str(a.get_nr()) + "  " + str(a.get_atomtype()) + " x=" + str(a.get_x())+ " y=" + str(a.get_y())+ " z=" + str(a.get_z()))
            
    def print_branches(self):
        for b in self.branch_list:
            print(b.get_donor())
            print(b.get_acceptor())
    
    #get atom object by his atom_number
    def find_atom(self,atom_number):
        atoms = self.get_atom_list()
        for atom in atoms:
            if (str(atom_number) == str(atom.get_nr())):
                return atom
            
    def print_bonds(self):
        for bond in self.get_bond_list():
            print(str(bond.get_a_1().get_nr()) +  " " + bond.get_a_1().get_atomtype() + " " + str(bond.get_a_2().get_nr()) + " " + bond.get_a_2().get_atomtype())
    
    #pega os bonds que um atomo contem      
    def get_bonds_by_atom(self,atom_number):
        bonds = []
        for bond in self.get_bond_list():
            if(str(bond.get_a_1().get_nr()) == str(atom_number) or str(bond.get_a_2().get_nr()) == str(atom_number)):
                bonds.append(bond)
        return bonds
    
    #Checa se o atom pertence a 2 bonds
    def is_common_atom(self,atom,b1,b2):
        if (b1.contains_atom(atom) and b2.contains_atom(atom)):
            return True
        return False

    def get_bond_list(self):
        return self.__bond_list


    def get_angle_list(self):
        return self.__angle_list


    def get_dihedral_list(self):
        return self.__dihedral_list


    def get_branch_list(self):
        return self.__branch_list


    def set_bond_list(self, value):
        self.__bond_list = value


    def set_angle_list(self, value):
        self.__angle_list = value


    def set_dihedral_list(self, value):
        self.__dihedral_list = value


    def set_branch_list(self, value):
        self.__branch_list = value
        
    def get_atom_list(self):
        return self.__atom_list


    def set_atom_list(self, value):
        self.__atom_list = value

    def get_pair_list(self):
        return self.__pair_list

    def set_pair_list(self, value):
        self.__pair_list = value
    
    def get_exclusion_list(self):
        return self.__exclusion_list

    def set_exclusion_list(self, value):
        self.__exclusion_list = value

    def get_nonbonded_pairs(self):
        return self.__nonbonded_pairs

    def set_nonbonded_pairs(self, value):
        self.__nonbonded_pairs = value

    atom_list = property(get_atom_list, set_atom_list, None, None)
    bond_list = property(get_bond_list, set_bond_list, None, None)
    angle_list = property(get_angle_list, set_angle_list, None, None)
    dihedral_list = property(get_dihedral_list, set_dihedral_list, None, None)
    branch_list = property(get_branch_list, set_branch_list, None, None)
    pair_list = property(get_pair_list, set_pair_list, None, None)
    exclusion_list = property(get_exclusion_list, set_exclusion_list, None, None)
    nonbonded_pairs = property(get_nonbonded_pairs, set_nonbonded_pairs, None, None)
    