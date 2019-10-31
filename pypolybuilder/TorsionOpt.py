'''
Created on Jan 22, 2019

@author: bruno
'''

from __future__ import division
from Printer import Printer
import sys
import time
import numpy as np
import math
import copy
from Polymer import Polymer
# from Dendrimer import Dendrimer
from Utils import unit_vector
from Utils import vectorFromAtoms
from vbga import VBGA
from random import  uniform


Matrix = None


class Zmatrix(object):
    atoms = []
    boundto = []
    bondvalue = []
    angleto = []
    angvalue = []
    dihedralto = []
    dihvalue = []
    def __init__(self, atoms, boundto, bond, angleto, ang, dihedralto, dih):
        self.atoms = atoms
        self.boundto = boundto
        self.bondvalue = bond
        self.angleto = angleto
        self.angvalue = ang
        self.dihedralto = dihedralto
        self.dihvalue = dih

    def printMatrix(self):
        print (len(self.atoms))
        print('{:>4s}'.format(self.atoms[0].get_atomtype()[0]))
        print('{:>4s} {:>4d} {:>11.5f}'.format(self.atoms[1].get_atomtype()[0],
                                                   self.boundto[0], float(self.bondvalue[0])))
        print('{:>4s} {:>4d} {:>11.5f} {:>4d} {:>11.5f}'.format(self.atoms[2].get_atomtype()[0],
                                                   self.boundto[1], float(self.bondvalue[1]),
                                                    self.angleto[0], float(self.angvalue[0])))
        for i in range(3, len(self.atoms)):
            print('{:>4s} {:>4d} {:>11.5f} {:>4d} {:>11.5f} {:>4d} {:>11s}'.format(self.atoms[i].get_atomtype()[0],
                                                                    self.boundto[i-1], float(self.bondvalue[i-1]),
                                                                    self.angleto[i-2], float(self.angvalue[i-2]),
                                                                    self.dihedralto[i-3],str(self.dihvalue[i-3])))

    def randomize(self):
        for i in range(0, len(self.dihvalue)):
            if self.dihvalue[i] == "NULL":
                a = np.random.randint(12)
                newangle = float(a*30)
                self.dihvalue[i] = newangle

    def movableTorsions(self):
        # Improve this at a later stage to include restrictions (maybe need info about dihedrals)
        indexList = []
        for i in range(0, len(self.dihvalue)):
            if self.dihvalue[i] == "NULL":
                indexList.append(i)
        return (indexList)




    def converToXYZ(self):
        natoms = len(self.atoms)
        #print (natoms)
        #print ('XYZ')

        xyzarr = np.zeros([natoms, 3])
        if (natoms > 1):
            xyzarr[1] = [self.bondvalue[0], 0.0, 0.0]

        if (natoms > 2):
            i = self.boundto[1] - 1
            j = self.angleto[0] - 1
            r = self.bondvalue[1]
            theta = self.angvalue[0] * np.pi / 180.0
            x = r * np.cos(theta)
            y = r * np.sin(theta)
            a_i = xyzarr[i]
            b_ij = xyzarr[j] - xyzarr[i]
            if (b_ij[0] < 0):
                x = a_i[0] - x
                y = a_i[1] - y
            else:
                x = a_i[0] + x
                y = a_i[1] + y
            xyzarr[2] = [x, y, 0.0]

        for n in range(3, natoms):
            r = self.bondvalue[n - 1]
            theta = self.angvalue[n - 2] * np.pi / 180.0
            phi = self.dihvalue[n - 3] * np.pi / 180.0

            sinTheta = np.sin(theta)
            cosTheta = np.cos(theta)
            sinPhi = np.sin(phi)
            cosPhi = np.cos(phi)

            x = r * cosTheta
            y = r * cosPhi * sinTheta
            z = r * sinPhi * sinTheta

            i = self.boundto[n - 1] - 1
            j = self.angleto[n - 2] - 1
            k = self.dihedralto[n - 3] - 1
            a = xyzarr[k]
            b = xyzarr[j]
            c = xyzarr[i]

            #print("DEBUG: ", i, " ", j," ", k)

            ab = b - a
            bc = c - b
            bc = bc / np.linalg.norm(bc)
            nv = np.cross(ab, bc)
            nv = nv / np.linalg.norm(nv)
            ncbc = np.cross(nv, bc)

            new_x = c[0] - bc[0] * x + ncbc[0] * y + nv[0] * z
            new_y = c[1] - bc[1] * x + ncbc[1] * y + nv[1] * z
            new_z = c[2] - bc[2] * x + ncbc[2] * y + nv[2] * z
            xyzarr[n] = [new_x, new_y, new_z]


        return(xyzarr)

class TorsionOpt(object):
    def __init__(self, top):
        self.top = top

###### busca otimizada em dois loops
    def topToZmatrix(self):
        Zatom = []
        ZbondTo = []
        ZangleTo = []
        ZdihedralTo = []
        ZbondValue = []
        ZangleValue = []
        ZdihValue = []
        # lets first get the elements and atoms
        count = 0
        bonddist = [1, 1, 1.09, 1.12, 1.23, 1.25, 1.32, 1.33, 1.33, 1.33, 1.34, 1.34, 1.36, 1.38, 1.39, 1.39, 1.4, 1.43,
                    1.43, 1.435, 1.47, 1.48, 1.48, 1.48, 1.5, 1.52, 1.53, 1.61, 1.63, 1.78, 1.78, 1.83, 1.87, 1.98, 2,
                    2.04, 2.21, 1, 1.1, 1.758, 1.53, 1.93799, 1.76, 1.265, 1.35, 1.63299, 2.33839, 2.90283, 2.79388,
                    2.91189, 2.077, 2.87407]
        bondangle = [90.000000, 90.000000, 96.000000, 100.000000, 103.000000, 104.000000, 108.000000, 109.500000,
                     109.500000, 109.500000, 109.500000, 109.500000, 109.500000, 109.600000, 111.000000, 113.000000,
                     115.000000, 115.000000, 115.000000, 116.000000, 116.000000, 117.000000, 120.000000, 120.000000,
                     120.000000, 120.000000, 120.000000, 120.000000, 120.000000, 121.000000, 122.000000, 123.000000,
                     124.000000, 125.000000, 125.000000, 126.000000, 126.000000, 126.000000, 132.000000, 155.000000,
                     180.000000, 109.500000, 107.570000, 111.300000, 97.400000, 106.750000, 108.530000, 109.500000,
                     107.600000, 109.500000, 110.300000, 111.400000, 117.200000, 121.400000]
        improper = [1000, 0.0, 35.26439, 0.0, 180.0, -35.26439]
        # import csv
        # c = csv.writer(open("MatrizZ.csv", "w"))
        f = open("MatrizZ.csv", "w", encoding="utf8")
        #### inicializa a matriz q vai ser impressa
        f.write("Atom,BondTo,BondValue,AngleTo,AngleValue,DihedralTo,DihValue")
        f.write("\n")

        if len(self.top.get_atom_list()) > 3:

            # print first atom
            # print('{:>4s}'.format(self.top.get_atom_list()[0].get_atomtype()[0]))
            Zatom.append(self.top.get_atom_list()[0])
            f.write(str(self.top.get_atom_list()[0]))
            f.write(str(self.top.get_atom_list()[0]))
            f.write(",")
            f.write("")
            f.write(",")
            f.write("")
            f.write(",")
            f.write("")
            f.write(",")
            f.write("")
            f.write(",")
            f.write("")
            f.write(",")
            f.write("")
            f.write("\n")

            # print second atom -- BOND
            if (self.top.is_bond(self.top.get_atom_list()[0], self.top.get_atom_list()[1])):
                bond = self.top.get_bond_from_atoms(self.top.get_atom_list()[0], self.top.get_atom_list()[1])
                indexbonddist = int(bond.get_param()[3:])  # Pra pegar so o numero e evitar o gb_
                # print('{:>4s} {:>4s} {:>11.5f}'.format(self.top.get_atom_list()[1].get_atomtype()[0],
                #                                       self.top.get_atom_list()[0].get_nr(), bonddist[indexbonddist]))
                Zatom.append(self.top.get_atom_list()[1])
                ZbondTo.append(int(self.top.get_atom_list()[0].get_nr()))
                ZbondValue.append(bonddist[indexbonddist] / 10)  # divide por 10 pra ser angstrom
                f.write(str(self.top.get_atom_list()[1]))
                f.write(",")
                f.write(str(int(self.top.get_atom_list()[0].get_nr())))
                f.write(",")
                f.write(str(bonddist[indexbonddist] / 10))
                f.write(",")
                f.write("")
                f.write(",")
                f.write("")
                f.write(",")
                f.write("")
                f.write(",")
                f.write("")
                f.write("\n")
            # print third atom -- BOND and ANGLE
            if (self.top.is_bond(self.top.get_atom_list()[0], self.top.get_atom_list()[2])):
                bond = self.top.get_bond_from_atoms(self.top.get_atom_list()[0], self.top.get_atom_list()[2])
                indexbonddist = int(bond.get_param()[3:])
                angle = self.top.find_angle(self.top.get_atom_list()[0], self.top.get_atom_list()[1],
                                            self.top.get_atom_list()[2])
                indexangle = int(angle.get_param()[3:])
                # print('{:>4s} {:>4s} {:>11.5f} {:>4s} {:>11.5f}'.format(self.top.get_atom_list()[2].get_atomtype()[0],
                #                                       self.top.get_atom_list()[0].get_nr(), bonddist[indexbonddist],
                #                                       self.top.get_atom_list()[1].get_nr(), bondangle[indexangle]))
                Zatom.append(self.top.get_atom_list()[2])
                ZbondTo.append(int(self.top.get_atom_list()[0].get_nr()))
                ZbondValue.append(bonddist[indexbonddist] / 10)
                ZangleTo.append(int(self.top.get_atom_list()[1].get_nr()))
                ZangleValue.append(bondangle[indexangle])

                f.write(str(self.top.get_atom_list()[2]))
                f.write(",")
                f.write(str(int(self.top.get_atom_list()[1].get_nr())))
                f.write(",")
                f.write(str(bonddist[indexbonddist] / 10))
                f.write(",")
                f.write(str(int(self.top.get_atom_list()[0].get_nr())))
                f.write(",")
                f.write(str(bondangle[indexangle]))
                f.write(",")
                f.write("")
                f.write(",")
                f.write("")
                f.write("\n")

            else:
                bond = self.top.get_bond_from_atoms(self.top.get_atom_list()[0], self.top.get_atom_list()[2])
                indexbonddist = int(bond.get_param()[3:])
                angle = self.top.find_angle(self.top.get_atom_list()[0], self.top.get_atom_list()[1],
                                            self.top.get_atom_list()[2])
                indexangle = int(angle.get_param()[3:])
                # print('{:>4s} {:>4s} {:>11.5f}'.format(self.top.get_atom_list()[2].get_atomtype()[0],
                #                                       self.top.get_atom_list()[1].get_nr(), bonddist[indexbonddist],
                #                                       self.top.get_atom_list()[0].get_nr(), bondangle[indexangle]))
                Zatom.append(self.top.get_atom_list()[2])
                ZbondTo.append(int(self.top.get_atom_list()[1].get_nr()))
                ZbondValue.append(bonddist[indexbonddist] / 10)
                ZangleTo.append(int(self.top.get_atom_list()[0].get_nr()))
                ZangleValue.append(bondangle[indexangle])

                f.write(str(self.top.get_atom_list()[2]))
                f.write(",")
                f.write(str(int(self.top.get_atom_list()[1].get_nr())))
                f.write(",")
                f.write(str(bonddist[indexbonddist] / 10))
                f.write(",")
                f.write(str(int(self.top.get_atom_list()[0].get_nr())))
                f.write(",")
                f.write(str(bondangle[indexangle]))
                f.write(",")
                f.write("")
                f.write(",")
                f.write("")
                f.write("\n")
            # print all other atoms

            for i in range(3, len(self.top.get_atom_list())):
                j = i - 1  # bond
                # k = i - 2 #angle
                l = 0  # dihedral
                found = False
                while j <= i:  # len(self.top.get_atom_list()):
                    # print ("DEBUG : ", i , " ", j)
                    if (self.top.is_bond(self.top.get_atom_list()[i], self.top.get_atom_list()[j]) and found == False):
                        k = j - 1
                        while k <= i:  # len(self.top.get_atom_list()):
                            # print ("DEBUG : ", i, " ", j, " ", k)
                            if (self.top.triplet_belongs_to_angle(self.top.get_atom_list()[i],
                                                                  self.top.get_atom_list()[j],
                                                                  self.top.get_atom_list()[k]) and found == False):
                                # l = k - 1
                                while l < i:  # len(self.top.get_atom_list()):
                                    if (self.top.quartet_belongs_to_dihedral(self.top.get_atom_list()[i],
                                                                             self.top.get_atom_list()[j],
                                                                             self.top.get_atom_list()[k],
                                                                             self.top.get_atom_list()[
                                                                                 l]) and found == False):
                                        bond = self.top.get_bond_from_atoms(self.top.get_atom_list()[i],
                                                                            self.top.get_atom_list()[j])
                                        indexbonddist = int(bond.get_param()[3:])
                                        angle = self.top.find_angle(self.top.get_atom_list()[i],
                                                                    self.top.get_atom_list()[j],
                                                                    self.top.get_atom_list()[k])
                                        # print (angle, " ", i, " ", j, " ", k)
                                        indexangle = int(angle.get_param()[3:])
                                        dihedral = self.top.get_dihedral_from_atoms(self.top.get_atom_list()[i],
                                                                                    self.top.get_atom_list()[j],
                                                                                    self.top.get_atom_list()[k],
                                                                                    self.top.get_atom_list()[l])
                                        # print (dihedral, " ", i, " ", j, " ", k, " ", l)
                                        indexdihedral = int(dihedral.get_param()[3:])
                                        dihedralvalue = None
                                        if (dihedral.get_func() == 1):
                                            # dihedralvalue = torsion[indexdihedral]
                                            dihedralvalue = "NULL"
                                            # print ('{:>4s} {:>4d} {:>11.5f} {:>4d} {:>11.5f} {:>4d} {:>11s}'.format(
                                            #    self.top.get_atom_list()[i].get_atomtype()[0],
                                            #    int(self.top.get_atom_list()[j].get_nr()), bonddist[indexbonddist],
                                            #    int(self.top.get_atom_list()[k].get_nr()), bondangle[indexangle],
                                            #    int(self.top.get_atom_list()[l].get_nr()), dihedralvalue))
                                            Zatom.append(self.top.get_atom_list()[i])
                                            ZbondTo.append(int(self.top.get_atom_list()[j].get_nr()))
                                            ZbondValue.append(bonddist[indexbonddist] / 10)
                                            ZangleTo.append(int(self.top.get_atom_list()[k].get_nr()))
                                            ZangleValue.append(bondangle[indexangle])
                                            ZdihedralTo.append(int(self.top.get_atom_list()[l].get_nr()))
                                            ZdihValue.append(dihedralvalue)

                                            # cria um arquivo .csv com a matrizZ. o nome do arquivo é MatrizZ.csv
                                            f.write(str(self.top.get_atom_list()[i]))
                                            f.write(",")
                                            f.write(str(int(self.top.get_atom_list()[j].get_nr())))
                                            f.write(",")
                                            f.write(str(bonddist[indexbonddist] / 10))
                                            f.write(",")
                                            f.write(str(int(self.top.get_atom_list()[k].get_nr())))
                                            f.write(",")
                                            f.write(str(bondangle[indexangle]))
                                            f.write(",")
                                            f.write(str(int(self.top.get_atom_list()[l].get_nr())))
                                            f.write(",")
                                            f.write(str(dihedralvalue))
                                            f.write("\n")
                                        else:
                                            dihedralvalue = improper[indexdihedral]
                                            # print ('{:>4s} {:>4d} {:>11.5f} {:>4d} {:>11.5f} {:>4d} {:>11.5f}'.format(self.top.get_atom_list()[i].get_atomtype()[0],
                                            #    int(self.top.get_atom_list()[j].get_nr()), bonddist[indexbonddist], int(self.top.get_atom_list()[k].get_nr()), bondangle[indexangle],
                                            #    int(self.top.get_atom_list()[l].get_nr()),  dihedralvalue))
                                            Zatom.append(self.top.get_atom_list()[i])
                                            ZbondTo.append(int(self.top.get_atom_list()[j].get_nr()))
                                            ZbondValue.append(bonddist[indexbonddist] / 10)
                                            ZangleTo.append(int(self.top.get_atom_list()[k].get_nr()))
                                            ZangleValue.append(bondangle[indexangle])
                                            ZdihedralTo.append(int(self.top.get_atom_list()[l].get_nr()))
                                            ZdihValue.append(dihedralvalue)

                                            # cria um arquivo .csv com a matrizZ. o nome do arquivo é MatrizZ.csv

                                            f.write(str(self.top.get_atom_list()[i]))
                                            f.write(",")
                                            f.write(str(int(self.top.get_atom_list()[j].get_nr())))
                                            f.write(",")
                                            f.write(str(bonddist[indexbonddist] / 10))
                                            f.write(",")
                                            f.write(str(int(self.top.get_atom_list()[k].get_nr())))
                                            f.write(",")
                                            f.write(str(bondangle[indexangle]))
                                            f.write(",")
                                            f.write(str(int(self.top.get_atom_list()[l].get_nr())))
                                            f.write(",")
                                            f.write(str(dihedralvalue))
                                            f.write("\n")

                                        # print ("DEBUG: ", self.top.get_atom_list()[i].get_nr(), " ", self.top.get_atom_list()[j].get_nr(), " ",
                                        #       self.top.get_atom_list()[k].get_nr(), " ",
                                        #       self.top.get_atom_list()[l].get_nr(), " ", dihedral.get_func(), " i: ", i)
                                        found = True
                                        break
                                    l = l + 1
                                break
                            k = k - 1
                        break
                    j = j - 1
            # print(type(Zatom[1]))
            # print(type(ZbondTo[1]))
            # print(type(ZbondValue[1]))
            # print(type(ZangleTo[1]))
            # print(type(ZangleValue[1]))
            # print(type(ZdihedralTo[1]))
            # print(type(ZdihValue[1]))
            return (Zatom, ZbondTo, ZbondValue, ZangleTo, ZangleValue, ZdihedralTo, ZdihValue)

###### Importa a matriz Z
    def ImportZmatrix(self):
        import pandas as pd
        Data = pd.read_csv("MatrizZ.csv", index_col=None)
        # print(Data)
        Atom = list(Data.Atom)
        BondTo = list(Data.BondTo)
        BondValue = list(Data.BondValue)
        AngleTo = list(Data.AngleTo)
        AngleValue = list(Data.AngleValue)
        DihedralTo = list(Data.DihedralTo)
        DihValue = list(Data.DihValue)
        Zatom = []
        ZbondTo = []
        ZbondValue = []
        ZangleTo = []
        ZangleValue = []
        ZdihedralTo = []
        ZdihValue = []
        for i in range(0, len(Atom)):
            Zatom.append(self.top.get_atom_list()[i])
            if i >= 1:
                ZbondTo.append(int(BondTo[i]))
                ZbondValue.append(BondValue[i])
            if i >= 2:
                ZangleTo.append(int(AngleTo[i]))
                ZangleValue.append(AngleValue[i])
            if i >= 3:
                ZdihedralTo.append(int(DihedralTo[i]))
                if type(DihValue[i])== int:
                    ZdihValue.append(DihedralTo[i])
                else:
                    dihedralvalue = "NULL"
                    ZdihValue.append(dihedralvalue)


        #print(Zatom)
        #print(ZbondTo)
        #print(ZbondValue)
        #print(ZangleTo)
        #print(ZangleValue)
        #print(ZdihedralTo)
        #print(ZdihValue)

        return (Zatom, ZbondTo, ZbondValue, ZangleTo, ZangleValue, ZdihedralTo, ZdihValue)


    def addXYZtoTOP(self,xyz):
        n = 0
        for atom in self.top.get_atom_list():
            atom.set_x(xyz[n][0])
            atom.set_y(xyz[n][1])
            atom.set_z(xyz[n][2])
            n = n + 1

    def run(self):
        start = time.time()
        print("Converting Top to Z-matrix ...")
        a, b, c, d, e, f, g = self.topToZmatrix()
        end = time.time()
        print("Elapsed time: ", end - start)
        global Matrix
        Matrix = Zmatrix(a, b, c, d, e, f, g)
        # GLOBAL MATRIX:

        #Matrix.printMatrix()
        #Matrix.randomize()
        #fitness = Matrix.calcRgy()
        #print ("Fitness: ", fitness)

        myGA = VBGA(Individual, calcRgy, DefaultSelectionMethod, UniformCrossoverMethod, 10.0, DefaultMutation, 20.0, populationSize=25)
        myGA.run(200)
        xyz = myGA.getBest().matrix.converToXYZ()
        self.addXYZtoTOP(xyz)


class Individual:

    def __init__(self):
        self.matrix = copy.deepcopy(Matrix)
        self.fitValue = 0

    def randomize(self):
        for i in range(0,len(self.matrix.dihvalue)):
            if self.matrix.dihvalue[i] == "NULL":
                a = np.random.randint(24)
                newangle = float(a*15)
                self.matrix.dihvalue[i] = newangle



def DefaultSelectionMethod(population, selectionSize=10):
    population.sort(key=lambda x: x.fitValue, reverse=True)
    return population[:selectionSize]  # return 0,1,2....selectionSize

def DefaultCrossoverMethod(father, mother, childOne, childTwo):
    child.x = (father.x + mother.x) / 2
    child.y = (father.y + mother.y) / 2

    return child

def UniformCrossoverMethod(father, mother, childOne, childTwo):
    for i in range(0,len(father.matrix.dihvalue)):
        childOne.matrix.dihvalue[i] = (father.matrix.dihvalue[i] + mother.matrix.dihvalue[i])/2
        childTwo.matrix.dihvalue[i] = (father.matrix.dihvalue[i] + mother.matrix.dihvalue[i])/3
    return [childOne, childTwo]

def DefaultMutation(individual):
    probability = 0.1
    for i in range(0,len(individual.matrix.dihvalue)):
        dice = uniform(0,1)
        #print("DEBUG -- prob and dice ", probability, " ", dice)
        if dice < probability:
            a = np.random.randint(12)
            newangle = float(a * 30)
            individual.matrix.dihvalue[i] = newangle
    return individual


def calcRgy(individual):
    natoms = len(individual.matrix.atoms)
    xyzarr = individual.matrix.converToXYZ()
    x, y, z = 0.0, 0.0, 0.0
    for i in range(natoms):
        x += xyzarr[i][0]
        y += xyzarr[i][1]
        z += xyzarr[i][2]
        # print ('{:<4s}\t{:>11.5f}\t{:>11.5f}\t{:>11.5f}'.format(self.atoms[i].get_atomtype()[0], xyzarr[i][0], xyzarr[i][1], xyzarr[i][2]))
    centerOfGeom = np.array([x, y, z]) / float(natoms)
    rcum = 0.0
    for i in range(natoms):
        rcum += np.linalg.norm(xyzarr[i] - centerOfGeom)
    return (rcum / natoms)

def calcDistSum(individual):
    natoms = len(individual.matrix.atoms)
    xyzarr = individual.matrix.converToXYZ()
    x, y, z = 0.0, 0.0, 0.0
    distSum = 0.0
    for i in range(natoms):
        x_i = xyzarr[i][0]
        y_i = xyzarr[i][1]
        z_i = xyzarr[i][2]
        for j in range(i,natoms):
            x_j = xyzarr[j][0]
            y_j = xyzarr[j][1]
            z_j = xyzarr[j][2]
            dist = (x_i - x_j) + (y_i - y_j) + (z_i - z_j)
            dist = dist * dist
            distSum += dist
    return distSum

