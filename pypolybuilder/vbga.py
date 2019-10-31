#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Sep  1 12:13:44 2018

@author: vitor
"""

import copy
from random import uniform
                
class VBGA:
    def __init__(self, individualClass,fitnessMethod,selectionMethod,crossoverMethod,crossoverRate, mutationMethod, mutationRate, populationSize):
        self.populationSize = populationSize
        self.individualClass = individualClass
        self.fitnessMethod = fitnessMethod
        self.selectionMethod = selectionMethod
        self.crossoverMethod = crossoverMethod
        self.crossoverRate = crossoverRate
        self.mutationRate = mutationRate
        self.mutationMethod = mutationMethod
        self.population = []
       
    def createRandomIndividual(self):
        ind = self.individualClass()
        ind.randomize()
        return ind
        
    def createInitialPopulation(self):
        self.population = []
        for i in range(0,self.populationSize):
            ind = self.createRandomIndividual()
            self.population.append(ind)
    
    def evaluateIndividual(self, individual):
        individual.fitValue = self.fitnessMethod(individual)
    
    def evaluatePopulation(self):
        for individual in self.population:
            self.evaluateIndividual(individual)

    def averageFitness(self):
        max = 0
        for individual in self.population:
            #print(individual.fitValue)
            if individual.fitValue > max:
                max = individual.fitValue
        return max

    def applySelectionMethod(self):
        selectedPopulation = self.selectionMethod(self.population)
        return selectedPopulation

    def fillPopulation(self, selectedPopulation):
        nextPopulation = copy.deepcopy(selectedPopulation)
        for i in range(len(selectedPopulation),self.populationSize):
            ind = self.createRandomIndividual()
            nextPopulation.append(ind)
        
        self.population = nextPopulation

    def applyCrossoverMethod(self, selected):
        populationCross = copy.deepcopy(selected)
        for i in range(0,len(selected)):
            for j in range(i+1, len(selected)):
                dice = uniform(0,100)
                if dice < self.crossoverRate:
                    father = selected[i]
                    mother = selected[j]
                    childOne = self.individualClass()
                    childTwo = self.individualClass()
                    children = self.crossoverMethod(father, mother, childOne, childTwo)
                    children[0] = self.applyMutationMethod(children[0])
                    children[1] = self.applyMutationMethod(children[1])
                    populationCross = populationCross + children
        #print("DEBUG -- population size: ", len(populationCross))
        return populationCross

    def applyMutationMethod(self, individual):
        dice = uniform(0,100)
        if dice < self.mutationRate:
            individual = self.mutationMethod(individual)
        return individual

    def getBest(self):
        max = 0
        best = None
        for individual in self.population:
            # print(individual.fitValue)
            if individual.fitValue > max:
                best = individual
                max = individual.fitValue
        return best

        

    def run(self, maxIterations=100):
        self.createInitialPopulation()
        for i in range (0,maxIterations):
            self.evaluatePopulation()
            selected = self.applySelectionMethod()
            populationCross = self.applyCrossoverMethod(selected)
            self.fillPopulation(populationCross)
            #self.fillPopulation(selected)
            avFitness = self.averageFitness()
            #print(avFitness)