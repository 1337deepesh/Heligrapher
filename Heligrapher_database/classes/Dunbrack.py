#This class performs grants easy access to Dunbrack's rotamer library:

import csv
import sys
import copy
import math
import numpy
import string
import random

class Dunbrack_library:
    
    ########################################
    #           INITIALISE CLASS           #
    ########################################
    def __init__(self, dunbrack_filename):
        #initialize variables:
        self.dunbrack_filename = dunbrack_filename
        self.dunbrack_dictionary = {}
        
        #read Dunbrack's database in .csv format:
        dunbrack_file = open(dunbrack_filename)
        dunbrack_csv = csv.reader(dunbrack_file, delimiter=',')
        dunbrack_data = []
        for row in dunbrack_csv:
            dunbrack_data.append(row)
        
        #convert data into a dictionary for easy access:
        #This dictionary stores residue name, phi and psi angles as a single string:
        for i in range(0, len(dunbrack_data)):
            if dunbrack_data[i][0] in self.dunbrack_dictionary:
                self.dunbrack_dictionary[str(dunbrack_data[i][0])].append({'probability': float(dunbrack_data[i][1]), 'chi':[float(j) for j in dunbrack_data[i][2:6]], 'SD':[float(j) for j in dunbrack_data[i][6:10]]})
            else:
                self.dunbrack_dictionary[str(dunbrack_data[i][0])] = [{'probability': float(dunbrack_data[i][1]), 'chi':[float(j) for j in dunbrack_data[i][2:6]], 'SD':[float(j) for j in dunbrack_data[i][6:10]]}]
    
    ########################################
    #   PRIVATE FUNCTION: weighted_choice  #
    ########################################
    def weighted_choice(self, weights):
        rnd = random.random() * sum(weights)
        for i, w in enumerate(weights):
            rnd -= w
            if rnd < 0:
                return i
    
    ########################################
    #        PUBLIC FUNCTION: Deliver      #
    ########################################
    #output all probability, chi, and SD vales for a given residue +phi/psi angles:
    def Deliver(self, residue_name, phi, psi):
        
        #ignore user call in cases of GLY/ALA. Simply send list of zeros:
        if residue_name == "GLY" or residue_name == "ALA":
            return [0, 0, 0, 0]
        
        #calculate closest angles to that of the input:
        phi_flag = float("inf")
        closest_phi = 0
        for i in range(-180, 181, 10):
            if abs(phi-i) < phi_flag:
                phi_flag = abs(phi-i)
                closest_phi = i
        
        psi_flag = float("inf")
        closest_psi = 0
        for i in range(-180, 181, 10):
            if abs(psi-i) < psi_flag:
                psi_flag = abs(psi-i)
                closest_psi = i
        
        #create string for dictionary search:
        search_string=residue_name+'% 4d' % closest_phi+'% 4d' % closest_psi
        
        #search database, output results as a list of sub-dictionaries:
        return self.dunbrack_dictionary[search_string]
    
    ########################################
    #         PUBLIC FUNCTION: Best        #
    ########################################
    #output chi values corresponding to highest probability rotamer:
    def Best(self, residue_name, phi, psi):
        
        #ignore user call in cases of GLY/ALA. Simply send list of zeros:
        if residue_name == "GLY" or residue_name == "ALA":
            return [0, 0, 0, 0]
        
        #call Deliver to output list of all chi angles:
        rotamers = self.Deliver(residue_name, phi, psi)
        
        #find highest probability rotamer:
        maximum_weight = 0
        maximum_address = 0
        for i in range(0, len(rotamers)):
            if rotamers[i]['probability'] > maximum_weight:
                maximum_weight = rotamers[i]['probability']
                maximum_address = i
        
        #output chi values:
        return rotamers[maximum_address]['chi']

    ########################################
    #     PUBLIC FUNCTION: Surprise_me     #
    ########################################
    #output randomly chosen chi vales for a given residue's phi/psi angles,
    #based on a probability function:
    def Surprise_me(self, residue_name, phi, psi):
        
        #ignore user call in cases of GLY/ALA. Simply send list of zeros:
        if residue_name == "GLY" or residue_name == "ALA":
            return [0, 0, 0, 0]
        
        #call Deliver to output list of all chi angles:
        rotamers = self.Deliver(residue_name, phi, psi)
        
        #make weights list for all chi/SD values:
        weights_list = []
        for i in range(0, len(rotamers)):
            weights_list.append(rotamers[i]['probability'])
        
        #choose a rotamer, output its chi values:
        return rotamers[self.weighted_choice(weights_list)]['chi']
    
    ########################################
    #     PUBLIC FUNCTION: Surprise_SD     #
    ########################################
    #output randomly chosen altered chi vales for a given residue's phi/psi angles,
    #based on a probability function. Chi values are randomly altered based on their SD:
    def Surprise_SD(self, residue_name, phi, psi):
        
        #ignore user call in cases of GLY/ALA. Simply send list of zeros:
        if residue_name == "GLY" or residue_name == "ALA":
            return [0, 0, 0, 0]
        
        #call Deliver to output list of all chi angles:
        rotamers = self.Deliver(residue_name, phi, psi)
        
        #make weights list for all chi/SD values:
        weights_list = []
        for i in range(0, len(rotamers)):
            weights_list.append(rotamers[i]['probability'])
        
        #store chi/SD values:
        chi= rotamers[self.weighted_choice(weights_list)]['chi']
        SD = rotamers[self.weighted_choice(weights_list)]['SD']
        #calculate new random chi values, chosen from gives gaussian parameters:
        chi_output = [0, 0, 0, 0]
        for i in range(0,len(chi_output)):
            #ensure that the output does not fall out of the +/-180 degree range:
            while(1):
                chi_output[i] = random.gauss(chi[i], SD[i])
                if chi_output[i] < -180 or chi_output[i] > 180:
                    continue
                else:
                    break

        #return altered chi values:
        return chi_output

















