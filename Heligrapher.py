#!/usr/bin/python
#Script started on: Tuesday 19 May 2015 06:01:32 PM IST 
#This script will generate graphs for all input antimicrobial helical structures.
#These graphs can then be used for subgraph comparison & scoring during SA-design.

import os
import csv
import sys
import copy
import math
import numpy
import shutil
import string
import random
import shutil
import getopt
import subprocess

from PIL import Image, ImageDraw, ImageFont

########################################
#       FUNCTION: GRAPH-PLOTTING       #
########################################

#create graph of score_log (recycle code from dinoplot):
def grapher(XY_data, OP_image, IP_database):
    
    #create a 'diamond' that will serve the standard graph-point:
    diamond_blue = [[255, 255, 255, 255, 255, 255],
                    [255, 255, 255, 255, 255, 255],
                    [255, 255, 255, 255, 255, 255],
                    [255, 255, 255, 255, 255, 255],
                    [255, 255, 255, 255, 255, 255],
                    [255, 255, 255, 255, 255, 255]]
    diamond_null = [[255, 255,   0,   0, 255, 255],
                    [255,   0,   0,   0,   0, 255],
                    [  0,   0,   0,   0,   0,   0],
                    [  0,   0,   0,   0,   0,   0],
                    [255,   0,   0,   0,   0, 255],
                    [255, 255,   0,   0, 255, 255]]
    
    diamond = [diamond_null, diamond_null, diamond_blue]
    
    W = 1000
    H = 1000
    
    img = Image.new("RGB", (W, H), "white")
    draw = ImageDraw.Draw(img)
    
    #create output image margin:
    for i in range(99,901):
        color = (0,0,0)
        draw.point((i,99), fill=color)
        draw.point((i,100), fill=color)
        draw.point((i,901), fill=color)
        draw.point((i,902), fill=color)
        draw.point((99,i), fill=color)
        draw.point((98,i), fill=color)
        draw.point((901,i), fill=color)
        draw.point((902,i), fill=color)

    #calculate axes range:
    X_min = float("inf")
    X_max = float("inf")*-1
    Y_min = float("inf")
    Y_max = float("inf")*-1

    for i in range(0, len(XY_data)):
        if XY_data[i][0] > X_max:
            X_max = float(XY_data[i][0])
        if XY_data[i][0] < X_min:
            X_min = float(XY_data[i][0])
        if XY_data[i][1] > Y_max:
            Y_max = float(XY_data[i][1])
        if XY_data[i][1] < Y_min:
            Y_min = float(XY_data[i][1])
    
    #calculate all points as fractions of min/max, convert to pixel location:
    #also plot on graph:
    for i in range(0, len(XY_data)):
        Y_point = (((XY_data[i][0]-X_min)/(X_max-X_min)*800)+100)
        X_point = 1000-(((XY_data[i][1]-Y_min)/(Y_max-Y_min)*800)+100)
        for j in range(0,len(diamond[0])):
            for k in range(0,len(diamond[0][0])):
                if diamond[0][j][k] == 0:
                    K = X_point+j-math.ceil(len(diamond[0])/2)
                    J = Y_point+k-math.ceil(len(diamond[0][0])/2)
                    color = (diamond[0][j][k], diamond[1][j][k], diamond[2][j][k])
                    draw.point((J,K), fill=color)

    #draw axes labels, ranges:
    font = ImageFont.truetype(IP_database+"/formatting/arial.ttf", 30)
    draw.text((445,950), "iteration", fill=(0,0,0), font=font)
    draw.text((45,420), "s", fill=(0,0,0), font=font)
    draw.text((45,450), "c", fill=(0,0,0), font=font)
    draw.text((45,480), "o", fill=(0,0,0), font=font)
    draw.text((45,510), "r", fill=(0,0,0), font=font)
    draw.text((45,540), "e", fill=(0,0,0), font=font)
    draw.text((310,50), "Simulated annealing (score log)", fill=(0,0,0), font=font)
    draw.text((100,930), str(X_min), fill=(128,128,128), font=font)
    draw.text((900,930), str(X_max), fill=(128,128,128), font=font)
    draw.text((10, 890), '% 4.2f' % Y_min, fill=(128,128,128), font=font)
    draw.text((10,90), '% 4.2f' % Y_max, fill=(128,128,128), font=font)

    #output image:
    img.save(OP_image, "png")

########################################
#       FUNCTION: CLASH-DETECTOR       #
########################################

#check for steric clashes between any 2 residues:
#PRE_function variables are generated for faster execution
#source: Seeliger, Daniel, and Bert L. de Groot. "Atomic contacts in protein structures. A detailed analysis of atomic radii, packing, and overlaps." Proteins: Structure, Function, and Bioinformatics 68.3 (2007): 595-601.
SF=1
R = {'H0':1.19*SF,   'HAR':1.14*SF,  'HA':1.03*SF,   'H':1.05*SF,
     'HC':0.58*SF,   'HDR':0.67*SF,  'C':1.43*SF,    'CA':1.48*SF,
     'CH1E':1.92*SF, 'CH2E':1.89*SF, 'CH3E':1.81*SF, 'CR1E':1.81*SF,
     'CR1W':1.76*SF, 'C5':1.76*SF,   'C5W':1.86*SF,  'CW':1.74*SF,
     'CH2G':1.76*SF, 'CH2P':1.47*SF, 'CY':1.87*SF,   'CY2':1.63*SF,
     'CF':1.83*SF,   'CDR':1.69*SF,  'CR1H':1.75*SF, 'CRHH':1.63*SF,
     'O':1.41*SF,    'OC':1.33*SF,   'OH1':1.31*SF,  'NH1':1.37*SF,
     'NH2':1.45*SF,  'NH3':1.35*SF,  'NC2':1.45*SF,  'NHS':1.40*SF,
     'SM':1.79*SF,   'S':1.83*SF}

#determine radii of all atoms of all amino acids:
GLY_radii = {' N  ':R['NH1'],  ' CA ':R['CH2G'],
             ' C  ':R['C'],    ' O  ':R['O'], 
             ' H  ':R['H'],    ' HA2':R['HA'],
             ' HA3':R['HA']}

ALA_radii = {' N  ':R['NH1'],  ' CA ':R['CA'],
             ' C  ':R['C'],    ' O  ':R['O'],
             ' CB ':R['CH3E'], ' H  ':R['H'],
             ' HA ':R['HA'],   ' HB1':R['H0'],
             ' HB2':R['H0'],   ' HB3':R['H0']}

VAL_radii = {' N  ':R['NH1'],  ' CA ':R['CA'], 
             ' C  ':R['C'],    ' O  ':R['O'], 
             ' CB ':R['CH1E'], ' CG1':R['CH3E'], 
             ' CG2':R['CH3E'], ' H  ':R['H'], 
             ' HA ':R['HA'],   ' HB ':R['H0'], 
             'HG11':R['H0'],   'HG12':R['H0'], 
             'HG13':R['H0'],   'HG21':R['H0'], 
             'HG22':R['H0'],   'HG23':R['H0']}

LEU_radii = {' N  ':R['NH1'],  ' CA ':R['CA'], 
             ' C  ':R['C'],    ' O  ':R['O'], 
             ' CB ':R['CH2E'], ' CG ':R['CH1E'], 
             ' CD1':R['CH3E'], ' CD2':R['CH3E'], 
             ' H  ':R['H'],    ' HA ':R['HA'], 
             ' HB2':R['H0'],   ' HB3':R['H0'], 
             ' HG ':R['H0'],   'HD11':R['H0'], 
             'HD12':R['H0'],   'HD13':R['H0'], 
             'HD21':R['H0'],   'HD22':R['H0'], 
             'HD23':R['H0']}

ILE_radii = {' N  ':R['NH1'],  ' CA ':R['CA'], 
             ' C  ':R['C'],    ' O  ':R['O'], 
             ' CB ':R['CH1E'], ' CG1':R['CH2E'], 
             ' CG2':R['CH3E'], ' CD1':R['CH3E'], 
             ' H  ':R['H'],    ' HA ':R['HA'], 
             ' HB ':R['H0'],   'HG12':R['H0'], 
             'HG13':R['H0'],   'HG21':R['H0'], 
             'HG22':R['H0'],   'HG23':R['H0'], 
             'HD11':R['H0'],   'HD12':R['H0'], 
             'HD13':R['H0']}

MET_radii = {' N  ':R['NH1'],  ' CA ':R['CA'], 
             ' C  ':R['C'],    ' O  ':R['O'], 
             ' CB ':R['CH2E'], ' CG ':R['CH2E'], 
             ' SD ':R['SM'],   ' CE ':R['CH3E'], 
             ' H  ':R['H'],    ' HA ':R['HA'],  
             ' HB2':R['H0'],   ' HB3':R['H0'], 
             ' HG2':R['H0'],   ' HG3':R['H0'], 
             ' HE1':R['H0'],   ' HE2':R['H0'], 
             ' HE3':R['H0']}

PHE_radii = {' N  ':R['NH1'],  ' CA ':R['CA'], 
             ' C  ':R['C'],    ' O  ':R['O'], 
             ' CB ':R['CH2E'], ' CG ':R['CF'], 
             ' CD1':R['CR1E'], ' CD2':R['CR1E'], 
             ' CE1':R['CR1E'], ' CE2':R['CR1E'], 
             ' CZ ':R['CR1E'], ' H  ':R['H'], 
             ' HA ':R['HA'],   ' HB2':R['H0'], 
             ' HB3':R['H0'],   ' HD1':R['HAR'], 
             ' HD2':R['HAR'],  ' HE1':R['HAR'], 
             ' HE2':R['HAR'],  ' HZ ':R['HAR']}

TYR_radii = {' N  ':R['NH1'],  ' CA ':R['CA'], 
             ' C  ':R['C'],    ' O  ':R['O'], 
             ' CB ':R['CH2E'], ' CG ':R['CY'],
             ' CD1':R['CR1E'], ' CD2':R['CR1E'], 
             ' CE1':R['CR1E'], ' CE2':R['CR1E'],
             ' CZ ':R['CY2'],  ' OH ':R['OH1'], 
             ' H  ':R['H'],    ' HA ':R['HA'], 
             ' HB2':R['H0'],   ' HB3':R['H0'], 
             ' HD1':R['HAR'],  ' HD2':R['HAR'], 
             ' HE1':R['HAR'],  ' HE2':R['HAR'], 
             ' HH ':R['H']}

TRP_radii = {' N  ':R['NH1'],  ' CA ':R['CA'],
             ' C  ':R['C'],    ' O  ':R['O'],
             ' CB ':R['CH2E'], ' CG ':R['C5W'], 
             ' CD1':R['CR1E'], ' CD2':R['CR1E'], 
             ' NE1':R['NH1'],  ' CE2':R['CR1W'], 
             ' CE3':R['CR1W'], ' CZ2':R['CR1W'], 
             ' CZ3':R['CR1W'], ' CH2':R['CR1E'], 
             ' H  ':R['H'],    ' HA ':R['HA'], 
             ' HB2':R['H0'],   ' HB3':R['H0'], 
             ' HD1':R['HAR'],  ' HE1':R['H'], 
             ' HE3':R['HAR'],  ' HZ2':R['HAR'], 
             ' HZ3':R['HAR'],  ' HH2':R['HAR']}

PRO_radii = {' N  ':R['NH1'],  ' CA ':R['CA'], 
             ' C  ':R['C'],    ' O  ':R['O'], 
             ' CB ':R['CH2P'], ' CG ':R['CH2P'], 
             ' CD ':R['CH2P'], ' HA ':R['HA'], 
             ' HB2':R['H0'],   ' HB3':R['H0'], 
             ' HG2':R['H0'],   ' HG3':R['H0'], 
             ' HD2':R['H0'],   ' HD3':R['H0']}

SER_radii = {' N  ':R['NH1'],  ' CA ':R['CA'], 
             ' C  ':R['C'],    ' O  ':R['O'], 
             ' CB ':R['CH2E'], ' OG ':R['OH1'], 
             ' H  ':R['H'],    ' HA ':R['HA'], 
             ' HB2':R['H0'],   ' HB3':R['H0'], 
             ' HG ':R['H']}

CYS_radii = {' N  ':R['NH1'],  ' CA ':R['CA'], 
             ' C  ':R['C'],    ' O  ':R['O'], 
             ' CB ':R['CH2E'], ' SG ':R['S'], 
             ' H  ':R['H'],    ' HA ':R['HA'], 
             ' HB2':R['H0'],   ' HB3':R['H0'], 
             ' HG ':R['H0']}

THR_radii = {' N  ':R['NH1'],  ' CA ':R['CA'], 
             ' C  ':R['C'],    ' O  ':R['O'], 
             ' CB ':R['CH1E'], ' OG1':R['OH1'], 
             ' CG2':R['CH3E'], ' H  ':R['H'], 
             ' HA ':R['HA'],   ' HB ':R['H0'], 
             ' HG1':R['H'],    'HG21':R['H0'], 
             'HG22':R['H0'],   'HG23':R['H0']}
    
ASN_radii = {' N  ':R['NH1'],  ' CA ':R['CA'], 
             ' C  ':R['C'],    ' O  ':R['O'], 
             ' CB ':R['CH2E'], ' CG ':R['C'], 
             ' OD1':R['O'],    ' ND2':R['NH2'], 
             ' H  ':R['H'],    ' HA ':R['HA'], 
             ' HB2':R['H0'],   ' HB3':R['H0'], 
             'HD21':R['H'],    'HD22':R['H']}

GLN_radii = {' N  ':R['NH1'],  ' CA ':R['CA'], 
             ' C  ':R['C'],    ' O  ':R['O'], 
             ' CB ':R['CH2E'], ' CG ':R['CH2E'], 
             ' CD ':R['C'],    ' OE1':R['O'], 
             ' NE2':R['NH2'],  ' H  ':R['H'], 
             ' HA ':R['HA'],   ' HB2':R['H0'], 
             ' HB3':R['H0'],   ' HG2':R['H0'], 
             ' HG3':R['H0'],   'HE21':R['H'], 
             'HE22':R['H']}

ASP_radii = {' N  ':R['NH1'],  ' CA ':R['CA'], 
             ' C  ':R['C'],    ' O  ':R['O'], 
             ' CB ':R['CH2E'], ' CG ':R['C'], 
             ' OD1':R['O'],    ' OD2':R['O'], 
             ' H  ':R['H'],    ' HA ':R['HA'], 
             ' HB2':R['H0'],   ' HB3':R['H0']}

GLU_radii = {' N  ':R['NH1'],  ' CA ':R['CA'], 
             ' C  ':R['C'],    ' O  ':R['O'], 
             ' CB ':R['CH2E'], ' CG ':R['CH2E'], 
             ' CD ':R['C'],    ' OE1':R['O'], 
             ' OE2':R['O'],    ' H  ':R['H'], 
             ' HA ':R['HA'],   ' HB2':R['H0'], 
             ' HB3':R['H0'],   ' HG2':R['H0'], 
             ' HG3':R['H0']}

LYS_radii = {' N  ':R['NH1'],  ' CA ':R['CA'], 
             ' C  ':R['C'],    ' O  ':R['O'],
             ' CB ':R['CH2E'], ' CG ':R['CH2E'], 
             ' CD ':R['CH2E'], ' CE ':R['CH2E'], 
             ' NZ ':R['NH3'],  ' H  ':R['H'], 
             ' HA ':R['HA'],   ' HB2':R['H0'], 
             ' HB3':R['H0'],   ' HG2':R['H0'], 
             ' HG3':R['H0'],   ' HD2':R['H0'], 
             ' HD3':R['H0'],   ' HE2':R['H0'], 
             ' HE3':R['H0'],   ' HZ1':R['HC'], 
             ' HZ2':R['HC'],   ' HZ3':R['HC']}

ARG_radii = {' N  ':R['NH1'],  ' CA ':R['CA'], 
             ' C  ':R['C'],    ' O  ':R['O'],
             ' CB ':R['CH2E'], ' CG ':R['CH2E'], 
             ' CD ':R['CDR'],  ' NE ':R['NC2'], 
             ' CZ ':R['CR1E'], ' NH1':R['NH2'], 
             ' NH2':R['NH2'],  ' H  ':R['H'], 
             ' HA ':R['HA'],   ' HB2':R['H0'], 
             ' HB3':R['H0'],   ' HG2':R['H0'], 
             ' HG3':R['H0'],   ' HD2':R['HDR'], 
             ' HD3':R['HDR'],  ' HE ':R['HC'], 
             'HH11':R['HC'],   'HH12':R['HC'], 
             'HH21':R['HC'],   'HH22':R['HC']}

HIS_radii = {' N  ':R['NH1'],  ' CA ':R['CA'], 
             ' C  ':R['C'],    ' O  ':R['O'],
             ' CB ':R['CH2E'], ' CG ':R['CR1E'], 
             ' ND1':R['NHS'],  ' CD2':R['CR1H'], 
             ' CE1':R['CRHH'], ' NE2':R['NHS'], 
             ' H  ':R['H'],    ' HA ':R['HA'], 
             ' HB2':R['H0'],   ' HB3':R['H0'], 
             ' HD2':R['H0'],   ' HE1':R['H0']}

PTN_radii = {'GLY':GLY_radii, 'ALA':ALA_radii, 'VAL':VAL_radii, 'LEU':LEU_radii, 
             'ILE':ILE_radii, 'MET':MET_radii, 'PHE':PHE_radii, 'TYR':TYR_radii, 
             'TRP':TRP_radii, 'PRO':PRO_radii, 'SER':SER_radii, 'CYS':CYS_radii, 
             'THR':THR_radii, 'ASN':ASN_radii, 'GLN':GLN_radii, 'ASP':ASP_radii, 
             'GLU':GLU_radii, 'LYS':LYS_radii, 'ARG':ARG_radii, 'HIS':HIS_radii}

#ligand radii are based on approximations/assumptions taken from protein radii:
#whenever a choice is presented, the smallest possible radius is used.
#it is assumed that AutoDock will penalise residues too close to the ligand:
LIG_radii = {'H':R['HC'], 'C':R['C'],  'N':R['NH3'],  'O':R['OH1'],
             'S':R['SM'], 'P':R['SM'], 'X':R['C']}

def clasher(residue_1, residue_2):
    
    #set atomic radii to a global variable.
    #This imports the radii given right before this function:
    global PTN_radii 
    global LIG_radii
    clash_flag = 0

    residue_list = ["ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"]
    ligand_list = ['H', 'C', 'N', 'O', 'S', 'P']

    #identify residue names:
    residue_1_name = residue_1[0][17:20]
    residue_2_name = residue_2[0][17:20]

    for i in range(0,len(residue_1)):
        if len(str(residue_1[i])) != 81:
            print "warning: <clasher> has encountered an incomplete residue (1)."
            print "now exiting <clasher> function..."
            return "WARNING"
        residue_1_coords = []
        residue_1_coords.append(float(residue_1[i][30:38]))
        residue_1_coords.append(float(residue_1[i][38:46]))
        residue_1_coords.append(float(residue_1[i][46:54]))
        #determine whether selected atom (i) belongs to backbone:
        residue_1_atom = copy.copy(residue_1[i][12:16])
        residue_1_backbone_flag = 0
        if residue_1_atom == " N  " or residue_1_atom == " CA " or residue_1_atom == " C  " or residue_1_atom == " O  ":
            residue_1_backbone_flag = 1
        
        for j in range(0,len(residue_2)):
            if len(str(residue_2[j])) != 81:
                print "warning: <clasher> has encountered an incomplete residue (2)."
                print "now exiting <clasher> function..."
                return "WARNING"
            #determine whether selected atom (j) belongs to backbone:
            residue_2_atom = copy.copy(residue_2[j][12:16])
            residue_2_backbone_flag = 0
            if residue_2_atom == " N  " or residue_2_atom == " CA " or residue_2_atom == " C  " or residue_2_atom == " O  ":
                residue_2_backbone_flag = 1
            #IGNORE backbone-backbone clashes:
            if residue_1_backbone_flag != 1 or residue_2_backbone_flag != 1:  
                #calculate the minimum acceptable inter-atomic radius:
                #read atomic X,Y,Z co-ordinates:
                residue_2_coords = []
                residue_2_coords.append(float(residue_2[j][30:38]))
                residue_2_coords.append(float(residue_2[j][38:46]))
                residue_2_coords.append(float(residue_2[j][46:54]))
                #calculate interatomic distance:
                distance = math.sqrt(pow(residue_1_coords[0]-residue_2_coords[0],2) +pow(residue_1_coords[1]-residue_2_coords[1],2) +pow(residue_1_coords[2]-residue_2_coords[2],2))

                residue_2_atom = residue_2[j][12:16]
                #extract atomic radii from double dictionary (radii):
                #check if inputted residue is protein or ligand:
                if residue_1_name in residue_list:
                    #identify atom types:
                    residue_1_atom = residue_1[i][12:16]
                    residue_1_radius = PTN_radii[residue_1_name][residue_1_atom]
                else:
                    residue_1_atom = residue_1[i][77]
                    if residue_1_atom not in ligand_list:
                        residue_1_atom = 'X'
                    residue_1_radius = LIG_radii[residue_1_atom]
                if residue_2_name in residue_list:
                    #identify atom types:
                    residue_2_atom = residue_2[j][12:16]
                    residue_2_radius = PTN_radii[residue_2_name][residue_2_atom]
                else:
                    residue_2_atom = residue_2[j][77]
                    if residue_2_atom not in ligand_list:
                        residue_2_atom = 'X'
                    residue_2_radius = LIG_radii[residue_2_atom]
                #calculate minimum permissible interatomic distance:
                min_distance = residue_1_radius +residue_2_radius
                if distance <= min_distance:
                    clash_flag = 1
    #print "CLASH: "+str(clash_flag)+" | "+residue_1_name+":"+residue_1_atom+"| "+residue_2_name+":"+residue_2_atom
    return clash_flag

########################################
#        FUNCTION: JUGAAD-FIXER        #
########################################

#check the class-based PDB-structure for errors (jugaad fix): 
def jugaad_fixer(PDBhelix):
    
    #scan throught 'PDBhelix' to spot any lines that contain only '0':
    for residue in range(0,len(PDBhelix)):
        error_flag = False
        for line in range(0, len(PDBhelix[residue].IP_residue)):
            if len(str(PDBhelix[residue].IP_residue[line])) != 81:
                error_flag = True
                break
        #correct a residue ONLY if a '0' line is spotted:
        if error_flag == True:
            recovery_list = []
            for line in range(0, len(PDBhelix[residue].IP_residue)):
                if len(str(PDBhelix[residue].IP_residue[line])) == 81:
                    recovery_list.append(PDBhelix[residue].IP_residue[line])
            PDBhelix[residue].IP_residue = copy.deepcopy(recovery_list)

########################################
#       FUNCTION: GRAPH-CREATION       #
########################################

def graph_maker(PDB_input, output_directory, database, verbose_flag):
    
    #extract names & coordinates for all PDB data:
    protein = extract_PDBdata(PDB_input)
    
    #run DSSP, determine helicity & hydrogen bonding pattern:
    subprocess.call([database+"/third_party_software/dssp-2.0.4-linux-amd64", "-i", PDB_input, "-o", output_directory+"/DSSP_output.txt"])

    #list all helical residues in DSSP file:
    DSSP_file = open(output_directory+"/DSSP_output.txt")
    DSSP_data = DSSP_file.readlines()
    DSSP_graph = {}
    helical_residues = []
        
    start_flag = False
    for i in range(0,len(DSSP_data)):
        if DSSP_data[i][0:15] == "  #  RESIDUE AA":
            start_flag = True
            continue
        if start_flag == True and DSSP_data[i][16] == 'H':
            helical_residues.append(int(DSSP_data[i][6:10]))

    #create DSSP-based helical graph as adjacency list:
    start_flag = False
    for i in range(0,len(DSSP_data)):
        if DSSP_data[i][0:15] == "  #  RESIDUE AA":
            start_flag = True
            continue
        if start_flag == True and DSSP_data[i][16] == 'H':
            residue = int(DSSP_data[i][6:10])
            DSSP_graph[residue] = {}
            
            #create 'NH->O' edge:
            DSSP_graph[residue]['NH->O'] = None
            if DSSP_data[i][43:45] == '-4':
                partner_residue_index = helical_residues.index(residue)-4
                if partner_residue_index < 0:
                    if verbose_flag == True:
                        print "    note: DSSP-predicted H-bonding partner not found."
                        print "    ignoring -4 'NH->O' for residue "+str(residue)+" and continuing search."
                        print "    "
                else:
                    partner_residue = helical_residues[partner_residue_index]
                    DSSP_graph[residue]['NH->O'] = partner_residue

            #create 'O->NH' edge:
            DSSP_graph[residue]['O->NH'] = None
            if DSSP_data[i][54:56] == ' 4':
                partner_residue_index = helical_residues.index(residue)+4
                if partner_residue_index > len(helical_residues)-1:
                    if verbose_flag == True:
                        print "    note: DSSP-predicted H-bonding partner not found."
                        print "    ignoring 4 'O->NH' for residue "+str(residue)+" and continuing search."
                        print "    "
                else:
                    partner_residue = helical_residues[partner_residue_index]
                    DSSP_graph[residue]['O->NH'] = partner_residue

            #create 'before' edge:
            DSSP_graph[residue]['before'] = None
            partner_residue_index = helical_residues.index(residue)-1
            if partner_residue_index < 0:
                if verbose_flag == True:
                    print "    note: helix 'before' neighbour not found."
                    print "    ignoring neighbour for residue "+str(residue)+" and continuing search."
                    print "    "
            else:
                partner_residue = helical_residues[partner_residue_index]
                DSSP_graph[residue]['before'] = partner_residue
            
            #create 'after' edge:
            DSSP_graph[residue]['after'] = None
            partner_residue_index = helical_residues.index(residue)+1
            if partner_residue_index > len(helical_residues)-1:
                if verbose_flag == True:
                    print "    note: helix 'after' neighbour not found."
                    print "    ignoring neighbour for residue "+str(residue)+" and continuing search."
                    print "    "
            else:
                partner_residue = helical_residues[partner_residue_index]
                DSSP_graph[residue]['after'] = partner_residue

            #mention residue name:
            DSSP_graph[residue]['name'] = DSSP_data[i][13]
            
            #add residue label:
            DSSP_graph[residue]['label'] = False

    return DSSP_graph

########################################
#    FUNCTION: FIND MAXCON SUBGRAPH    #
########################################

def maxcon_subgraphs(DB_graph, IP_graph, bijection_threshold):
    
    DB_keys = DB_graph.keys()
    IP_keys = IP_graph.keys()
    
    #solve MCS problem by brute force. Crude but effective.
    counter = 0
    all_node_pairs = []
    consolation_node_pairs = []
    for i in range(0,len(DB_keys)):
        for j in range(0,len(IP_keys)):
            node_pairs = []
            DFS(DB_graph, DB_keys[i], IP_graph, IP_keys[j], node_pairs)
            #add node-pairs to main list (all_node_pairs) only if length exceeds threshold:
            if len(node_pairs) >= bijection_threshold:
                for k in range(0,len(node_pairs)):
                    if node_pairs[k] not in all_node_pairs:
                        all_node_pairs.append(node_pairs[k])
            #'consolation_node_pairs' are scores to help SA algo move in right graphical direction:
            elif len(node_pairs) >= 2:
                for k in range(0,len(node_pairs)):
                    if node_pairs[k] not in consolation_node_pairs:
                        consolation_node_pairs.append(node_pairs[k])
    
    #extract unique database nodes from node-pairs.
    node_list = []
    for i in range(0,len(all_node_pairs)):
        if all_node_pairs[i][0] not in node_list:
            node_list.append(all_node_pairs[i][0])

    return len(node_list)+len(consolation_node_pairs)*0.001

########################################
#   FUNCTION: OUTPUT MAXCON SUBGRAPHS  #
########################################

def OP_subgraphs(DB_graph, IP_graph, bijection_threshold, database):
    
    DB_keys = DB_graph.keys()
    IP_keys = IP_graph.keys()
    
    #solve MCS problem by brute force. Crude but effective.
    counter = 0
    all_node_pairs = []
    OP_node_pairs = []
    
    for i in range(0,len(DB_keys)):
        for j in range(0,len(IP_keys)):
            node_pairs = []
            DFS(DB_graph, DB_keys[i], IP_graph, IP_keys[j], node_pairs)
            #add node-pairs to main list (all_node_pairs) only if length exceeds threshold:
            if len(node_pairs) >= bijection_threshold:
                sorted_node_pairs = sorted(node_pairs)
                if sorted_node_pairs not in all_node_pairs:
                    interim_node_pairs = []
                    for residue in range(0,len(sorted_node_pairs)):
                        interim_node_pairs.append([  sorted_node_pairs[residue][0], DB_graph[sorted_node_pairs[residue][0]]['name'], sorted_node_pairs[residue][1], IP_graph[sorted_node_pairs[residue][1]]['name']  ])
                    all_node_pairs.append(sorted_node_pairs)
                    OP_node_pairs.append(interim_node_pairs)
    
    return OP_node_pairs
    
########################################
#     FUNCTION: DEPTH-FIRST SEARCH     #
########################################

def DFS(graph1, graph1_node, graph2, graph2_node, node_pairs):
    
    #perform actual depth-first-search:
    DFS_algorithm(graph1, graph1_node, graph2, graph2_node, node_pairs)
    #cleanup 'labels' after depth-first-search:
    #also count number of 'True's:
    graph1_keys = graph1.keys()
    for i in range(0, len(graph1_keys)):
        graph1[graph1_keys[i]]['label'] = False
    
    #END function
    return None

def DFS_algorithm(graph1, graph1_node, graph2, graph2_node, node_pairs):
    
    #             UP       DOWN     LEFT      RIGHT:
    direction = ['NH->O', 'O->NH', 'before', 'after']
    
    #IF graph1_node & graph2_node are NOT blank:
    if graph1_node != None and graph2_node != None:
        #IF graph1_node has NOT been visited before:
        if graph1[graph1_node]['label'] != True:
            #IF graph1_node & graph2_node share the same residue name:
            if graph1[graph1_node]['name'] == graph2[graph2_node]['name']:
                #label graph1_node:
                graph1[graph1_node]['label'] = True
                #add graph1/graph2 pair to node_pairs for later tallying:
                node_pairs.append([graph1_node, graph2_node])
                #recurse to all neighbouring nodes:
                for i in range(0,len(direction)):
                    DFS_algorithm(graph1, graph1[graph1_node][direction[i]], graph2, graph2[graph2_node][direction[i]], node_pairs)
    
    #END function
    return None

########################################
#      FUNCTION: SimAnneal SCORE       #
########################################
def SA_scorer(IP_graph, graph_list, bijection_threshold):
    SA_score = 0
    for i in range(0,len(graph_list)):
        SA_score += maxcon_subgraphs(IP_graph, graph_list[i], bijection_threshold)
    
    return SA_score*-1
    
########################################
#     FUNCTION: PDB-DATA EXTRACTION    #
########################################

def extract_PDBdata(PDB_input):
    
    PDB_file = open(PDB_input)
    PDB_data = PDB_file.readlines()
    old_residue_ID = PDB_data[0][21:27]
    new_residue_ID = []
    residue_store = []
    protein = {}
    #fill 'protein' dictionary with all residue numbers:
    for i in range(0, len(PDB_data)):
        if PDB_data[i][0:6] == "ATOM  ":
            protein[int(PDB_data[i][22:26])] = {}
    #fill 'protein' dictionary with all details for respective residue numbers:
    for i in range(0, len(PDB_data)):
        if PDB_data[i][0:6] == "ATOM  " and (PDB_data[i][12:16] == " N  " or PDB_data[i][12:16] == " CA " or PDB_data[i][12:16] == " C  " or PDB_data[i][12:16] == " O  "):
            #3-letter residue name:
            protein[int(PDB_data[i][22:26])]['name'] = PDB_data[i][17:20]
            #atom type:
            protein[int(PDB_data[i][22:26])]['atom'] = PDB_data[i][12:16]
            #3D coordinates:
            protein[int(PDB_data[i][22:26])]['XYZ'] = [float(PDB_data[i][30:38]), float(PDB_data[i][38:46]), float(PDB_data[i][46:54])]
    
    return protein

########################################
#      FUNCTION: COOLING SCHEDULE      #
########################################

#Select cooling schedule based on user input:
def cooling_scheduler(i, cycles, cooling_schedule, start_temperature, end_temperature):
    i = float(i)
    T0 = float(start_temperature)
    Tn = float(end_temperature)
    N = float(cycles)
    
    if cooling_schedule == 0:
        Ti = T0 -i*(T0-Tn)/N
    
    if cooling_schedule == 1:
        Ti = T0*(Tn/T0)**(i/N)
    
    if cooling_schedule == 2:
        A = (T0-Tn)*float(N+1)/N
        B = T0 -A
        Ti = A/(i+1) +B
    
    if cooling_schedule == 3:
        print "warning: cooling_schedule '3' does not work as described"
        print "switching to cooling_schedule '5'"
        cooling_schedule = 5
    
    if cooling_schedule == 4:
        Ti = (T0-Tn)/(1+math.exp(0.01*(i-N/2))) +Tn;
    
    if cooling_schedule == 5:
        Ti = 0.5*(T0 -Tn)*(1+math.cos(i*math.pi/N)) +Tn
    
    if cooling_schedule == 6:
        Ti = 0.5*(T0-Tn)*(1-math.tanh(i*10/N-5)) +Tn;
    
    if cooling_schedule == 7:
        Ti = (T0-Tn)/math.cosh(i*10/N) +Tn;
    
    if cooling_schedule == 8:
        A = (1/N)*math.log(T0/Tn)
        Ti = T0*math.exp(-A*i)

    if cooling_schedule == 9:
        A = (1/N**2)*math.log(T0/Tn)
        Ti = T0*math.exp(-A*i**2);
    
    return Ti

########################################
#     FUNCTION: ACCEPTANCE_DECISION    #
########################################
def acceptance_decision(current_score, previous_score, temperature, probability_flag):
    if current_score < previous_score:
        return "accept"
    else:
        if previous_score == 0:
            previous_score = -0.001
        if current_score == 0:
            current_score = -0.001
        score_ratio = current_score/previous_score
        if probability_flag == True:
            acceptance_threshold = temperature*(  1/(1+1000*math.exp(-12*score_ratio))  )
        else:
            acceptance_threshold = temperature
        #acceptance_threshold = temperature
        if random.random() < acceptance_threshold:
            return "accept"
        else:
            return "reject"
     
########################################
#   FUNCTION: CHECK IF PROGRAM EXISTS  #
########################################

def program_check(program):
    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            path = path.strip('"')
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file

    return None


########################################
#    FUNCTION: SIMULATED ANNEALING     #
########################################
def simulated_anneal(IP_graph, graph_list, cycles, cooling_schedule, bijection_threshold, probability_flag, design_data, verbose_flag):

    #list all amino-acid residues:
    residue_list = ['G', 'A', 'V', 'L', 'I', 'M', 'F', 'Y', 'W', 'P', 'S', 'C', 'T', 'N', 'Q', 'D', 'E', 'K', 'R', 'H']
    OLD_score = float("inf")
    BEST_score = OLD_score
    BEST_sequence = []
    score_log = {}
    score_log["iteration"] = []
    score_log["temperature"] = []
    score_log["score"] = []
    score_log["sequence"] =[]
    
    for i in range(0,cycles):
        
        #I'm using cooling-schedules from: http://www.btluke.com/simanf1.html
        #temperature ranges from 1:0.01
        temperature = cooling_scheduler(i, cycles, cooling_schedule, 1, 1/cycles)
        
        #progress reports and progress logging:
        if verbose_flag == True:
            print "    iteration: "+str(i)+" score: "+str(OLD_score)+" temperature: "+str(temperature)

        IP_graph_keys = IP_graph.keys()
        current_sequence = []
        for j in range(0,len(IP_graph_keys)):
            current_sequence.append(IP_graph[IP_graph_keys[j]]['name'])
        current_sequence = ''.join(current_sequence)
        score_log["iteration"].append(i)
        score_log["temperature"].append(temperature)
        score_log["score"].append(OLD_score)
        score_log["sequence"].append(current_sequence)
        
        #mutate a randomly-selected residue:
        chosen_residue = random.choice(IP_graph_keys)
        residue_storage = IP_graph[chosen_residue]['name']
        IP_graph[chosen_residue]['name'] = random.choice(design_data[chosen_residue])

        #score mutant:
        NEW_score = SA_scorer(IP_graph, graph_list, bijection_threshold)

        #determine whether to perform mutation:
        decision = acceptance_decision(NEW_score, OLD_score, temperature, probability_flag)
        if decision == "accept":
            #update OLD_score:
            OLD_score = NEW_score
            #update BEST_score if necessary. also store corresponding sequence:
            if NEW_score < BEST_score:
                BEST_score = NEW_score
                BEST_sequence = []
                BEST_graph = {}
                for j in range(0,len(IP_graph_keys)):
                    BEST_sequence.append(IP_graph[IP_graph_keys[j]]['name'])
                BEST_sequence = ''.join(BEST_sequence)
                BEST_graph = copy.deepcopy(IP_graph)
        if decision == "reject":
            #revert to previous residue:
            IP_graph[chosen_residue]['name'] = residue_storage
    
    #END function:
    return [BEST_graph, BEST_sequence, BEST_score, score_log]

########################################
#       DATA INPUT AND PROCESSING      #
########################################

#standard usage/error message:
if len(sys.argv) < 7 or len(sys.argv) > 23:
    print "usage: Heligrapher.py -i <Input_data.txt> -o <output_directory> -d <database>"
    print "MANDATORY FLAGS:"
    print "  -i/--input: path to Input_data.txt (helical parameters)"
    print "  -o/--output: choose the path for all output files"
    print "  -d/--database: path to Heligrapher database."
    print "OPTIONAL FLAGS:"
    print "  -n/--iterations: number of simulated annealing iterations (default=1000)."
    print "  -a/--annealer: choose an annealing schedule (0 to 9, default 7)"
    print "  -s/--seed: define a random-number seed for reproducible runs."
    print "  -b/--bijection: choose a bijection threshold for residue-subgraph comparison (default=5)."
    print "  -p/--probability: enable acceptance probability function (Y/N, default N)."
    print "  -v/--verbose: print progesss/warnings on terminal (Y/N, default Y)"
    print "  -S/--structure: output a PDB structure of the best designed sequence (Y/N, default N)"
    print "  -P/--pymol: output a .png image using pymol that colours maximally intersected"
    print "              residues green, minimally intersected residues red (Y/N, default N)"
    sys.exit()

#definitions:
opts, args = getopt.getopt(sys.argv[1:], "i:o:d:n:a:s:b:p:v:S:P:", ["input=", "output=", "database=", "iterations=", "annealer=", "seed=", "bijection=", "probability=","verbose=", "structure=", "pymol="])

#read all input arguments:
cycles = 1000
cooling_schedule = 7
seed = random.randint(0, sys.maxint)
bijection_threshold = 5
verbose_flag = True
probability_flag = False
structure_flag = False
pymol_flag = False

for o,a in opts:
    if o == "-i" or o == "--input":
        input_data = a
    if o == "-o" or o == "--output":
        output_directory = a
    if o == "-d" or o == "--database":
        database = a
    if o == "-n" or o == "--iterations":
        cycles = int(a)
    if o == "-a" or o == "--annealer":
        cooling_schedule = int(a)
    if o == "-s" or o == "--seed":
        seed = a
    if o == "-b" or o == "--bijection":
        bijection_threshold = int(a)
    if o == "-p" or o == "--probability":
        probability_flag = a
    if o == "-v" or o == "--verbose":
        verbose_flag = a
    if o == "-S" or o == "--structure":
        structure_flag = a
    if o == "-P" or o == "--pymol":
        pymol_flag = a

if verbose_flag == 'Y':
    verbose_flag = True
if verbose_flag == 'N':
    verbose_flag = False

if probability_flag == 'Y':
    probability_flag = True
if probability_flag == 'N':
    probability_flag = False

if structure_flag == 'Y':
    structure_flag = True
if structure_flag == 'N':
    structure_flag = False

if pymol_flag == 'Y':
    pymol_flag = True
if pymol_flag == 'N':
    pymol_flag = False

#Use the random-number seed to initialize random-number generation:
random.seed(seed)

#make output_directory (if absent):
if os.path.isdir(output_directory) == False:
    os.makedirs(output_directory)

#check contents of output_directory. Terminate script if not completely empty:
if len(os.listdir(output_directory)) != 0:
    print "error: --output directory is not empty."
    print "Heligrapher.py will now terminate..."
    sys.exit()

#read data on the length/characteristics of desired antimicrobial helix:
Helix_length = 0
design_data = {}
with open(input_data, 'rb') as csvfile:
    line_reader = csv.reader(csvfile, delimiter=',')
    for row in line_reader:
        if "HELIX_LENGTH" in row[0]:
            Helix_length = int(row[1])
        if "PIKAA" in row[0]:
            design_data[int(row[1])] = row[2].replace(" ", "")

if Helix_length < 8 or Helix_length > 40:
    print "error: HELIX_LENGTH must be between 8 and 40 residues"
    sys.exit()

#ensure that all input residue numbers fall within the helical length:
design_data_keys = design_data.keys()
for i in range(0,len(design_data_keys)):
    if design_data_keys[i] <= 0:
        print "error: residue "+str(design_data_keys[i])+" should be >= 1"
        sys.exit()
    if design_data_keys[i] > Helix_length:
        print "error: residue "+str(design_data_keys[i])+" should be <= helix_length("+str(Helix_length)+")"
        sys.exit()

#fill in any missing keys. Use all residues:
for residue in range(1,Helix_length+1):
    if residue not in design_data_keys:
        design_data[residue] = "GAVLIMFYTWSCTNQDEKRH"

#create an input graph for the desired antimicrobial helix:
IP_graph = {}
for residue in range(1,Helix_length+1):
    IP_graph[residue] = {}
    IP_graph[residue]['name'] = 'X'
    IP_graph[residue]['label'] = False
    if residue >= 2:
        IP_graph[residue]['before'] = residue-1
    else:
        IP_graph[residue]['before'] = None
    if residue < Helix_length:
        IP_graph[residue]['after'] = residue+1
    else:
        IP_graph[residue]['after'] = None
    if residue >= 5:
        IP_graph[residue]['NH->O'] = residue-4
    else:
        IP_graph[residue]['NH->O'] = None
    if residue <= Helix_length-4:
        IP_graph[residue]['O->NH'] = residue+4
    else:
        IP_graph[residue]['O->NH'] = None

#create graphs for database PDB files:
graph_list = []
AM_helix_path = database+"/antimicrobial_helices/"
AM_helices = [ f for f in os.listdir(AM_helix_path) if os.path.isfile(os.path.join(AM_helix_path,f)) ]
for i in range(0,len(AM_helices)):
    graph_list.append(graph_maker(AM_helix_path+AM_helices[i], output_directory, database, verbose_flag))

########################################
#   RUN SIMULATED ANNEALING PROTOCOL   #
########################################

if verbose_flag == True:
    print "initiating annealing protocol..."
    print ""
    
[BEST_graph, BEST_sequence, BEST_score, score_log] = simulated_anneal(IP_graph, graph_list, cycles, cooling_schedule, bijection_threshold, probability_flag, design_data, verbose_flag)

if verbose_flag == True:
    print ""
    print "designed helix sequence: "+BEST_sequence+" | score: "+str(BEST_score)
    print ""

########################################
#              DATA OUTPUT             #
########################################

#output BEST_sequence.fasta:
BEST_sequence_file = open(str(output_directory)+'/BEST_sequence.fasta', 'w+')
BEST_sequence_file.write(">designed_helix\n")
BEST_sequence_file.write("%s\n" % BEST_sequence)
BEST_sequence_file.close()

#output score_log.csv:
score_log_file = open(str(output_directory)+'/score_log.txt', 'w+')
score_log_keys = score_log.keys()
score_log_file.write(str(score_log_keys[0])+", "+str(score_log_keys[1])+", "+str(score_log_keys[2])+", "+str(score_log_keys[3])+"\n")
for i in range(0, len(score_log[score_log_keys[0]])):
    iteration = score_log[score_log_keys[0]][i]
    temperature = score_log[score_log_keys[1]][i]
    score = score_log[score_log_keys[2]][i]
    sequence = score_log[score_log_keys[3]][i]
    score_log_file.write(str(iteration)+", "+str(temperature)+", "+str(score)+", "+str(sequence)+"\n")

score_log_file.close()

#output score_log.png. format XY_data:
XY_data = []
for i in range(0, len(score_log[score_log_keys[0]])):
    if score_log['score'][i] < float('inf'):
        XY_data.append([score_log['iteration'][i], score_log['score'][i]])

grapher(XY_data, str(output_directory)+"/score_log.png", str(database))

#output BEST_sequence.pdb:
#check if the structure has fully annealed:
if structure_flag == True and 'X' in BEST_sequence:
    
    print "    error: BEST_sequence not fully annealed (contains 'X' residues)."
    print "    the 'BEST_sequence.pdb' file cannot be written."
    
elif structure_flag == True:
    
    #Import class 'PDB_residue'. Useful for rotamer generation.
    execfile(database+'/classes/PDB_residue.py')
    
    #Import class 'Dunbrack'. Useful for easy access to Dunbrack's rotamer library.
    execfile(database+'/classes/Dunbrack.py')
    
    #use a dictionary for quick access to naive residues based on residue name:
    naive_residue = {}
    residue_list = ["ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"]
    
    for i in range(0, len(residue_list)):
        naive_file = open(database+"/PDBnaive_files/"+residue_list[i]+"_naiveH.pdb")
        naive_residue[residue_list[i]] = naive_file.readlines()
    
    #create a PDB-format helix of desired length:
    PDBhelix_file = open(database+"/ideal_helix/helix_40.pdb")
    PDBhelix_data = PDBhelix_file.readlines()
    old_residue_ID = PDBhelix_data[0][21:27]
    new_residue_ID = []
    residue_store = []
    PDBhelix = []
    residue_counter = 0
    for i in range(0, len(PDBhelix_data)):
        if PDBhelix_data[i][0:6] == "ATOM  ":
            #extract individual amino acid residues, store as 2D lists:
            new_residue_ID = copy.copy(PDBhelix_data[i][21:27])
            if new_residue_ID != old_residue_ID or i == len(PDBhelix_data)-1:
                if i == len(PDBhelix_data)-1:
                    residue_store.append(copy.copy(PDBhelix_data[i]))
                PDBhelix.append(PDB_residue(residue_store, naive_residue[residue_store[0][17:20]]))
                residue_store = []
                residue_counter += 1
                #stop inputting residues when the user-defined 'Helix_length' has been reached:
                if residue_counter == Helix_length:
                    break
            residue_store.append(copy.copy(PDBhelix_data[i]))
            old_residue_ID = new_residue_ID
    
    #extract Dunbrack's rotamer library (2010). Store in a dictionary-based class:
    Dunbrack = Dunbrack_library(database+"/dunbrack/dunbrack_rotamers_2010.csv")
    
    #output structure. Mutate residues (surprise_SD) on IP poly-alanine helix to best sequence (from simulated annealing):
    residue_dictionary = {'G':'GLY', 'A':'ALA', 'V':'VAL', 'L':'LEU', 'I':'ILE', 'M':'MET', 'F':'PHE', 'Y':'TYR', 'W':'TRP', 'P': 'PRO', 'S':'SER', 'C':'CYS', 'T':'THR', 'N':'ASN', 'Q':'GLN', 'D':'ASP', 'E':'GLU', 'K':'LYS', 'R':'ARG', 'H':'HIS'}
    
    #mutate the previously-created class-based PDB-structure:
    for residue in range(0, len(PDBhelix)):
        [PDBhelix[residue].phi, PDBhelix[residue].psi] = [-65.0, -35.0]
        mutant = residue_dictionary[BEST_sequence[residue]]
        chi = Dunbrack.Surprise_SD(mutant, PDBhelix[residue].phi, PDBhelix[residue].psi)
        PDBhelix[residue].Alter(naive_residue[mutant], PDBhelix[residue].IP_residue, chi)

    #check the class-based PDB-structure for errors (jugaad fix): 
    jugaad_fixer(PDBhelix)
    
    #check class-based PDB-structure for clashes:
    resolution_attempts = 20
    for i in range(0,resolution_attempts):
        if verbose_flag == True:
            print "    clash-resolving attempt "+str(i+1)+"/"+str(resolution_attempts)
        clash_flag = 0
        for res1 in range(0,len(PDBhelix)):
            for res2 in range(0,len(PDBhelix)):
                #if a steric clash is detected:
                if res1 != res2 and clasher(PDBhelix[res1].IP_residue, PDBhelix[res2].IP_residue)==1:
                    if verbose_flag == True:
                        print "    note: steric clash detected between residues "+str(res1)+" and "+str(res2)
                    if random.random() < 0.5:
                        resX = res1
                    else:
                        resX = res2
                    #correct a randomly-chosen clashing residue:
                    mutant = residue_dictionary[BEST_sequence[resX]]
                    chi = Dunbrack.Surprise_SD(mutant, PDBhelix[resX].phi, PDBhelix[resX].psi)
                    PDBhelix[resX].Alter(naive_residue[mutant], PDBhelix[resX].IP_residue, chi)
                    #perform 'jugaad_fix' to remove lines consisting only of '0':
                    jugaad_fixer(PDBhelix)
                    clash_flag += clash_flag
        
        if clash_flag == 0:
            if verbose_flag == True:
                print "    all clashes have been resolved"
            break

    if clash_flag != 0 and verbose_flag == True:
        print "    warning: all clashes have NOT been resolved" 
    
    #output class-based PDB structure (external to python):
    BEST_helix = open(str(output_directory)+'/BEST_helix.pdb', 'w+')
    for residue in range(0, len(PDBhelix)):
        for line in range(0, len(PDBhelix[residue].IP_residue)):
            BEST_helix.write("%s" % PDBhelix[residue].IP_residue[line])
    
    BEST_helix.close()

#perform graph-searches to determine what graphs the BEST_helix draws data from:
subgraph_match_counter = [0]*Helix_length
graphs_used = []
subgraphs_used = []
for i in range(0,len(graph_list)):
    subgraph_matches = OP_subgraphs(graph_list[i], BEST_graph, bijection_threshold, database)
    if len(subgraph_matches) != 0:
        graphs_used.append(AM_helices[i])
        subgraph_string = AM_helices[i]
        for node in subgraph_matches[0]:
            subgraph_string = subgraph_string+' '+str(node[0])+'/'+str(node[2])+'/'+node[1]
        subgraphs_used.append(subgraph_string)
        for node in subgraph_matches[0]:
            subgraph_match_counter[node[2]-1] += 1

#output common_graphs.txt:
common_graph = open(str(output_directory)+'/common_graph.txt', 'w+')

common_graph.write("graph overlaps per residue:\n\n")
for i in range(0,Helix_length):
    common_graph.write("%s %d\n" % (BEST_sequence[i], subgraph_match_counter[i]))
    
common_graph.write("\n")

common_graph.write("all subgraphs are listed below in this format:\n")
common_graph.write("database_node_number/designed_node_number/residue_type\n\n")
for line in subgraphs_used:
    common_graph.write(line+"\n")
common_graph.write("\n")

common_graph.write("%d database peptide-graphs were used to create BEST_helix:\n\n" % len(graphs_used))
print_counter = 0
while True:
    for j in range(0,4):
        if print_counter == len(graphs_used):
            break
        common_graph.write("%s " % graphs_used[print_counter])
        print_counter += 1
    common_graph.write("\n")
    if print_counter == len(graphs_used):
        break
        
common_graph.close()  

#output common_graphs.pymol:
common_graph_pymol = open(str(output_directory)+'/common_graph.pymol', 'w+')
common_graph_pymol.write("show cartoon\n")
common_graph_pymol.write("show sticks\n")
common_graph_pymol.write("hide (h.)\n")
for i in range(0,Helix_length):
    greenery = float(subgraph_match_counter[i])/float(max(subgraph_match_counter))
    common_graph_pymol.write("set_colour custom_colour_%d = [%f, %f, 0]\n" % (i+1,1-greenery, greenery))
    common_graph_pymol.write("color custom_colour_%d, resi %d\n" % (i+1,i+1))
    
common_graph_pymol.write("color blue, (name N*)\n")
common_graph_pymol.write("color hotpink, (name O*)\n")
common_graph_pymol.write("set bg_rgb = [1,1,1]\n")
common_graph_pymol.write("set ray_trace_mode = 1\n")
common_graph_pymol.write("ray 800,800\n")
common_graph_pymol.write("png %s, dpi=300\n" % (str(output_directory)+'/common_graph.png'))
common_graph_pymol.write("quit")
common_graph_pymol.close()

#output common_graphs.png:
if pymol_flag == True:
    if program_check("pymol") != False:
        subprocess.call(["pymol", output_directory+'/BEST_helix.pdb', "-u", output_directory+'/common_graph.pymol'])
    else:
        print ""
        print "pymol not found. image 'common_graph.png' will not be printed."

#END program:
if verbose_flag == True:
    print ""
    print "terminating script..."

sys.exit()
    




