#This class performs operations on each amino acid residue within an inputted PDB file:

import csv
import sys
import copy
import math
import numpy
import string
import random

class PDB_residue:

    #initialize class-wide variables:
    pocket_list = []

    ########################################
    #           INITIALISE CLASS           #
    ########################################
    def __init__(self, IP_residue, naive_residue):

        #initialise first batch of variables:
        self.backup_residue = []
        self.backup_name = []
        self.backup_naive_residue = []

        self.backbone_residue = []
        self.IP_residue = IP_residue
        self.naive_residue = naive_residue
        self.best_residue = []
        self.coordinates = []
        self.name = str(self.naive_residue[0][17:20])
        self.pocket_flag = "NO"
        self.clash_zone = []
        
        self.phi = 0
        self.psi = 0

    ########################################
    #         PUBLIC FUNCTION: Rama        #
    ########################################
    #calculate phi and psi angles for a given residue:
    def Rama(self):

        #SUBFUNCTION: calculate angle between 2 vectors (sign-independent)
        def Dot_angle(vector_1, vector_2):
            dotproduct = vector_1*vector_2
            dotproduct = dotproduct[0][0]+dotproduct[0][1]+dotproduct[0][2]
            length_1 = math.sqrt(math.pow(vector_1[0][0],2)+math.pow(vector_1[0][1],2)+math.pow(vector_1[0][2],2))
            length_2 = math.sqrt(math.pow(vector_2[0][0],2)+math.pow(vector_2[0][1],2)+math.pow(vector_2[0][2],2))
            return float(180)/math.pi*math.acos(float(dotproduct)/float(length_1*length_2))

        #definitions:
        IP_residue = self.IP_residue
        residue_name = self.name
        rama_angles = [0, 0]

        #create dictionaries for all amino acids:
        residue_dictionary = {}
        for i in range(0, len(IP_residue)):
            residue_dictionary[str(IP_residue[i][12:16])] = [float(IP_residue[i][30:38]), float(IP_residue[i][38:46]), float(IP_residue[i][46:54])]
        if residue_name == "PRO":
            NH_vector = numpy.matrix(residue_dictionary[" CD "]) -numpy.matrix(residue_dictionary[" N  "])
        else:
            NH_vector = numpy.matrix(residue_dictionary[" H  "]) -numpy.matrix(residue_dictionary[" N  "]) 
        NCa_vector = numpy.matrix(residue_dictionary[" CA "]) -numpy.matrix(residue_dictionary[" N  "]) 
        CaC_vector = numpy.matrix(residue_dictionary[" C  "]) -numpy.matrix(residue_dictionary[" CA "]) 
        CO_vector = numpy.matrix(residue_dictionary[" O  "]) -numpy.matrix(residue_dictionary[" C  "]) 
            
        NHCa_cross = numpy.cross(NH_vector, NCa_vector)
        NCaC_cross = numpy.cross(NCa_vector, CaC_vector)
        CaCO_cross = numpy.cross(CO_vector, CaC_vector)
            
        rama_angles[0] = Dot_angle(NHCa_cross, NCaC_cross)
        rama_angles[1] = Dot_angle(NCaC_cross, CaCO_cross)
        #determine direction of phi dihedral angle. Sign accordingly:
        phi_sign_cross = numpy.cross(NCa_vector, NHCa_cross)
        phi_sign_angle = Dot_angle(phi_sign_cross, NCaC_cross)
        if rama_angles[0] > 90 and phi_sign_angle < 90:
            rama_angles[0] *=  1
        if rama_angles[0] < 90 and phi_sign_angle < 90:
            rama_angles[0] *=  1
        if rama_angles[0] < 90 and phi_sign_angle > 90:
            rama_angles[0] *=  -1
        if rama_angles[0] > 90 and phi_sign_angle > 90:
            rama_angles[0] *=  -1
        #determine direction of psi dihedral angle. Sign accordingly:
        psi_sign_cross = numpy.cross(CaC_vector, NCaC_cross)
        psi_sign_angle = Dot_angle(psi_sign_cross, CaCO_cross)
        if rama_angles[1] > 90 and psi_sign_angle < 90:
            rama_angles[1] *=  1
        if rama_angles[1] < 90 and psi_sign_angle < 90:
            rama_angles[1] *=  1
        if rama_angles[1] < 90 and psi_sign_angle > 90:
            rama_angles[1] *=  -1
        if rama_angles[1] > 90 and psi_sign_angle > 90:
            rama_angles[1] *=  -1
            
        #output phi/psi angles:    
        self.phi = rama_angles[0]
        self.psi = rama_angles[0]

    ########################################
    #       PRIVATE FUNCTION: Kabsch       #
    ########################################
    #superpose two coordinate sets while minimizing RMSD:
    #NOTE: all outputs are in the form of Python matrices:
    def Kabsch(self, P, Q):
        
        NULL_vector = [-1, -1, -1]
        #standard error messages:
        if len(P) != len(Q):
            print "error: P and Q must be of the same size"
            print "Kabsch function will now terminate"
            return NULL_vector
        
        P_dimensions = len(P)
        Q_dimensions = len(Q)
        if P_dimensions != Q_dimensions:
            print "error: P and Q must have the same number of dimensions"
            print "Kabsch function will now terminate"
            return NULL_vector
        
        for i in range(0, P_dimensions):
            if len(P[i]) != len(Q[0]):
                print "error: P & Q must maintain a constant number of points"
                print "Kabsch function will now terminate"
                return NULL_vector
        
        for i in range(0, Q_dimensions):
            if len(Q[i]) != len(P[0]):
                print "error: Q must maintain a constant number of points in all dimensions"
                print "Kabsch function will now terminate"
                return NULL_vector
        
        #dimensions in space:
        D = P_dimensions
        #number of points:
        N = len(P[1])
        #'m' factor. May be expanded later to incorporate weights:
        m = [float(i)/float(N) for i in [1]*N]
        m = copy.copy([m])
        #calculate centroid of P:
        p0 = numpy.dot(P, numpy.transpose(m))
        #calculate centroid of Q:
        q0 = numpy.dot(Q, numpy.transpose(m))
        #row vector of N ones  
        v1 = [1]*N
        v1 = copy.copy([v1])
        #translate P to origin:
        P = P -numpy.dot(p0, v1)
        #translate Q to origin:
        Q = Q -numpy.dot(q0, v1)
        
    	#C is a covariance matrix of the coordinates
    	#C = P*diag(m)*Q' 
    	#but this is inefficient, involving an N*N matrix, while typically D << N.
    	#so we use another way to compute Pdm = P*diag(m):
        Pdm = [copy.copy([0]*N) for i in range(0,D)]
        for i in range(0, N):
            for j in range(0, D):
                Pdm[j][i] = m[0][i]*P[j][i]
        
        C = numpy.dot(Pdm, numpy.transpose(Q))
        #singular value decomposition:
        V,S,W = numpy.linalg.svd(C)
        #transpose W:
        W = numpy.transpose(W)
        I = numpy.identity(D)
        #more numerically stable than using (det(C) < 0):
        if numpy.linalg.det(numpy.dot(V, numpy.transpose(W))) < 0:
            I[D-1][D-1] = -1
        
        U = numpy.dot(W, numpy.dot(I, numpy.transpose(V)))
        r = q0 -numpy.dot(U, p0)
        
        #P, Q already centered:
        Diff = numpy.dot(U, P) -Q
        
        #RMSD = sqrt(sum(sum(Diff.*Diff))/N):
        RMSD = 0
        for i in range(0,N):
            Diff_column = [0]*D
            for j in range(0,D):
                Diff_column[j] = copy.copy(Diff[j][i])
            RMSD = RMSD +m[0][i]*numpy.dot(numpy.transpose(Diff_column), Diff_column)
        RMSD = math.sqrt(RMSD)
        
        #bundle all results for export from the function:
        OP = [U, r, RMSD]
        return OP

    ########################################
    #   PRIVATE FUNCTION: Atomic_orderer   #
    ########################################
    #order residue-atoms according to given input:
    def Atomic_orderer(self, atomic_order):
        
        #initialize a blank list to hold atom_names & 3D coordinates:
        self.coordinates = copy.copy([0]*len(atomic_order))
        #initialize a blank list to hold complete PDB-data:
        residue_store = copy.copy([0]*len(atomic_order))
        
        #extract residue coordinates in order specified:
        for i in range(0,len(self.naive_residue)):
            for j in range(0,len(atomic_order)):
                if str(self.naive_residue[i][12:16]) == atomic_order[j]:
                    self.coordinates[j] = [str(self.naive_residue[i][12:16]), 
                                           float(self.naive_residue[i][30:38]), 
                                           float(self.naive_residue[i][38:46]), 
                                           float(self.naive_residue[i][46:54])]
                    residue_store[j] = copy.copy(self.naive_residue[i])
        
        #also reorder the naive residue:
        self.naive_residue = copy.copy(residue_store)
    
    ########################################
    #       PRIVATE FUNCTION: Rotator      #
    ########################################
    #perform a rotation of points around an arbitrary axis:
    #NOTE: direction of axis vector is from axis[0] to axis[1].
    #NOTE: this function works on 3 dimensions or less.
    #NOTE: angles must be in radian.
    def Rotator(self, rotate_points, axis, angle):
    
        #turn a right-handed angle into a left-handed angle:
        angle = -angle
        
        axis = numpy.matrix(axis)
        axis.astype(float)
        rotate_points = numpy.matrix(rotate_points)
        rotate_points.astype(float)
        
        #extract translation matrix from axis coordinates:
        translation_matrix = -axis[0]
        
        #translate axis to origin:
        axis[0] = axis[0] +translation_matrix
        axis[1] = axis[1] +translation_matrix
        
        #convert axis to unit vector:
        magnitude = 0
        for i in range(0,len(numpy.array(axis[1])[0])):
            magnitude = magnitude +math.pow(axis[1,i],2)
        magnitude = math.sqrt(magnitude)
        axis = copy.deepcopy(axis/float(magnitude))
        
        #perform translation on all coordinates:
        for i in range(0,len(rotate_points)):
            rotate_points[i] = rotate_points[i] +translation_matrix
        
        #create a rotation matrix based on 'angle' term:
        cos_A = math.cos(angle)
        sin_A = math.sin(angle)
        Ux = axis[1,0]
        Uy = axis[1,1]
        Uz = axis[1,2]
        R = [[cos_A +Ux*Ux*(1-cos_A),
              Ux*Uy*(1-cos_A) -Uz*sin_A,
              Ux*Uz*(1-cos_A) +Uy*sin_A],
              
             [Uy*Ux*(1-cos_A) +Uz*sin_A,
              cos_A +Uy*Uy*(1-cos_A),
              Uy*Uz*(1-cos_A) -Ux*sin_A],
              
             [Uz*Ux*(1-cos_A) -Uy*sin_A,
              Uz*Uy*(1-cos_A) +Ux*sin_A,
              cos_A +Uz*Uz*(1-cos_A)]]
        R = numpy.matrix(R)
        
        #perform rotations on coordinates:
        rotate_points = rotate_points*R
        
        #translate points back to their original location:
        for i in range(0,len(rotate_points)):
            rotate_points[i] = rotate_points[i] -translation_matrix
        
        #output coordinates rotated along a given axis:
        return rotate_points
    
    ########################################
    #     PRIVATE FUNCTION: Chi_rotator    #
    ########################################
    #Chi rotator will perform rotamer creation on individual dihedral angles:
    def Chi_rotator(self, axis, rotate_names, chi_angle):
        
        #extract axis values:
        for i in range(0, len(self.coordinates)):
            if self.coordinates[i][0] == axis[0]:
                axis[0] = [self.coordinates[i][1], 
                           self.coordinates[i][2], 
                           self.coordinates[i][3]]
            if self.coordinates[i][0] == axis[1]:
                axis[1] = [self.coordinates[i][1], 
                           self.coordinates[i][2], 
                           self.coordinates[i][3]]
        
        #extract rotatable points values:
        rotate_points=copy.copy([0]*len(rotate_names))
        for i in range(0, len(self.coordinates)):
            for j in range(0,len(rotate_names)):
                if self.coordinates[i][0] == rotate_names[j]:
                    rotate_points[j] = [self.coordinates[i][1], 
                                        self.coordinates[i][2], 
                                        self.coordinates[i][3]]
        
        #perform rotation about CA-CB axis:
        rotate_points = self.Rotator(rotate_points, axis, chi_angle)
        
        #integrate rotated points with self.coordinates:
        for i in range(0, len(self.coordinates)):
            for j in range(0,len(rotate_names)):
                if self.coordinates[i][0] == rotate_names[j]:
                    self.coordinates[i][1] = rotate_points[j,0]
                    self.coordinates[i][2] = rotate_points[j,1]
                    self.coordinates[i][3] = rotate_points[j,2]
    
    ########################################
    #       PRIVATE FUNCTION: Rotamer      #
    ########################################
    #make a rotamer for any selected amino acid residue, based only on chi angles:
    def Rotamer(self, chi_values):
        
        chi = chi_values

        if self.name == "ALA":
            #rearrange atomic points to this order:
            atomic_order = [" N  ", " CA ", " C  ", " CB ", " HA ", " HB1", " HB2", " HB3"]
            self.Atomic_orderer(atomic_order)
            #END of 'IF == "ALA"' code block.
        
        if self.name == "ARG":
            #rearrange atomic points to this order:
            atomic_order = [" N  ", " CA ", " C  ", " CB ", " CG ", " CD ", " NE ", " CZ ", " NH1", " NH2", " HA ", " HB2", " HB3", " HG2", " HG3", " HD2", " HD3", " HE ", "HH11", "HH12", "HH21", "HH22"]
            self.Atomic_orderer(atomic_order)
            #rotate atoms to form chi-1 angle:
            axis = [" CA ", " CB "]
            rotate_names = [" CG ", " CD ", " NE ", " CZ ", " NH1", " NH2", " HB2", " HB3", " HG2", " HG3", " HD2", " HD3", " HE ", "HH11", "HH12", "HH21", "HH22"]
            self.Chi_rotator(axis, rotate_names, chi[0])
            #rotate atoms to form chi-2 angle:
            axis = [" CB ", " CG "]
            rotate_names = [" CD ", " NE ", " CZ ", " NH1", " NH2", " HG2", " HG3", " HD2", " HD3", " HE ", "HH11", "HH12", "HH21", "HH22"]
            self.Chi_rotator(axis, rotate_names, chi[1])
            #rotate atoms to form chi-3 angle:
            axis = [" CG ", " CD "]
            rotate_names = [" NE ", " CZ ", " NH1", " NH2", " HD2", " HD3", " HE ", "HH11", "HH12", "HH21", "HH22"]
            self.Chi_rotator(axis, rotate_names, chi[2])
            #rotate atoms to form chi-4 angle:
            axis = [" CD ", " NE "]
            rotate_names = [" CZ ", " NH1", " NH2",  " HE ", "HH11", "HH12", "HH21", "HH22"]
            self.Chi_rotator(axis, rotate_names, chi[3])
            #END of 'IF == "ARG"' code block.
        
        if self.name == "ASN":
            #rearrange atomic points to this order:
            atomic_order = [" N  ", " CA ", " C  ", " CB ", " CG ", " OD1", " ND2", " HA ", " HB2", " HB3", "HD21", "HD22"]
            self.Atomic_orderer(atomic_order)
            #rotate atoms to form chi-1 angle:
            axis = [" CA ", " CB "]
            rotate_names = [" CG ", " OD1", " ND2"," HB2", " HB3", "HD21", "HD22"]
            self.Chi_rotator(axis, rotate_names, chi[0])
            #rotate atoms to form chi-2 angle:
            axis = [" CB ", " CG "]
            rotate_names = [" OD1", " ND2", "HD21", "HD22"]
            self.Chi_rotator(axis, rotate_names, chi[1])
            #END of 'IF == "ASN"' code block.
        
        if self.name == "ASP":
            #rearrange atomic points to this order:
            atomic_order = [" N  ", " CA ", " C  ", " CB ", " CG ", " OD1", " OD2", " HA ", " HB2", " HB3"]
            self.Atomic_orderer(atomic_order)
            #rotate atoms to form chi-1 angle:
            axis = [" CA ", " CB "]
            rotate_names = [" CG ", " OD1", " OD2", " HB2", " HB3"]
            self.Chi_rotator(axis, rotate_names, chi[0])
            #rotate atoms to form chi-2 angle:
            axis = [" CB ", " CG "]
            rotate_names = [" OD1", " OD2"]
            self.Chi_rotator(axis, rotate_names, chi[1])
            #END of 'IF == "ASP"' code block.
        
        if self.name == "CYS":
            #rearrange atomic points to this order:
            atomic_order = [" N  ", " CA ", " C  ", " CB ", " SG ", " HA ", " HB2", " HB3", " HG "]
            self.Atomic_orderer(atomic_order)
            #rotate atoms to form chi-1 angle:
            axis = [" CA ", " CB "]
            rotate_names = [" SG ", " HB2", " HB3", " HG "]
            self.Chi_rotator(axis, rotate_names, chi[0])
            #END of 'IF == "CYS"' code block.
        
        if self.name == "GLN":
            #rearrange atomic points to this order:
            atomic_order = [" N  ", " CA ", " C  ", " CB ", " CG ", " CD ", " OE1", " NE2", " HA ", " HB2", " HB3", " HG2", " HG3", "HE21", "HE22"]
            self.Atomic_orderer(atomic_order)
            #rotate atoms to form chi-1 angle:
            axis = [" CA ", " CB "]
            rotate_names = [" CG ", " CD ", " OE1", " NE2", " HB2", " HB3", " HG2", " HG3", "HE21", "HE22"]
            self.Chi_rotator(axis, rotate_names, chi[0])
            #rotate atoms to form chi-2 angle:
            axis = [" CB ", " CG "]
            rotate_names = [" CD ", " OE1", " NE2", " HG2", " HG3", "HE21", "HE22"]
            self.Chi_rotator(axis, rotate_names, chi[1])
            #rotate atoms to form chi-3 angle:
            axis = [" CG ", " CD "]
            rotate_names = [" OE1", " NE2", "HE21", "HE22"]
            self.Chi_rotator(axis, rotate_names, chi[2])
            #END of 'IF == "GLN"' code block.
        
        if self.name == "GLU":
            #rearrange atomic points to this order:
            atomic_order = [" N  ", " CA ", " C  ", " CB ", " CG ", " CD ", " OE1", " OE2", " HA ", " HB2", " HB3", " HG2", " HG3"]
            self.Atomic_orderer(atomic_order)
            #rotate atoms to form chi-1 angle:
            axis = [" CA ", " CB "]
            rotate_names = [" CG ", " CD ", " OE1", " OE2", " HB2", " HB3", " HG2", " HG3"]
            self.Chi_rotator(axis, rotate_names, chi[0])
            #rotate atoms to form chi-2 angle:
            axis = [" CB ", " CG "]
            rotate_names = [" CD ", " OE1", " OE2", " HG2", " HG3"]
            self.Chi_rotator(axis, rotate_names, chi[1])
            #rotate atoms to form chi-3 angle:
            axis = [" CG ", " CD "]
            rotate_names = [" OE1", " OE2"]
            self.Chi_rotator(axis, rotate_names, chi[2])
            #END of 'IF == "GLN"' code block.
        
        if self.name == "GLY":
            #rearrange atomic points to this order:
            atomic_order = [" N  ", " CA ", " C  ", " HA2", " HA3"]
            self.Atomic_orderer(atomic_order)
            #END of 'IF == "GLY"' code block.
        
        if self.name == "HIS":
            #rearrange atomic points to this order:
            atomic_order = [" N  ", " CA ", " C  ", " CB ", " CG ", " ND1", " CD2", " CE1", " NE2", " HA ", " HB2", " HB3", " HD2", " HE1"]
            self.Atomic_orderer(atomic_order)
            #rotate atoms to form chi-1 angle:
            axis = [" CA ", " CB "]
            rotate_names = [" CG ", " ND1", " CD2", " CE1", " NE2", " HB2", " HB3", " HD2", " HE1"]
            self.Chi_rotator(axis, rotate_names, chi[0])
            #rotate atoms to form chi-2 angle:
            axis = [" CB ", " CG "]
            rotate_names = [" ND1", " CD2", " CE1", " NE2", " HD2", " HE1"]
            self.Chi_rotator(axis, rotate_names, chi[1])
            #END of 'IF == "HIS"' code block.
        
        if self.name == "ILE":
            #rearrange atomic points to this order:
            atomic_order = [" N  ", " CA ", " C  ", " CB ", " CG1", " CG2", " CD1", " HA ", " HB ", "HG12", "HG13", "HG21", "HG22", "HG23", "HD11", "HD12", "HD13"]
            self.Atomic_orderer(atomic_order)
            #rotate atoms to form chi-1 angle:
            axis = [" CA ", " CB "]
            rotate_names = [" CG1", " CG2", " CD1", " HB ", "HG12", "HG13", "HG21", "HG22", "HG23", "HD11", "HD12", "HD13"]
            self.Chi_rotator(axis, rotate_names, chi[0])
            #rotate atoms to form chi-2 angle:
            axis = [" CB ", " CG1"]
            rotate_names = [" CD1", "HG12", "HG13", "HD11", "HD12", "HD13"]
            self.Chi_rotator(axis, rotate_names, chi[1])
            #END of 'IF == "ILE"' code block.
        
        if self.name == "LEU":
            #rearrange atomic points to this order:
            atomic_order = [" N  ", " CA ", " C  ", " CB ", " CG ", " CD1", " CD2", " HA ", " HB2", " HB3", " HG ", "HD11", "HD12", "HD13", "HD21", "HD22", "HD23"]
            self.Atomic_orderer(atomic_order)
            #rotate atoms to form chi-1 angle:
            axis = [" CA ", " CB "]
            rotate_names = [" CG ", " CD1", " CD2", " HB2", " HB3", " HG ", "HD11", "HD12", "HD13", "HD21", "HD22", "HD23"]
            self.Chi_rotator(axis, rotate_names, chi[0])
            #rotate atoms to form chi-2 angle:
            axis = [" CB ", " CG "]
            rotate_names = [" CD1", " CD2", " HG ", "HD11", "HD12", "HD13", "HD21", "HD22", "HD23"]
            self.Chi_rotator(axis, rotate_names, chi[1])
            #END of 'IF == "LEU"' code block.
        
        if self.name == "LYS":
            #rearrange atomic points to this order:
            atomic_order = [" N  ", " CA ", " C  ", " CB ", " CG ", " CD ", " CE ", " NZ ", " HA ", " HB2", " HB3", " HG2", " HG3", " HD2", " HD3", " HE2", " HE3", " HZ1", " HZ2", " HZ3"]
            self.Atomic_orderer(atomic_order)
            #rotate atoms to form chi-1 angle:
            axis = [" CA ", " CB "]
            rotate_names = [" CG ", " CD ", " CE ", " NZ ", " HB2", " HB3", " HG2", " HG3", " HD2", " HD3", " HE2", " HE3", " HZ1", " HZ2", " HZ3"]
            self.Chi_rotator(axis, rotate_names, chi[0])
            #rotate atoms to form chi-2 angle:
            axis = [" CB ", " CG "]
            rotate_names = [" CD ", " CE ", " NZ ", " HG2", " HG3", " HD2", " HD3", " HE2", " HE3", " HZ1", " HZ2", " HZ3"]
            self.Chi_rotator(axis, rotate_names, chi[1])
            #rotate atoms to form chi-3 angle:
            axis = [" CG ", " CD "]
            rotate_names = [" CE ", " NZ ", " HD2", " HD3", " HE2", " HE3", " HZ1", " HZ2", " HZ3"]
            self.Chi_rotator(axis, rotate_names, chi[2])
            #rotate atoms to form chi-4 angle:
            axis = [" CD ", " CE "]
            rotate_names = [" NZ ", " HE2", " HE3", " HZ1", " HZ2", " HZ3"]
            self.Chi_rotator(axis, rotate_names, chi[3])
            #END of 'IF == "LYS"' code block.
        
        if self.name == "MET":
            #rearrange atomic points to this order:
            atomic_order = [" N  ", " CA ", " C  ", " CB ", " CG ", " SD ", " CE ", " HA ", " HB2", " HB3", " HG2", " HG3", " HE1", " HE2", " HE3"]
            self.Atomic_orderer(atomic_order)
            #rotate atoms to form chi-1 angle:
            axis = [" CA ", " CB "]
            rotate_names = [" CG ", " SD ", " CE ", " HB2", " HB3", " HG2", " HG3", " HE1", " HE2", " HE3"]
            self.Chi_rotator(axis, rotate_names, chi[0])
            #rotate atoms to form chi-2 angle:
            axis = [" CB ", " CG "]
            rotate_names = [" SD ", " CE ", " HG2", " HG3", " HE1", " HE2", " HE3"]
            self.Chi_rotator(axis, rotate_names, chi[1])
            #rotate atoms to form chi-3 angle:
            axis = [" CG ", " SD "]
            rotate_names = [" CE ", " HE1", " HE2", " HE3"]
            self.Chi_rotator(axis, rotate_names, chi[2])
            #END of 'IF == "MET"' code block.
        
        if self.name == "PHE":
            #rearrange atomic points to this order:
            atomic_order = [" N  ", " CA ", " C  ", " CB ", " CG ", " CD1", " CD2", " CE1", " CE2", " CZ ", " HA ", " HB2", " HB3", " HD1", " HD2", " HE1", " HE2", " HZ "]
            self.Atomic_orderer(atomic_order)
            #rotate atoms to form chi-1 angle:
            axis = [" CA ", " CB "]
            rotate_names = [" CG ", " CD1", " CD2", " CE1", " CE2", " CZ ", " HB2", " HB3", " HD1", " HD2", " HE1", " HE2", " HZ "]
            self.Chi_rotator(axis, rotate_names, chi[0])
            #rotate atoms to form chi-2 angle:
            axis = [" CB ", " CG "]
            rotate_names = [" CD1", " CD2", " CE1", " CE2", " CZ ", " HD1", " HD2", " HE1", " HE2", " HZ "]
            self.Chi_rotator(axis, rotate_names, chi[1])
            #END of 'IF == "PHE"' code block.
        if self.name == "PRO":
            
            #rearrange atomic points to this order:
            atomic_order = [" N  ", " CA ", " C  ", " CB ", " CG ", " CD ", " HA ", " HB2", " HB3", " HG2", " HG3", " HD2", " HD3"]
            self.Atomic_orderer(atomic_order)
            #rotate atoms to form chi-1 angle:
            axis = [" CA ", " CB "]
            rotate_names = [" CG ", " CD ", " HB2", " HB3", " HG2", " HG3", " HD2", " HD3"]
            self.Chi_rotator(axis, rotate_names, chi[0])
            #rotate atoms to form chi-2 angle:
            axis = [" CB ", " CG "]
            rotate_names = [" CD ", " HG2", " HG3", " HD2", " HD3"]
            self.Chi_rotator(axis, rotate_names, chi[1])
            #END of 'IF == "PRO"' code block.
        
        if self.name == "SER":
            #rearrange atomic points to this order:
            atomic_order = [" N  ", " CA ", " C  ", " CB ", " OG ", " HA ", " HB2", " HB3", " HG "]
            self.Atomic_orderer(atomic_order)
            #rotate atoms to form chi-1 angle:
            axis = [" CA ", " CB "]
            rotate_names = [" OG ", " HB2", " HB3", " HG "]
            self.Chi_rotator(axis, rotate_names, chi[0])
            #END of 'IF == "SER"' code block.
        
        if self.name == "THR":
            
            #rearrange atomic points to this order:
            atomic_order = [" N  ", " CA ", " C  ", " CB ", " OG1", " CG2", " HA ", " HB ", " HG1", "HG21", "HG22", "HG23"]
            self.Atomic_orderer(atomic_order)
            #rotate atoms to form chi-1 angle:
            axis = [" CA ", " CB "]
            rotate_names = [" OG1", " CG2", " HB ", " HG1", "HG21", "HG22", "HG23"]
            self.Chi_rotator(axis, rotate_names, chi[0])
            #END of 'IF == "THR"' code block.
        
        if self.name == "TRP":
            #rearrange atomic points to this order:
            atomic_order = [" N  ", " CA ", " C  ", " CB ", " CG ", " CD1", " CD2", " NE1", " CE2", " CE3", " CZ2", " CZ3", " CH2", " HA ", " HB2", " HB3", " HD1", " HE1", " HE3", " HZ2", " HZ3", " HH2"]
            self.Atomic_orderer(atomic_order)
            #rotate atoms to form chi-1 angle:
            axis = [" CA ", " CB "]
            rotate_names = [" CG ", " CD1", " CD2", " NE1", " CE2", " CE3", " CZ2", " CZ3", " CH2", " HB2", " HB3", " HD1", " HE1", " HE3", " HZ2", " HZ3", " HH2"]
            self.Chi_rotator(axis, rotate_names, chi[0])
            #rotate atoms to form chi-2 angle:
            axis = [" CB ", " CG "]
            rotate_names = [" CD1", " CD2", " NE1", " CE2", " CE3", " CZ2", " CZ3", " CH2", " HD1", " HE1", " HE3", " HZ2", " HZ3", " HH2"]
            self.Chi_rotator(axis, rotate_names, chi[1])
            #END of 'IF == "TRP"' code block.
        
        if self.name == "TYR":
            #rearrange atomic points to this order:
            atomic_order = [" N  ", " CA ", " C  ", " CB ", " CG ", " CD1", " CD2", " CE1", " CE2", " CZ ", " OH ", " HA ", " HB2", " HB3", " HD1", " HD2", " HE1", " HE2", " HH "]
            self.Atomic_orderer(atomic_order)
            #rotate atoms to form chi-1 angle:
            axis = [" CA ", " CB "]
            rotate_names = [" CG ", " CD1", " CD2", " CE1", " CE2", " CZ ", " OH ", " HB2", " HB3", " HD1", " HD2", " HE1", " HE2", " HH "]
            self.Chi_rotator(axis, rotate_names, chi[0])
            #rotate atoms to form chi-2 angle:
            axis = [" CB ", " CG "]
            rotate_names = [" CD1", " CD2", " CE1", " CE2", " CZ ", " OH ", " HD1", " HD2", " HE1", " HE2", " HH "]
            self.Chi_rotator(axis, rotate_names, chi[1])
            #END of 'IF == "TYR"' code block.
        
        if self.name == "VAL":
            #rearrange atomic points to this order:
            atomic_order = [" N  ", " CA ", " C  ", " CB ", " CG1", " CG2", " HA ", " HB ", "HG11", "HG12", "HG13", "HG21", "HG22", "HG23"]
            self.Atomic_orderer(atomic_order)
            #rotate atoms to form chi-1 angle:
            axis = [" CA ", " CB "]
            rotate_names = [" CG1", " CG2", " HB ", "HG11", "HG12", "HG13", "HG21", "HG22", "HG23"]
            self.Chi_rotator(axis, rotate_names, chi[0])
            #END of 'IF == "VAL"' code block.
        
        #incorporate new coordinates into 'naive_residue':
        for i in range(0, len(self.naive_residue)):
            self.naive_residue[i] = self.naive_residue[i][0:30]+'% 8.3f' % self.coordinates[i][1]+'% 8.3f' % self.coordinates[i][2]+'% 8.3f' % self.coordinates[i][3]+self.naive_residue[i][54:81]

    ########################################
    #     PRIVATE FUNCTION: Backbone_H     #
    ########################################
    #create a backbone hydrogen atom for given (N, CA, C) coordinates:
    def Backbone_H(self, backend_residue):
        #maintain an internal reference of what a standard backone looks like:
        standard_backbone_H = [[32.744,  60.369,  37.844],  #' N  '
                               [33.089,  61.557,  38.616],  #' CA '
                               [33.625,  59.400,  37.624],  #' C  ' (previous)
                               [31.803,  60.356,  37.505]]  #' H  '
        standard_backbone = copy.copy(standard_backbone_H[0:3])
        
        #create a backbone subset from the residue provided:
        backbone_list = [' N  ', ' CA ']
        IP_backbone = copy.copy([0]*(len(backbone_list)+1))
        for i in range(0, len(self.IP_residue)):
            for j in range(0, len(backbone_list)):
                if self.IP_residue[i][12:16] == backbone_list[j]:
                    IP_backbone[j] = [float(self.IP_residue[i][30:38]), 
                                      float(self.IP_residue[i][38:46]), 
                                      float(self.IP_residue[i][46:54])]
        for i in range(0, len(backend_residue)):
            if backend_residue[i][12:16] == ' C  ':
                    IP_backbone[2] = [float(backend_residue[i][30:38]), 
                                      float(backend_residue[i][38:46]), 
                                      float(backend_residue[i][46:54])]
        
        #pass both (std/IP) backbones through Kabsch's algo:
        [rotation_matrix, translation, RMSD] = self.Kabsch(numpy.transpose(standard_backbone), numpy.transpose(IP_backbone))
        if RMSD >= 0.5:
            print "error: backbone RMSD exceeds 0.5 angstrom. residue mismatch apparent"
            print "'PDB_residue.py' class <Backbone_H> terminating..."
            sys.exit()
        
        #perform rotation/translation on backbone matrix:
        standard_backbone_H = numpy.matrix(numpy.dot(rotation_matrix, numpy.transpose(standard_backbone_H)))
        standard_backbone_H += translation
        standard_backbone_H = numpy.transpose(standard_backbone_H)
        
        #output a PDB-format line containing the backbone hydrogen atom:
        return "ATOM      0  H  "+self.IP_residue[0][16:30]+str('% 8.3f' % standard_backbone_H[3,0]+'% 8.3f' % standard_backbone_H[3,1]+'% 8.3f' % standard_backbone_H[3,2])+"  1.00  0.00           H  \n"

    ########################################
    #        PUBLIC FUNCTION: Alter        #
    ########################################
    #Place a new rotamer upon the protein backbone:
    def Alter(self, IP_naive_residue, backend_residue, chi_values):
        
        #convert chi_values from degree to radians:
        for i in range(0, len(chi_values)):
            chi_values[i] = chi_values[i]*math.pi/180
        
        #replace internal 'naive_residue' with user's choice:
        self.naive_residue = copy.copy(IP_naive_residue)
        original_name = self.name
        self.name = str(self.naive_residue[0][17:20])
        current_name = self.name
        
        #check for proline backbones, alter backbone-H accordingly:
        if current_name == "PRO" and original_name != "PRO":
            #remove the backbone H-atom:
            print "note: XXX to PRO mutation requested. Will now remove backbone_H atom."
            temp_residue = []
            for i in range(0, len(self.IP_residue)):
                if self.IP_residue[i][12:16] != " H  ":
                   temp_residue.append(copy.copy(self.IP_residue[i]))
            self.IP_residue = copy.copy(temp_residue)
        
        if current_name != "PRO" and original_name == "PRO":
            #insert a backbone H-atom:
            print "note: PRO to XXX mutation requested. Will now insert backbone_H atom."
            H_line = self.Backbone_H(backend_residue)
            self.IP_residue.append(copy.copy(H_line))
        
        #create the requested rotamer from naive_residue:
        self.Rotamer(chi_values)

        #Superpose rotated naive residue upon backbone:
        backbone_names = [" N  ", " CA ", " C  "]
        #extract backbone coordinates (IP_residue):
        backbone_IP = copy.copy([0]*len(backbone_names))
        for i in range(0, len(self.IP_residue)):
            for j in range(0, len(backbone_names)):
                if str(self.IP_residue[i][12:16]) == backbone_names[j]:
                    backbone_IP[j] = [float(self.IP_residue[i][30:38]), 
                                      float(self.IP_residue[i][38:46]), 
                                      float(self.IP_residue[i][46:54])]
        
        #extract backbone coordinates (naive_residue). Since these coordinates are 
        #already ordered, I can extract the first 3 points:
        backbone_naive = copy.copy([0]*len(backbone_names))
        for i in range(0, len(backbone_names)):
            backbone_naive[i] = [float(self.naive_residue[i][30:38]), 
                                 float(self.naive_residue[i][38:46]), 
                                 float(self.naive_residue[i][46:54])]
        
        #perform Kabsch superposition:
        [rotation_matrix, translation, RMSD] = self.Kabsch(numpy.transpose(backbone_naive), numpy.transpose(backbone_IP))
        if RMSD >= 0.5:
            print "error: backbone RMSD exceeds 0.5 angstrom. residue mismatch apparent"
            print "'PDB_residue.py' class <Alter> terminating..."
            sys.exit()
        
        #extract all naive_residue coordinates.
        naive_sidechain = []
        for i in range(0, len(self.naive_residue)):
            naive_sidechain.append([float(self.naive_residue[i][30:38]), 
                                    float(self.naive_residue[i][38:46]), 
                                    float(self.naive_residue[i][46:54])])
             
        #perform rotation/translation on sidechain_naive matrix:
        naive_sidechain = numpy.matrix(numpy.dot(rotation_matrix, numpy.transpose(naive_sidechain)))
        naive_sidechain += translation
        naive_sidechain = numpy.transpose(naive_sidechain)
        #NOTE: naive_sidechain coordinates are congruent to self.naive_residue coordinates:
        
        #retain original backbone coordinates, but change residue name:
        if current_name == "PRO":
            current_backbone_names = [" N  ", " CA ", " C  ", " O  "]
        else:
            current_backbone_names = [" N  ", " CA ", " C  ", " O  ", " H  "]
        OP_residue = copy.copy([0]*len(current_backbone_names))
        for i in range(0, len(self.IP_residue)):
            for j in range(0, len(current_backbone_names)):
                if self.IP_residue[i][12:16] == current_backbone_names[j]:
                    OP_residue[j] = "ATOM      0 "+current_backbone_names[j]+" "+self.name+self.IP_residue[i][20:81]
        
        #incorporate naive_sidechain coordinates into 'IP_residue':
        for i in range(3, len(self.naive_residue)):
            OP_residue.append("ATOM      0 "+self.naive_residue[i][12:20]+str(self.IP_residue[0][20:30])+'% 8.3f' % naive_sidechain[i,0]+'% 8.3f' % naive_sidechain[i,1]+'% 8.3f' % naive_sidechain[i,2]+"  1.00  0.00          "+self.naive_residue[i][76:78]+"  \n")

        #change self.IP_residue to altered residue:
        self.IP_residue = copy.copy(OP_residue)

    ########################################
    #        PUBLIC FUNCTION: Backup       #
    ########################################
    def Backup(self):
        self.backup_residue = copy.copy(self.IP_residue)
        self.backup_name = self.name
        self.backup_naive_residue = self.naive_residue
        
    ########################################
    #        PUBLIC FUNCTION: Reload       #
    ########################################
    def Restore(self):
        self.IP_residue = copy.copy(self.backup_residue)
        self.name = self.backup_name
        self.naive_residue = self.backup_naive_residue



