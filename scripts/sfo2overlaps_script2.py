#!/usr/bin/env python
# -*-coding:Utf-8 -*#
                                
import sys, os, re


__author__ = "Amal Zine El Aabidine"

usage = """

Deduce overlaps from single-single overlaps output by SFO, script2.

"""

#cutoff file_single
NUMBER_SINGLES =  int(sys.argv[1]) 
NUMBER_PAIRES = int(sys.argv[2]) 


NUMBER_DEBUT_REVERSE = NUMBER_SINGLES + NUMBER_PAIRES 
SAME_ORDER = 1
DIFFERENT_ORDER = 2
ss  = "\ts\ts\t\n"
sp  = "\ts\tp\t\n"
pp  = "\tp\tp\n"
MPP = "\t-\t+\t+\t" 
MPM = "\t-\t+\t-\t"
MMP = "\t-\t-\t+\t"
PM  = "\t+\t-\t"
ps  = "\tp\ts\t\n"
MP  = "\t-\t+\t"
PP  = "\t+\t+\t"
## Path to save output files: current working directory
#PATH = os.path.split(os.path.abspath(os.path.realpath(sys.argv[3])))[0]
PATH = os.getcwd()

#input files sorted following the first and then second fields numerically
S_S_I = open(PATH+"/S_S_I_sorted","r")
S_S_N = open(PATH+"/S_S_N_sorted","r")
  
S_F_I_sorted = open(PATH+"/S_F_I_sorted","r") 
S_R_I_sorted = open(PATH+"/S_R_I_sorted","r")

S_F_N_sorted = open(PATH+"/S_F_N_sorted","r")
S_R_N_sorted = open(PATH+"/S_R_N_sorted","r")

F_S_I_sorted = open(PATH+"/F_S_I_sorted","r")
R_S_I_sorted = open(PATH+"/R_S_I_sorted","r")


F_S_N_sorted = open(PATH+"/F_S_N_sorted","r") 
R_S_N_sorted = open(PATH+"/R_S_N_sorted","r")  

F_R_I_sorted = open(PATH+"/F_R_I_sorted","r")   
R_F_I_sorted = open(PATH+"/R_F_I_sorted","r")

F_P_N_sorted = open(PATH+"/F_P_N_sorted","r")
R_P_N_sorted = open(PATH+"/R_P_N_sorted","r")  

S_F_I_SYMETRIC_F_S_I_sorted = open(PATH+"/S_F_I_SYMETRIC_F_S_I_sorted","r")			   			  
S_F_N_SYMETRIC_F_S_N_sorted = open(PATH+"/S_F_N_SYMETRIC_F_S_N_sorted","r")			  
S_R_I_SYMETRIC_R_S_I_sorted = open(PATH+"/S_R_I_SYMETRIC_R_S_I_sorted","r")
S_R_N_SYMETRIC_R_S_N_sorted = open(PATH+"/S_R_N_SYMETRIC_R_S_N_sorted","r")
R_F_N_SYMETRIC_F_R_N_sorted = open(PATH+"/R_F_N_SYMETRIC_F_R_N_sorted","r")

# Outputs files 
S_S_I_formatted = open(PATH+"/S_S_I_formatted9","w")
S_S_N_formatted = open(PATH+"/S_S_N_formatted9","w")

SFI_SRI_formatted = open(PATH+"/SFI_SRI_formatted9","w")

SFN_SRN_formatted = open(PATH+"/SFN_SRN_formatted9","w")

FSI_RSI_formatted = open(PATH+"/FSI_RSI_formatted9","w")
FSN_RSN_formatted = open(PATH+"/FSN_RSN_formatted9","w")

FFN_RRN_formatted = open(PATH+"/FFN_RRN_formatted9","w")
FRI_RFI_formatted = open(PATH+"/FRI_RFI_formatted9","w")


##################################################################################### FUNCTIONS  #####################################################################################

def liste_to_integer(string) :
  string = int(str(string).strip("'[]"))
  return string
#####################################################################################

def extract_ID_from_list_line(list_line):
    string_line = str(list_line).strip("'[]")
    ID = int(string_line.split("\\t")[1])
    return ID
#####################################################################################  
def absolue(string):
  return str(abs(int(string)))
#####################################################################################
def maximum(str1,str2):
  return str(max(int(str1),int(str2)))

#####################################################################################
def is_gauche(string):
  if int(string) >= 0 :
    boolean = True  
  if int(string) < 0 :
    boolean = False

  return boolean
#####################################################################################
def calcul_ordre_overlap(POS1,POS2):
  """
  calcul order of overlaps based on the signs of POS1 and POS2. We calculated the ORD only  for pair-pair  overlaps
  """
  if (int(POS1) >=0 and int(POS2) >= 0) or  (int(POS1) <= 0 and int(POS2) <= 0) :
    ORD = SAME_ORDER
    
  if (int(POS1) > 0 and int(POS2) < 0) or  (int(POS1) < 0 and int(POS2) > 0):

    ORD = DIFFERENT_ORDER
  return ORD 
#####################################################################################

def split_overlap_in_intersection(overlap):
  overlap = str(overlap).strip("[]'")
  liste = overlap.split("\\t")
  return liste
#####################################################################################

def is_reverse(ID_number):
  boolean = False
  if int(ID_number) >= NUMBER_DEBUT_REVERSE :
    boolean= True 
  return boolean

#####################################################################################
def is_forward(ID_number):
  boolean = False
  ID_number=int(ID_number)
  if ID_number < NUMBER_DEBUT_REVERSE and ID_number > NUMBER_SINGLES :
    boolean = True 
  return boolean

#####################################################################################
def is_single(ID_number):
  boolean = False
  ID_number=int(ID_number)
  if ID_number <=  NUMBER_SINGLES :
    boolean = True 
  return boolean

#####################################################################################
def extract_set_F(listeofoverlaps): 
  
  overlaps_F = []
  for overlap in listeofoverlaps:    
    if overlap !=[]:     
      liste = split_overlap_in_intersection(overlap) 
      if is_forward(liste[1]) == True :
	    overlaps_F.append(int(liste[1]))
      else:
	  print "extract_set_F: There is a probleme with the sequence ID %d. It must be comprised between %d (>) and  %d (<=)" % (int(liste[1]),NUMBER_SINGLES,NUMBER_DEBUT_REVERSE)
	  
  return set(overlaps_F)
#####################################################################################
def extract_set_R(liste_line_overlaps): 
 
  overlaps_R = []
  for overlap in liste_line_overlaps:   
    if overlap !=[]:           
      liste = split_overlap_in_intersection(overlap)
      if  is_reverse(liste[1]) == True:
	  overlaps_R.append(int(liste[1]) - NUMBER_PAIRES)
      else:
	  print "extract_set_R There is a probleme with the sequence ID %d. It must be  higher than  %d."  % (int(liste[1]), NUMBER_DEBUT_REVERSE)

  return set(overlaps_R)



#####################################################################################

def extract_set_S(liste_line_overlaps): 
  
  overlaps_single = []
  for overlap in liste_line_overlaps:   
    if overlap !=[]:      
      liste = split_overlap_in_intersection(overlap)
      if  is_single(liste[1]) == True:
	overlaps_single.append(int(liste[1]))		
      else:
	  print "There is a probleme with the sequence ID %d. It must be lesser than  %d."  % (int(liste[1]), NUMBER_SINGLES)	
  return set(overlaps_single)

#####################################################################################

def extract_pos_line(string): 

  string=str(string).strip("'[]")
  line_data = string.split("\\t")
  POS = int(line_data[3])
  
  return POS, line_data

#####################################################################################

def ID1_inf_ID2_inverted_SP(POS_RI,POS_FI): 
  boolean = False
  if POS_RI >= 0 and  POS_FI > 0 and  POS_RI < POS_FI:
    boolean = True
  if POS_RI < 0 and  POS_FI > 0  :    
    boolean = True
  return boolean

#####################################################################################

def ID1_inf_ID2_inverted_PS(POS_RI,POS_FI): 
  boolean = False
  if POS_FI >= 0 and  POS_RI < 0 :
    boolean = True
  if POS_FI < 0 and  POS_RI < 0 and  abs(POS_FI)< abs(POS_RI) :    
    boolean = True
  return boolean
#####################################################################################

def ID1_inf_ID2_normal_SP(POS_FN, POS_RN): 
  boolean = False
  if POS_FN >= 0 and  POS_RN > 0 and  POS_FN < POS_RN:
    boolean = True
  if POS_FN < 0 and  POS_RN > 0  :    
    boolean = True

  return boolean
#####################################################################################
def ID1_inf_ID2_normal_PS(POS_FN, POS_RN): 
  boolean = False
  if POS_FN >= 0 and  POS_RN < 0 :
    boolean = True
  if POS_FN < 0 and  POS_RN < 0 and  abs(POS_FN)< abs(POS_RN) :    
    boolean = True

  return boolean
#####################################################################################
def extract_single_singleN(S_S_N_sorted):  
  """
  Extract single single overlaps in normal orientation and write it to file SSN 
  """
  line = S_S_N_sorted.readline()
  while line != "" :
    list_info=line.split("\t")
    if is_gauche(list_info[3]) == True:
      print "part 0"
      S_S_N_formatted.write(list_info[0]+"\t"+list_info[1]+ "\t" +absolue(list_info[3]) +"\t" +"-"+ MPP + maximum(list_info[8],list_info[9]) +"\t"+"-"  + "\t"+list_info[5]+ "\t"+ "-"+ ss)
      line = S_S_N.readline()
    if is_gauche(list_info[3]) == False:
      print "part 1"
      S_S_N_formatted.write(list_info[1]+"\t"+list_info[0]+ "\t" +absolue(list_info[3]) +"\t" +"-"+ MPP + maximum(list_info[8],list_info[9]) +"\t"+"-"  + "\t"+list_info[5]+ "\t"+ "-"+ ss)
      line = S_S_N.readline()      
      
  S_S_N_formatted.close()  
  S_S_N_sorted.close() 
#####################################################################################
def extract_single_singleI(S_S_I_sorted):
  """
  Extract single single overlaps in inverted orientation and write it to file SSI 
  """  
  line = S_S_I_sorted.readline()
  while line != "" :
    list_info=line.split("\t")

    if is_gauche(list_info[3]) == True:  
      print "part 2"
      S_S_I_formatted.write(list_info[0]+"\t"+list_info[1]+ "\t" + absolue(list_info[3])+"\t" +"-"+ MPM + maximum(list_info[8],list_info[9]) +"\t"+"-" + "\t"+list_info[5]+ "\t"+ "-"+ ss)
      line = S_S_I_sorted.readline()
      
    if is_gauche(list_info[3]) == False:   
      print "part 3"
      S_S_I_formatted.write(list_info[1]+"\t"+list_info[0]+ "\t" + absolue(list_info[3])+"\t" +"-"+ MMP + maximum(list_info[8],list_info[9]) +"\t"+"-" + "\t"+ list_info[5]+ "\t"+ "-"+ ss)
      line = S_S_I_sorted.readline()      
      
  S_S_I_formatted.close()  
  S_S_I_sorted.close()    
###################################################################################### 
def extract_SFI_SRI(result,SFI_SRI_formatted): 
   """
   Extract single paired overlaps in inverted  orientation and write it to file SFI_SRI
   """

   if result[1] !=[] and result[2] !=[] :
      
      overlaps_single_F = extract_set_F(result[1])     
      overlaps_single_R = extract_set_R(result[2])      
      intersection_F_R_Single = overlaps_single_R.intersection(overlaps_single_F)
      if list(intersection_F_R_Single) != []:
	extract_line_data_SFI_SRI(intersection_F_R_Single, SFI_SRI_formatted)
	
   if result[1] == result[2] == [] :
     print "il y a un probleme, les deux listes sont vide :" , result 
#####################################################################################    

def extract_line_data_SFI_SRI(intersection_F_R_Single,  SFI_SRI_formatted):    
  list(intersection_F_R_Single).sort()
  POS2_RI ="" 
  POS1_FI =""
 
  for element1 in intersection_F_R_Single:
    search = element1 + NUMBER_PAIRES
    
    for string1 in result[1]:
      
      if element1 == extract_ID_from_list_line(string1):
	  
	  POS1_FI, line1_FI = extract_pos_line(string1)

    for string2 in result[2]:
    
      if search == extract_ID_from_list_line(string2):
	   POS2_RI,line2_RI = extract_pos_line(string2)


    if POS2_RI !="" and POS1_FI !="" and ID1_inf_ID2_inverted_SP(POS2_RI,line2_RI) ==True :
      line1_FI[9]=line1_FI[9].replace("\\n","")
      line2_RI[9]=line2_RI[9].replace("\\n","")

      maximum_R = maximum(line2_RI[8],line2_RI[9])
      maximum_F = maximum(line1_FI[8],line1_FI[9])
      if is_gauche(line2_RI[3]) == True:
	print "part 4"
	SFI_SRI_formatted.write(line1_FI[0] +"\t" +line1_FI[1]+ "\t"+ absolue(line2_RI[3]) +"\t" + absolue(line1_FI[3]) +MPM + maximum_R  +"\t"+ maximum_F + "\t"+ line2_RI[5]+ "\t" + line1_FI[5] + sp)

      if is_gauche(line2_RI[3]) == False:
	print "part 5"
	SFI_SRI_formatted.write(line1_FI[1] +"\t" +line1_FI[0]+ "\t"+ absolue(line2_RI[3]) +"\t" + absolue(line1_FI[3]) +MMP  + maximum_R  +"\t"+ maximum_F + "\t"+ line2_RI[5]+ "\t" + line1_FI[5] + ps)
    
#####################################################################################
def extract_SFN_SRN(result,SFN_SRN_formatted):  
   if result[1] !=[] and result[2] !=[] :

      overlaps_single_F = extract_set_F(result[1])     
      overlaps_single_R = extract_set_R(result[2])            
      intersection_F_R_Single = overlaps_single_R.intersection(overlaps_single_F)
      if list(intersection_F_R_Single) != []:
        
	extract_line_data_SFN_SRN(intersection_F_R_Single, SFN_SRN_formatted)
	
   if result[1] == result[2] ==[] :
     print "There is a probleme, both liste are empty : " , result 
     
#####################################################################################
def extract_line_data_SFN_SRN(intersection_F_R_Single, SFN_SRN_formatted):    
  POS2_RN = "" 
  POS1_FN = ""
  for element1 in intersection_F_R_Single:
    search = element1 + NUMBER_PAIRES
  
    for string1 in result[1]:
      if element1 == extract_ID_from_list_line(string1):
	   POS1_FN, line1_FN = extract_pos_line(string1)
	    	    
    for string2 in result[2]:
      if search == extract_ID_from_list_line(string2):
	   POS2_RN, line2_RN =  extract_pos_line(string2)
    if POS2_RN != "" and POS1_FN != ""  : 
      line1_FN[9]=line1_FN[9].replace("\\n","")
      line2_RN[9]=line2_RN[9].replace("\\n","")
      maximum_F = maximum(line1_FN[8],line1_FN[9])
      maximum_R = maximum(line2_RN[8],line2_RN[9])
      
      if is_gauche(line1_FN[3]) == True:
	print "part 6"
	SFN_SRN_formatted.write(line1_FN[0] +"\t" +line1_FN[1]+ "\t" + absolue(line1_FN[3]) + "\t"+ absolue(line2_RN[3]) + MPP + maximum_F +"\t"+  maximum_R +"\t"+ line1_FN[5] +"\t"+  line2_RN[5]+  sp)
	
      if is_gauche(line1_FN[3]) == False:
	print "part 7"
	SFN_SRN_formatted.write(line1_FN[1] +"\t" +line1_FN[0]+ "\t" + absolue(line1_FN[3]) + "\t"+ absolue(line2_RN[3]) + MPP + maximum_F +"\t"+  maximum_R +"\t"+ line1_FN[5] +"\t"+  line2_RN[5]+  ps)


#####################################################################################
def extract_FSI_RSI(result,FSI_RSI_formatted):  

   if result[1] !=[] and result[2] !=[] :

      overlaps_single_F = extract_set_S(result[1])     
      overlaps_single_R = extract_set_S(result[2])      
      intersection_F_R_Single = overlaps_single_R.intersection(overlaps_single_F)

      if list(intersection_F_R_Single) != []:
	extract_line_data_FSI_RSI(intersection_F_R_Single, FSI_RSI_formatted)
	
   if result[1] == result[2] == [] :
     print "il y a un probleme , les deux listes sont vide :" , result 
      

def extract_line_data_FSI_RSI(intersection_F_R_Single, FSI_RSI_formatted):    
  list(intersection_F_R_Single).sort()
  POS2_RI = "" 
  POS1_FI = ""
  for element1 in intersection_F_R_Single:
    
    element1 = liste_to_integer(element1)
    for string1 in result[1]:
      if element1 == extract_ID_from_list_line(string1):
	    POS1_FI,line1_FI = extract_pos_line(string1)

	    	    
    for string2 in result[2]:
      if element1== extract_ID_from_list_line(string2):
	   POS2_RI,line2_RI = extract_pos_line(string2)

    if POS2_RI != "" and POS1_FI != "" and  ID1_inf_ID2_inverted_PS(POS2_RI, POS1_FI) ==True : 
      line1_FI[9]=line1_FI[9].replace("\\n","")
      line2_RI[9]=line2_RI[9].replace("\\n","")
      maximum_F = maximum(line1_FI[8],line1_FI[9])
      maximum_R = maximum(line2_RI[8],line2_RI[9]) 
      
      if is_gauche(line2_RI[3]) == True:    
	print "part 8"
	FSI_RSI_formatted.write(line1_FI[0] +"\t" +line1_FI[1]+ "\t" + absolue(line2_RI[3]) + "\t"+ absolue(line1_FI[3])  +MPM + maximum_R +"\t"+ maximum_F +"\t"+ line2_RI[5]+"\t"+line1_FI[5]   + ps)
      if is_gauche(line2_RI[3]) == False:
	print "part 9"
	FSI_RSI_formatted.write(line1_FI[1] +"\t" +line1_FI[0]+ "\t" + absolue(line2_RI[3]) + "\t"+ absolue(line1_FI[3])  +MMP + maximum_R +"\t"+ maximum_F +"\t"+ line2_RI[5]+"\t"+line1_FI[5] + sp)

#####################################################################################
def extract_FSN_RSN(result,FSN_RSN_formatted): 
   
   if result[1] !=[] and result[2] !=[] :

      overlaps_single_F = extract_set_S(result[1])     
      overlaps_single_R = extract_set_S(result[2])      
      
      intersection_F_R_Single = overlaps_single_R.intersection(overlaps_single_F)
      
      if list(intersection_F_R_Single) != []:
	extract_line_data_FSN_RSN(intersection_F_R_Single, FSN_RSN_formatted)
	
   if result[1] == result[2] == [] :
     print "il y a un probleme , les deux listes sont vide :" , result 
          
#####################################################################################
def extract_line_data_FSN_RSN(intersection_F_R_Single, FSN_RSN_formatted):    
  list(intersection_F_R_Single).sort()
  POS2_RN = "" 
  POS1_FN = ""
  for element1 in intersection_F_R_Single:
    
    element1 = liste_to_integer(element1)
    for string1 in result[1]:
      if element1 == extract_ID_from_list_line(string1):
	    POS1_FN,line1_FN = extract_pos_line(string1)
	    	    
    for string2 in result[2]:
      if element1== extract_ID_from_list_line(string2):
	   POS2_RN,line2_RN = extract_pos_line(string2)
    if POS2_RN != "" and POS1_FN != "" and ID1_inf_ID2_normal_PS(POS1_FN, POS2_RN) == True :
      line1_FN[9]=line1_FN[9].replace("\\n","")
      line2_RN[9]=line2_RN[9].replace("\\n","")
      maximum_F = maximum(line1_FN[8],line1_FN[9])
      maximum_R = maximum(line2_RN[8],line2_RN[9]) 
      
      if is_gauche(line1_FN[3]) == True:  
	  print "part 10"
	  FSN_RSN_formatted.write(line1_FN[0] +"\t" +line1_FN[1]+ "\t" + absolue(line1_FN[3]) + "\t"+ absolue(line2_RN[3]) +MPP + maximum_F +"\t"+ maximum_R  +"\t"+  line1_FN[5]  +"\t"+  line2_RN[5]+ ps)

      if is_gauche(line1_FN[3]) == False: 
	  print "part 11"
	  FSN_RSN_formatted.write(line1_FN[1] +"\t" +line1_FN[0]+ "\t" + absolue(line1_FN[3]) + "\t"+ absolue(line2_RN[3]) +MPP + maximum_F +"\t"+ maximum_R  +"\t"+ line1_FN[5]  +"\t"+  line2_RN[5]+ sp)



#####################################################################################

     
def extract_FPN_RFN(result,FFN_RRN_formatted): 

   if result[1] !=[] and result[2] !=[] :    
     overlaps_FFN = extract_set_FN(result[1])   
     overlaps_RRN = extract_set_RN(result[2])  

     intersection_FFN_RRN =  overlaps_FFN.intersection(overlaps_RRN)
     
     if list(intersection_FFN_RRN) != []:
	extract_line_data_FFN_RRN(intersection_FFN_RRN, result, FFN_RRN_formatted)
	
   if result[1] == result[2] == [] :
     print "il y a un probleme, les deux listes sont vide :" , result 
     
#####################################################################################   
def extract_set_FN(listeofoverlaps): 
  
  overlaps_F = []
  for overlap in listeofoverlaps:    
    if overlap !=[]:     
      liste = split_overlap_in_intersection(overlap)
      if is_forward(liste[1]) == True :
	    overlaps_F.append(int(liste[1]))
	  
  return set(overlaps_F)

#####################################################################################
def extract_set_RN(liste_line_overlaps): 
 
  overlaps_R = []
  for overlap in liste_line_overlaps:   
    if overlap !=[]:           
      liste = split_overlap_in_intersection(overlap)

      if  is_reverse(liste[1]) == True:
	  overlaps_R.append(int(liste[1]) - NUMBER_PAIRES)
  return set(overlaps_R)

#####################################################################################
def extract_line_data_FFN_RRN(intersection_FFN_RRN, result, FFN_RRN_formatted):    
  POS1_FFN ="" 
  POS2_RRN = ""
  for element1 in list(intersection_FFN_RRN):

    element1 = liste_to_integer(element1)
    search=element1+NUMBER_PAIRES

    for string1 in result[1]:
      if element1 == extract_ID_from_list_line(string1):
	    POS1_FFN,line1_FFN = extract_pos_line(string1)

 
    for string2 in result[2]:
      if search== extract_ID_from_list_line(string2):
	   POS2_RRN,line2_RRN = extract_pos_line(string2)

    ORD = calcul_ordre_overlap(POS1_FFN,POS2_RRN)

    if POS1_FFN !="" and POS2_RRN != ""  :
      line1_FFN[9]=line1_FFN[9].replace("\\n","")
      line2_RRN[9]=line2_RRN[9].replace("\\n","")
      maximum_F = maximum(line1_FFN[8],line1_FFN[9])
      maximum_R = maximum(line2_RRN[8],line2_RRN[9])
      
      if is_gauche(line1_FFN[3]) == True:  
	print "part 12"
	FFN_RRN_formatted.write(line1_FFN[0] +"\t" + line1_FFN[1]+  "\t" + absolue(line1_FFN[3]) + "\t" + absolue(line2_RRN[3]) +"\t" + str(ORD) + PP + maximum_F  + "\t" + maximum_R + "\t" +line1_FFN[5]+ "\t" + line2_RRN[5]+  pp)
	
      if is_gauche(line1_FFN[3]) == False:     
	print "part 13"
	FFN_RRN_formatted.write(line1_FFN[1] +"\t" + line1_FFN[0]+  "\t" +  absolue(line1_FFN[3]) + "\t" + absolue(line2_RRN[3]) +"\t" + str(ORD) + PP + maximum_F  + "\t" + maximum_R + "\t" +line1_FFN[5]+ "\t" + line2_RRN[5]+ pp)

#####################################################################################
def extract_FPI_RFI(result,FRI_RFI_formatted): 

   if result[1] !=[] and result[2] !=[] :    
      overlaps_FRI = extract_set_R(result[1])   
      overlaps_RFI = extract_set_F(result[2])     
      intersection_FRI_RFI = overlaps_FRI.intersection(overlaps_RFI)

      if list(intersection_FRI_RFI) != []:
	extract_line_data_FRI_RFI(intersection_FRI_RFI, FRI_RFI_formatted)
	
   if result[1] ==  result[2] == [] :
     print "il y a un probleme, les deux listes sont vide :" , result
     
#####################################################################################   
def extract_line_data_FRI_RFI(intersection_FRI_RFI,  FRI_RFI_formatted):  
  POS2_RFI = ""
  POS1_FRI = ""
  list(intersection_FRI_RFI).sort()
  for element1 in intersection_FRI_RFI:

    element1 = liste_to_integer(element1)
    search = element1 + NUMBER_PAIRES
      
    for string1 in result[1]: 
      if search == extract_ID_from_list_line(string1):
	    POS1_FRI,line1_FRI  = extract_pos_line(string1)

    for string2 in result[2]:
      if element1 == extract_ID_from_list_line(string2):
	   POS2_RFI,line2_RFI  = extract_pos_line(string2)

    ORD  = calcul_ordre_overlap(POS1_FRI,POS2_RFI) 

    if POS2_RFI != "" and POS1_FRI != ""  :
      line1_FRI[9]=line1_FRI[9].replace("\\n","")
      line2_RFI[9]=line2_RFI[9].replace("\\n","")
      maximum_R = maximum(line1_FRI[8],line1_FRI[9])
      maximum_F = maximum(line2_RFI[8],line2_RFI[9])
      
      if is_gauche(line1_FRI[3]) == True: 
	print "part 14"
	FRI_RFI_formatted.write(line1_FRI[0] +"\t" + line2_RFI[1] + "\t" + absolue(line1_FRI[3]) + "\t"  + absolue(line2_RFI[3]) +"\t" +str(ORD)+ PM + maximum_R+ "\t"+  maximum_F + "\t"+ line1_FRI[5] + "\t"+line2_RFI[5]+ pp)

      if is_gauche(line1_FRI[3]) == False: 
	print "part 15"
	FRI_RFI_formatted.write(line2_RFI[1] +"\t" + line1_FRI[0] + "\t" + absolue(line1_FRI[3]) + "\t"  + absolue(line2_RFI[3])  +"\t" +str(ORD)+ MP + maximum_R+ "\t"+  maximum_F + "\t"+  line1_FRI[5] + "\t"+line2_RFI[5]+ pp)

 
#####################################################################################

 
def extract_cases_overlaps(file_forward,file_reverse,file_output_formatted):
 
    
    line1 = file_forward.readline()
    line2 = file_reverse.readline()

    if line1 != "" and line2 != "":
      liste1 = line1.split("\t")
      liste2 = line2.split("\t")
      ID_F = int(liste1[0])
      ID_R = int(liste2[0])

      if ID_R >= NUMBER_DEBUT_REVERSE :
	    ID_R_cal = ID_R - NUMBER_PAIRES 
      else:
	    ID_R_cal = ID_R
      

      ID_ref = min(ID_F,ID_R_cal)
      
      EstNewID = True
      
      while line1 != "" and line2 != "":

	  liste1 = line1.split("\t")
	  liste2 = line2.split("\t")
	  
	  ID_F = int(liste1[0])
	  ID_R = int(liste2[0])  
	  
	  if ID_R >= NUMBER_DEBUT_REVERSE :
	    ID_R_cal = ID_R - NUMBER_PAIRES 
	  else:
	    ID_R_cal = ID_R


	  if ID_F == ID_ref == ID_R_cal  :

	      if EstNewID:	      
		  lretour = [ID_ref,[[line1]],[[line2]] ]
	      else:
		  lretour[1].append([line1])
		  lretour[2].append([line2])
		
	      EstNewID = False           
	      line1 = file_forward.readline()
	      line2 = file_reverse.readline()
	  

	  if ID_F == ID_ref  and ID_R_cal != ID_ref :
	    
	      if EstNewID:
		  lretour = [ID_ref ,[[line1]],[]] 
	      else:
		  lretour[1].append([line1])
		  
	      EstNewID = False
	      line1 = file_forward.readline()
	  

	  if ID_F != ID_ref  and ID_R_cal == ID_ref :
	      if EstNewID:
		  lretour = [ID_ref,[],[[line2]]] 
	      else:
		  lretour[2].append([line2])
	  
	      EstNewID = False            
	      line2 = file_reverse.readline()


	  if ID_F != ID_ref  and ID_R_cal != ID_ref :
	    
	      yield lretour

	      ID_ref  = min(ID_F,ID_R_cal)            
	      EstNewID = True
      
      if line1 != "":
	  while line1 != "":
	      liste1 = line1.split("\t")
	      ID_F =int(liste1[0])
	      
	      if ID_F  != ID_ref :
		  yield lretour
		  ID_ref  = ID_F 
		  lretour = [ID_ref ,[[line1]],[]]
	      else:
		  lretour[1].append([line1])
		  
	      line1 = file_forward.readline()
  
	  yield lretour
      
      if line2 != "":
	  while line2 != "":
	      liste2 = line2.split("\t")
	      
	      ID_R_cal = int(liste2[0])- NUMBER_PAIRES
	      
	      if ID_R_cal != ID_ref :
		  yield lretour
		  ID_ref  = ID_R_cal
		  lretour = [ID_ref,[],[[line2]]]
		  
	      else:
		  lretour[2].append([line2])
		  
	      line2 = file_reverse.readline()
	      
	  yield lretour

      yield lretour   
	
#####################################################################################	 
def extract_cases_overlaps_RFI_RFI2(file_reverse,FRI_RFI_formatted):

    
    liste_line= file_reverse.readlines()

    
    cmpt_line=-1

    for line1_RFI in liste_line:
      cmpt_line = cmpt_line +1
      print "-----------------line1  ", cmpt_line 

      print "-----------------\nline1  " ,line1_RFI

      line1_RFI_splitted = line1_RFI.split("\t")  
      ID_R1 = int(line1_RFI_splitted[0])
      ID_F2 = int(line1_RFI_splitted[1])
      ID_F1 = ID_R1 - NUMBER_PAIRES
      ID_R2 = ID_F2 + NUMBER_PAIRES 

      cmpt_line2=cmpt_line+1

      if not cmpt_line2 > len(liste_line)-1:
        
	line2_RFI = liste_line[cmpt_line2]
	line2_RFI_splitted = line2_RFI.split("\t")

        print int(line2_RFI_splitted[0]) > ID_R2
        
	while not int(line2_RFI_splitted[0]) > ID_R2  :
            print "-----------------line2  ", cmpt_line2
            
	    if int(line2_RFI_splitted[0]) == ID_R2 and  int(line2_RFI_splitted[1]) == ID_F1 :
		
		line2_RFI_splitted_symetric =line2_RFI_splitted[:len(line2_RFI_splitted)+1]
		line2_RFI_splitted_symetric[0] = line2_RFI_splitted[0]
		line2_RFI_splitted_symetric[1] = line2_RFI_splitted[1]
		line2_RFI_splitted_symetric[3] =  int(line2_RFI_splitted[4])
		line2_RFI_splitted_symetric[4] =  int(line2_RFI_splitted[3])
		line2_RFI_splitted_symetric[5] = line2_RFI_splitted[5]
		line2_RFI_splitted_symetric[6:] = line2_RFI_splitted[6:]
		
		line1_RFI_splitted[9]=line1_RFI_splitted[9].replace("\\n","")
		line2_RFI_splitted_symetric[9]=line2_RFI_splitted_symetric[9].replace("\\n","")
		maximum_F1 = maximum(line1_RFI_splitted[8],line1_RFI_splitted[9])
		maximum_F2= maximum(line2_RFI_splitted_symetric[8],line2_RFI_splitted_symetric[9])
		ORD  = calcul_ordre_overlap(line1_RFI_splitted[3],line2_RFI_splitted_symetric[3])
		if is_gauche(line1_RFI_splitted[3]) == True:     
		  print "part 18"
		  FRI_RFI_formatted.write(line1_RFI_splitted[1] +"\t" + line2_RFI_splitted_symetric[1] + "\t"   + absolue(line1_RFI_splitted[4]) + "\t"+ absolue(line2_RFI_splitted_symetric[4])  +"\t" +str(ORD)+ PM + maximum_F2+ "\t"+  maximum_F1 + "\t"+ line2_RFI_splitted_symetric[5] + "\t"+line1_RFI_splitted[5]+ pp)
    
		if is_gauche(line1_RFI_splitted[3]) == False:  
		  print "part 19"
		  FRI_RFI_formatted.write(line2_RFI_splitted_symetric[1] +"\t" + line1_RFI_splitted[1] + "\t" + absolue(line2_RFI_splitted_symetric[3]) + "\t"  + absolue(line1_RFI_splitted[3])  +"\t" +str(ORD)+ PM + maximum_F2+ "\t"+  maximum_F1 + "\t"+ line2_RFI_splitted_symetric[5] + "\t"+line1_RFI_splitted[5]+ pp)



	    
	    cmpt_line2=cmpt_line2+1
	    if not cmpt_line2 > len(liste_line)-1:	 	  
	      line2_RFI = liste_line[cmpt_line2]
	      line2_RFI_splitted = line2_RFI.split("\t")
	       


#####################################################################################
def extract_cases_overlaps_FRI_FRI2(file_forward,FRI_RFI_formatted):

    liste_line= file_forward.readlines()

    cmpt_line=-1

    for line1_FRI in liste_line:
      cmpt_line = cmpt_line +1
      
      line1_FRI_splitted = line1_FRI.split("\t")  
      ID_F1 = int(line1_FRI_splitted[0])
      ID_R2 = int(line1_FRI_splitted[1])
      ID_F2 = ID_R2 - NUMBER_PAIRES
      ID_R1 = ID_F1 + NUMBER_PAIRES 

      cmpt_line2=cmpt_line+1

      if not cmpt_line2 > len(liste_line)-1:
        
	line2_FRI = liste_line[cmpt_line2]
	line2_FRI_splitted = line2_FRI.split("\t")
       
	while not int(line2_FRI_splitted[0]) > ID_F2  :

	    if int(line2_FRI_splitted[0]) == ID_F2 and  int(line2_FRI_splitted[1]) == ID_R1 :

		  line1_FRI_splitted[9]=line1_FRI_splitted[9].replace("\\n","")
		  line2_FRI_splitted[9]=line2_FRI_splitted[9].replace("\\n","")
		  maximum_F1 = maximum(line1_FRI_splitted[8],line1_FRI_splitted[9])
		  maximum_F2= maximum(line2_FRI_splitted[8],line2_FRI_splitted[9])
		  line2_FRI_splitted_symetric =line2_FRI_splitted[:len(line2_FRI_splitted)+1]
		  line2_FRI_splitted_symetric[0] = line2_FRI_splitted[0]
		  line2_FRI_splitted_symetric[1] = line2_FRI_splitted[1]
		  line2_FRI_splitted_symetric[3] =  int(line2_FRI_splitted[4])
		  line2_FRI_splitted_symetric[4] =  int(line2_FRI_splitted[3])
		  line2_FRI_splitted_symetric[5] = line2_FRI_splitted[5]
		  line2_FRI_splitted_symetric[6:] = line2_FRI_splitted[6:]
		  

		  ORD  = calcul_ordre_overlap(line1_FRI_splitted[3],line2_FRI_splitted_symetric[3])


		  if is_gauche(line1_FRI_splitted[3]) == True:    
		    print "part 16"
		    FRI_RFI_formatted.write(line1_FRI_splitted[0] +"\t" + line2_FRI_splitted_symetric[0] + "\t" + absolue(line1_FRI_splitted[3]) + "\t"  + absolue(line2_FRI_splitted_symetric[3]) +"\t" +str(ORD)+ PM + maximum_F2+ "\t"+  maximum_F1 + "\t"+ line1_FRI_splitted[5] + "\t"+line2_FRI_splitted_symetric[5]+ pp)

		  if is_gauche(line1_FRI_splitted[3]) == False:  
		    print "part 17"
		    FRI_RFI_formatted.write(line2_FRI_splitted_symetric[0] +"\t" + line1_FRI_splitted[0] + "\t" + absolue(line1_FRI_splitted[3]) + "\t"  + absolue(line2_FRI_splitted_symetric[3])  +"\t" +str(ORD)+ MP + maximum_F2+ "\t"+  maximum_F1 + "\t"+  line1_FRI_splitted[5] + "\t"+line2_FRI_splitted_symetric[5]+ pp)
	    cmpt_line2=cmpt_line2+1
	    if not cmpt_line2 > len(liste_line)-1:	 	  
	      line2_FRI = liste_line[cmpt_line2]
	      line2_FRI_splitted = line2_FRI.split("\t")
	       
#####################################################################################  main  #####################################################################################


print "Starting  ..........."

print "1- ANALYSING S_S_I"
extract_single_singleI(S_S_I)

print "2- ANALYSING S_S_N"
extract_single_singleN(S_S_N)

  
print "3- ANALYSING  S_F_I and S_R_I"
for result in extract_cases_overlaps(S_F_I_sorted,S_R_I_sorted,SFI_SRI_formatted):
  extract_SFI_SRI(result,SFI_SRI_formatted) 
S_F_I_sorted.close()
S_R_I_sorted.close()

S_R_I_sorted = open(PATH+"/S_R_I_sorted","r")  
print "ANALYSING  S_F_I et S_R_I"
for result in extract_cases_overlaps(S_F_I_SYMETRIC_F_S_I_sorted,S_R_I_sorted,SFI_SRI_formatted):
  extract_SFI_SRI(result,SFI_SRI_formatted) 
S_F_I_SYMETRIC_F_S_I_sorted.close() 



  
print "4-1 ANALYSINGal S_F_N and S_R_N" 
for result in extract_cases_overlaps(S_F_N_sorted,S_R_N_sorted,SFN_SRN_formatted):
  extract_SFN_SRN(result,SFN_SRN_formatted)
S_F_N_sorted.close() 
S_F_N_sorted = open(PATH+"/S_F_N_sorted","r")

print "4-2 ANALYSING   S_F_N and S_R_N" 
for result in extract_cases_overlaps(S_F_N_sorted,S_R_N_SYMETRIC_R_S_N_sorted,SFN_SRN_formatted):
  extract_SFN_SRN(result,SFN_SRN_formatted)
S_R_N_SYMETRIC_R_S_N_sorted.close()  


  
print "5-1 ANALYSING F_P_N  and  R_P_N " 
for result in extract_cases_overlaps(F_P_N_sorted, R_P_N_sorted,FFN_RRN_formatted):    
  extract_FPN_RFN(result,FFN_RRN_formatted)
F_P_N_sorted.close()


F_P_N_sorted = open(PATH+"/F_P_N_sorted","r") 
print "5-2 ANALYSING  F_P_N  and R_P_N normal" 
for result in extract_cases_overlaps(F_P_N_sorted, R_F_N_SYMETRIC_F_R_N_sorted,FFN_RRN_formatted):    
  extract_FPN_RFN(result,FFN_RRN_formatted)
  

print "6-1 ANALYSING   F_R_I and R_F_I inversé"   
for result in extract_cases_overlaps(F_R_I_sorted, R_F_I_sorted,FRI_RFI_formatted):      
  extract_FPI_RFI(result,FRI_RFI_formatted)
F_R_I_sorted.close()
R_F_I_sorted.close()

print "6-2  ANALYSING e  F_R_I and R_F_I inversé" 
F_R_I_sorted = open(PATH+"/F_R_I_sorted","r")  
#extract_cases_overlaps_FRI_FRI2(F_R_I_sorted,FRI_RFI_formatted)
F_R_I_sorted.close()

#print "6-3  ANALYSING  F_R_I and R_F_I inversé" 
#F_R_I_sorted = open(PATH+"/F_R_I_sorted","r")  
#R_F_I_sorted = open(PATH+"/R_F_I_sorted","r")  
#extract_cases_overlaps_RFI_RFI2(R_F_I_sorted,FRI_RFI_formatted)


print "7-  ANALYSING  F_S_N  and R_S_N normal" 
for result in extract_cases_overlaps(F_S_N_sorted, R_S_N_sorted,FSN_RSN_formatted): 
  extract_FSN_RSN(result,FSN_RSN_formatted)
  
  
print "8-  ANALYSING   F_S_I  and R_S_I normal" 
for result in extract_cases_overlaps(F_S_I_sorted, R_S_I_sorted,FSI_RSI_formatted):  
  extract_FSI_RSI(result,FSI_RSI_formatted) 
  
