#!/usr/bin/env python
# -*-coding:Utf-8 -*#
                                 
import sys, os, re


__author__ = "Amal Zine El Aabidine"

usage = """

Deduce overlaps from single-single overlaps output by SFO, script1.

"""

NUMBER_SINGLES =  int(sys.argv[1]) 
NUMBER_PAIRES = int(sys.argv[2]) 


ORIENTATION_READ_NORMAL ="N"
ORIENTATION_READ_INVERSE ="I"

#PATH = os.path.split(os.path.abspath(os.path.realpath(sys.argv[3])))[0]
PATH = os.getcwd()
PATH=PATH+"/"
file_Initial = open(sys.argv[3],"r")

S_S_I = open(PATH+"S_S_I","w")
S_S_N = open(PATH+"S_S_N","w")
  
S_F_I = open(PATH+"S_F_I","w")
S_R_I = open(PATH+"S_R_I","w")

S_F_N = open(PATH+"S_F_N","w")
S_R_N = open(PATH+"S_R_N","w")

F_S_I = open(PATH+"F_S_I","w")
S_F_I_SYMETRIC_F_S_I =open(PATH+"S_F_I_SYMETRIC_F_S_I","w")  
			   
			   
F_S_N = open(PATH+"F_S_N","w")  
S_F_N_SYMETRIC_F_S_N =open(PATH+"S_F_N_SYMETRIC_F_S_N","w")  
			   
R_S_I = open(PATH+"R_S_I","w")
S_R_I_SYMETRIC_R_S_I =open(PATH+"S_R_I_SYMETRIC_R_S_I","w")  
			   
R_S_N = open(PATH+"R_S_N","w") 
S_R_N_SYMETRIC_R_S_N =open(PATH+"S_R_N_SYMETRIC_R_S_N","w") 

  
F_R_I = open(PATH+"F_R_I","w")
R_F_I = open(PATH+"R_F_I","w")



F_P_N = open(PATH+"F_P_N","w")  
R_P_N = open(PATH+"R_P_N","w")   

R_F_N_SYMETRIC_F_R_N =open(PATH+"R_F_N_SYMETRIC_F_R_N","w")  

def generate_file():


  print "extracting line"

  ligne = file_Initial.readline()
  while ligne !="" :
    ligne=ligne[:-1]
    ID_LIST=ligne.split("\t")
   
    ID1 = int(ID_LIST[0])
    ID2 = int(ID_LIST[1])
    ORI2 = ID_LIST[2]
    A = int(ID_LIST[3])
    B = int(ID_LIST[4])
    L = int(ID_LIST[5])
    
    if A <= 0 : 
      
      if B < 0: 
	 len_ID1 =  L + abs(B)
	 len_ID2 =  L + abs(A)
	 
      if B >= 0:
	 len_ID1 =  L 
	 len_ID2 =  L + abs(A)+ abs(B)	
	 
    if A > 0 : 
      
      if B <= 0: 
	 len_ID1 =  L + abs(A)+ abs(B)
	 len_ID2 =  L 
 
      if B > 0:	 
	 len_ID1 =  L + abs(A)
	 len_ID2 =  L + abs(B)


    perc_ID1 = (L *100)/len_ID1
    perc_ID2 = (L *100)/len_ID2
    
    

    boolean=True
    if boolean == True :
      
      LIGNE_FORMMATED_TO_WRITE = ligne+ "\t"+str(perc_ID1) +"\t"+str(perc_ID2)+"\t" +str(len_ID1) +"\t" +str(len_ID2)+"\t\n"
      ligne_splitted=ligne.split("\t")
      
      
      LIGNE_SYMETRIC_inverted_TO_WRITE = ligne_splitted[1]+ "\t"+ligne_splitted[0]+ "\t"+ligne_splitted[2]+ "\t"+str(int(ligne_splitted[4])*(1))+ "\t"+str(int(ligne_splitted[3])*(1))+"\t"+ligne_splitted[6]+"\t"+ligne_splitted[5]+"\t"+ligne_splitted[7]+ "\t"+str(perc_ID2) +"\t"+str(perc_ID1)+"\t" +str(len_ID2) +"\t" +str(len_ID1)+"\t\n"
      LIGNE_SYMETRIC_normal_TO_WRITE = ligne_splitted[1]+ "\t"+ligne_splitted[0]+ "\t"+ligne_splitted[2]+ "\t"+str(int(ligne_splitted[3])*(-1))+ "\t"+str(int(ligne_splitted[4])*(-1))+"\t"+ligne_splitted[5]+"\t"+ligne_splitted[6]+"\t"+ligne_splitted[7]+ "\t"+str(perc_ID2) +"\t"+str(perc_ID1)+"\t" +str(len_ID2) +"\t" +str(len_ID1)+"\t\n"
     
      if ID1 <= NUMBER_SINGLES: 
	
	if ID2 <= NUMBER_SINGLES: 
	  
	  if ORI2 == ORIENTATION_READ_NORMAL:
	    S_S_N.write(LIGNE_FORMMATED_TO_WRITE)
	    
	  if ORI2 == ORIENTATION_READ_INVERSE:
	    S_S_I.write(LIGNE_FORMMATED_TO_WRITE)
	    
	if ID2 > NUMBER_SINGLES and ID2 <= (NUMBER_SINGLES + NUMBER_PAIRES) :
	  
	  if ORI2 == ORIENTATION_READ_NORMAL:
	    S_F_N.write(LIGNE_FORMMATED_TO_WRITE)
	    
	  if ORI2 == ORIENTATION_READ_INVERSE :
	    S_F_I.write(LIGNE_FORMMATED_TO_WRITE)
	    
	    
	if ID2 > (NUMBER_SINGLES + NUMBER_PAIRES) :
	  
	  if ORI2 == ORIENTATION_READ_NORMAL:
	    S_R_N.write(LIGNE_FORMMATED_TO_WRITE)
	  
	  if ORI2 == ORIENTATION_READ_INVERSE :
	    S_R_I.write(LIGNE_FORMMATED_TO_WRITE)
	   
	    
      if ID1 > NUMBER_SINGLES and ID1 <= (NUMBER_SINGLES + NUMBER_PAIRES):
	  
	    if ID2 <= NUMBER_SINGLES:
	      
	      if ORI2 == ORIENTATION_READ_NORMAL:
		F_S_N.write(LIGNE_FORMMATED_TO_WRITE)
		S_F_N_SYMETRIC_F_S_N.write(LIGNE_SYMETRIC_normal_TO_WRITE)
		
	      if ORI2 == ORIENTATION_READ_INVERSE:
		F_S_I.write(LIGNE_FORMMATED_TO_WRITE)
		S_F_I_SYMETRIC_F_S_I.write(LIGNE_SYMETRIC_inverted_TO_WRITE)
	    if ID2 > NUMBER_SINGLES  : 
	      
		if ORI2 == ORIENTATION_READ_NORMAL and ID2 <= (NUMBER_SINGLES + NUMBER_PAIRES): 
		  F_P_N.write(LIGNE_FORMMATED_TO_WRITE)
		  
		if ORI2 == ORIENTATION_READ_INVERSE and ID2 > (NUMBER_SINGLES + NUMBER_PAIRES) :
		  F_R_I.write(LIGNE_FORMMATED_TO_WRITE)	 
		  
		 
				
      if ID1 > (NUMBER_SINGLES + NUMBER_PAIRES) :
      
	    if ID2 <= NUMBER_SINGLES:
	      
	      if ORI2 == ORIENTATION_READ_NORMAL:
		R_S_N.write(LIGNE_FORMMATED_TO_WRITE)
		S_R_N_SYMETRIC_R_S_N.write(LIGNE_SYMETRIC_normal_TO_WRITE)
	      if ORI2 == ORIENTATION_READ_INVERSE:
		R_S_I.write(LIGNE_FORMMATED_TO_WRITE)
		S_R_I_SYMETRIC_R_S_I.write(LIGNE_SYMETRIC_inverted_TO_WRITE)
	    if ID2 > NUMBER_SINGLES : 
		if ORI2 == ORIENTATION_READ_NORMAL and ID2 > (NUMBER_SINGLES + NUMBER_PAIRES):
		  R_P_N.write(LIGNE_FORMMATED_TO_WRITE)
		  R_F_N_SYMETRIC_F_R_N.write(LIGNE_SYMETRIC_normal_TO_WRITE)
		if ORI2 == ORIENTATION_READ_INVERSE  and ID2 <= (NUMBER_SINGLES + NUMBER_PAIRES) :
		  R_F_I.write(LIGNE_FORMMATED_TO_WRITE)	
		  
      
    ligne = file_Initial.readline()
   
	  
  
def sort_file():
  os.system("sort -n -k1,1n -k2,2n  %sS_S_I -o %sS_S_I_sorted" %(PATH,PATH))
  os.system("sort -n -k1,1n -k2,2n %sS_S_N -o %sS_S_N_sorted" %(PATH,PATH))
  
  os.system("sort -n -k1,1n -k2,2n %sS_F_I -o %sS_F_I_sorted" %(PATH,PATH))
  os.system("sort -n -k1,1n -k2,2n %sS_F_N -o %sS_F_N_sorted" %(PATH,PATH))
  
  os.system("sort -n -k1,1n -k2,2n %sS_R_I -o %sS_R_I_sorted" %(PATH,PATH))
  os.system("sort -n -k1,1n -k2,2n %sS_R_N -o %sS_R_N_sorted" %(PATH,PATH))
  
  os.system("sort -n -k1,1n -k2,2n %sF_S_I -o %sF_S_I_sorted" %(PATH,PATH))
  os.system("sort -n -k1,1n -k2,2n %sF_S_N -o %sF_S_N_sorted" %(PATH,PATH))

  os.system("sort -n -k1,1n -k2,2n %sR_S_I -o %sR_S_I_sorted" %(PATH,PATH))
  os.system("sort -n -k1,1n -k2,2n %sR_S_N -o %sR_S_N_sorted" %(PATH,PATH))
  
  os.system("sort -n -k1,1n -k2,2n %sF_R_I -o %sF_R_I_sorted" %(PATH,PATH))
  os.system("sort -n -k1,1n -k2,2n  %sF_P_N -o %sF_P_N_sorted" %(PATH,PATH))
  
  os.system("sort -n -k1,1n -k2,2n %sR_F_I -o %sR_F_I_sorted" %(PATH,PATH))
  os.system("sort -n -k1,1n -k2,2n %sR_P_N -o %sR_P_N_sorted" %(PATH,PATH))
  
  
  os.system("sort -n -k1,1n -k2,2n  %sS_F_I_SYMETRIC_F_S_I -o %sS_F_I_SYMETRIC_F_S_I_sorted" %(PATH,PATH))  
  os.system("sort -n -k1,1n -k2,2n %sS_F_N_SYMETRIC_F_S_N -o %sS_F_N_SYMETRIC_F_S_N_sorted" %(PATH,PATH))
  os.system("sort -n -k1,1n -k2,2n %sS_R_I_SYMETRIC_R_S_I -o %sS_R_I_SYMETRIC_R_S_I_sorted" %(PATH,PATH))
  os.system("sort -n -k1,1n -k2,2n %sS_R_N_SYMETRIC_R_S_N -o %sS_R_N_SYMETRIC_R_S_N_sorted" %(PATH,PATH))
  os.system("sort -n -k1,1n -k2,2n %sR_F_N_SYMETRIC_F_R_N -o %sR_F_N_SYMETRIC_F_R_N_sorted" %(PATH,PATH))  




def remove_file():
 
  os.system("rm %sS_S_I" %PATH)
  os.system("rm %sS_S_N " %PATH) 
  os.system("rm %sS_F_I " %PATH)
  os.system("rm %sS_F_N " %PATH)  
  os.system("rm %sS_R_I " %PATH)
  os.system("rm %sS_R_N " %PATH)  
  os.system("rm %sF_S_I " %PATH)
  os.system("rm %sF_S_N " %PATH)
  os.system("rm %sR_S_I " %PATH)
  os.system("rm %sR_S_N " %PATH)  
  os.system("rm %sF_R_I " %PATH)
  os.system("rm  %sF_P_N " %PATH) 
  os.system("rm %sR_F_I " %PATH)
  os.system("rm %sR_P_N " %PATH)
  os.system("rm  %sS_F_I_SYMETRIC_F_S_I " %PATH) 
  os.system("rm %sS_F_N_SYMETRIC_F_S_N " %PATH)
  os.system("rm %sS_R_I_SYMETRIC_R_S_I " %PATH)
  os.system("rm %sS_R_N_SYMETRIC_R_S_N " %PATH)
  os.system("rm %sR_F_N_SYMETRIC_F_R_N " %PATH)  
  

generate_file()
file_Initial.close()
S_S_I.close()
S_S_N.close() 
S_F_I.close()
S_F_N.close()
S_R_I.close()
S_R_N.close()
F_S_I.close()
F_S_N.close()
R_S_I.close()
R_S_N.close()
F_R_I.close()
F_P_N.close()
R_F_I.close()
R_P_N.close()
S_F_I_SYMETRIC_F_S_I.close()			   			  
S_F_N_SYMETRIC_F_S_N.close()			  
S_R_I_SYMETRIC_R_S_I.close()
S_R_N_SYMETRIC_R_S_N.close()
R_F_N_SYMETRIC_F_R_N.close()

sort_file()
remove_file()
