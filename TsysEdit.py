#!/usr/bin/env python3

# -*- coding: utf-8 -*-
"""
THIS SCRIPT IS BASED ON 

            Author:     Daewon Kim at Max-Planck-Institute for Radioastronomy (MPIfR)
            Goal:       To examine/modify Tsys and Gain data of the GMVA antennas
            Version:    2
            Tested:     GMVA 2021A session
            Contact:    dwkim@mpifr-bonn.mpg.de (or dwkimastro@gmail.com)
            Last update:  2022.09.04
            *NOTE*
            - It might be quite complex to use this program for those who are not familiar with VLBI antab files.
            - Data formats of the GMVA antennas are very complicated and different by antennas. 
            - If they have changed the formats, then this program will not work.
            - There are some functions to be updated.
            - The current version can be used for 3mm GMVA observations

AND WAS UPDATED BY

Author:     Lena Debbrecht at Max-Planck-Institute for Radioastronomy (MPIfR)
Tested:     GMVA c232a
Contact:    ldebbrecht@mpifr-bonn.mpg.de (or lenadebbrecht@web.de)
Last update:  2024/MAR
==============================================================================
"""





"""
TypeA: VLBA
TypeB: EF, ON, MH, etc and also KVN (since 2020 KVN has 16 columns instead of 2)
TypeC: PV, NN, GLT 

"""







import time
import os
import warnings
import re 
warnings.filterwarnings("ignore")
import matplotlib
from os import listdir
from os.path import isfile, join
import sys
import shutil
import logging
import copy
import numpy as np
import matplotlib.pyplot as plt
import copy
from scipy.interpolate import interp1d
import time
import logging

# import os
# import warnings
# warnings.filterwarnings("ignore")
# import matplotlib
# from os import listdir
# from os.path import isfile, join
# import sys
# import re

matplotlib.rcParams.update({'legend.borderaxespad': 0.3})
matplotlib.rcParams.update({'legend.fontsize': 9.5})
matplotlib.rcParams.update({'legend.framealpha': 0.5})
matplotlib.rcParams.update({'legend.borderaxespad': 0.3})
matplotlib.rcParams.update({'legend.fontsize': 9.5})
matplotlib.rcParams.update({'legend.framealpha': 0.5})



# from tsys_classesA import TypeA
# from tsys_classesB import TypeB
# from tsys_classesC import TypeC
# from Logger import Logger
# import Logger
"""
set interactive mode
for interactive mode False, no confirmation of the input by user is needed.
If problems occur set interactive mode on to control input 

set also wavelength for interactive mode, usually 3mm
"""
# gain same for each station
# GAIN KYS ELEV FREQ=84000,117000 DPFU=0.05253844
# GAIN KUS ELEV FREQ=84000,117000 DPFU=0.04842397
# GAIN KTN ELEV FREQ=84000,117000 DPFU=0.04981619
interactive_mode = False
if interactive_mode:
    antfile= None
# wavelength='7mm'
wavelength='3mm'
# pref = 'c232'
# pref = 'Test'



class Logger:
    def __init__(self):
        self.logger = None
        # self.pref = None
        # self.pref = 'c232'

    def genlog(self):
        i = 0
        # self.pref ="c232a"
        self.pref = 'c232'
        # self.pref = input("<Q> Give a prefix to the log and output file (e.g., c222EF, c221KT)::  ")
        logfile = f"{self.pref}_{i}_LOG.txt"

        while logfile in os.listdir():
            print("The file name exists! Adding a higher number to the end")
            i += 1
            logfile = f"{self.pref}_{i}_LOG.txt"

        print(logfile)
        Logform = "%(asctime)s :: %(message)s"
        logging.basicConfig(filename=os.getcwd() + "/" + logfile, level=logging.INFO, format=Logform)
        self.logger = logging.getLogger()

    def log(self, *args, logtype='info', sep=' '):
        getattr(self.logger, logtype)(sep.join(str(a) for a in args))
        print(*args)

    def newlog(self):
        self.logger.handlers[0].stream.close()
        self.logger.removeHandler(self.logger.handlers[0])
        del self.logger


class TypeA:
    """
    ONLY FOR VLBA antennas.
    output: self.dall, self.tsys, self.arg, self.antcode, self.gain, self.index, self.totcol, self.freq, and self.legend
    """
    def __init__(self, antabfile, gainfile, interactive_mode, wavelength):
        self.wavelength= wavelength
        self.interactive_mode = interactive_mode
        self.logger = Logger()
        self.logger.genlog()
        self.arg = None
        self.antcode = None      #
        self.antabfile= antabfile
        self.gainfile= gainfile
        self.antcodetrue = None      # e.g., 3mm, HN and SC are not available
        self.session = None
        self.dall = None      #
        self.tsysout = None
        self.tsys0 = None      #
        self.tsys1 = None      #
        self.tarr = None
        self.rcp = None
        self.lcp = None
        self.rlcp = None
        self.exmdat = None
        self.pref = None      #
        # self.pref = pref      #
        self.intp = None
        self.intpout = None
        self.gain = None      #
        self.arg = None      #
        self.index = None      #
        self.totcol = None      #
        self.columnorder = None
        self.freq1 = None      #
        self.freq2 = None      #
        self.legend = None
        self.smip = None
        self.bwidth0 = None
        self.bwidth1 = None
        self.cpdat = None
        self.tsysread()
        # self.vlbacon()


    def parserVLBA(self):

        gain_poly = ''
        gn4 = []
        for i in range(len(self.antcodetrue)):
            self.logger.log("\nAntenna -GAIN-  %.0f /" % int(i+1), len(self.antcodetrue), " -->", self.antcodetrue[i])
            self.logger.log("*********************************")
            if self.antcodetrue[i] not in ['BR','FD','HN','KP','LA','MK','NL','OV','PT','SC']:
                self.logger.log("For GBT or phased-VLA, give their Gain Info separately and manually..")
                continue
            gn3 = []
            for k,j in enumerate(self.z2, 0):
                if ('BAND' in j) and ('TIMERANG' in j) and (self.gn1[self.gn2] in j) and (self.antcodetrue[i] in self.z2[k+3]):
                    gn3.append(self.z2[k:k+4])
            print("print1", self.antcodetrue[i])
            print("print",gn3[-1])
            self.logger.log("The latest data of", self.antcodetrue[i], "(BELOW!) \n:", gn3[-1])


            gain_poly_match = re.search(r'GAIN.*?/', str(gn3[-1]))
            if gain_poly_match:
                gain_poly = gain_poly_match.group(0).strip()
                gain_poly_VLBA = gain_poly.replace('GAIN', f'GAIN {self.antcodetrue[i]}')


            self.logger.log("\n<Q> Find the GAIN Info. (including the POLY, see example below) and write it down here \n (e.g., GAIN PT ALTAZ DPFU=xx,xx POLY=xx,xx,xx /) \n **End the line with '/'")
            if self.interactive_mode: 
                while True:
                    provided_gain = input("\nIs this the right GAIN-POLY line? Do not forget the station code!\n(Press ENTER to confirm or input a different line): \n" + gain_poly_VLBA + "\n: ")
                    if provided_gain.strip() == "":
                        gn91 = gain_poly_VLBA
                        break
                    else:
                        gn91 = provided_gain
                        if not gn91.startswith('GAIN'):
                            print("Please provide a line starting with 'GAIN'.")
                        else:
                            break
            else:
                #  add which gain line was used
                ##############################################
                gn91 = gain_poly_VLBA
            gn91 = gn91 + '\n'
            gn4.append(gn91); del gn91
        self.gain = gn4


        #start INDEX 
        #  index vlba :  TSYS       antcode          INDEX='R1:2','L1:2','R3:4','L3:4','R5:6','L5:6','R7:8','L7:8' /
        # tn2 = []
        # for i in range(len(self.antcodetrue)):
        #     if self.interactive_mode:
        #         self.logger.log("\n--TSYS & INDEX-- Checking the order and number of receivers (also Bandwidth & Central Freqs.)")
        #         for i34 in self.tsys1[i][:18]:
        #             self.logger.log(i34)
        #         self.logger.log("\nAntenna -INDEX- %.0f /" % int(i+1), len(self.antcodetrue), " -->", self.antcodetrue[i])
        #         self.logger.log("*********************************")
        #         tn1 = input("<Q> In the output lines above, check the number of Tsys columns and the R/L order! \n (e.g., TSYS PT INDEX='R1','L1','R2','L2','R3','L3','R4','L4' /) \n (e.g., TSYS BR INDEX='R1:2','L1:2','R3:4','L3:4','R5:6','L5:6','R7:8','L7:8' /) \n **The frequencies should be increasing from top to bottom \n **End the line with '/' \n: ")
        #         tn1 = tn1 + '\n'
        #         tn2.append(tn1); del tn1
        #     else:
        #         self.logger.log("\nThe index line for " + ', '.join(self.antcodetrue) + " is probably:\nINDEX='R1:2','L1:2','R3:4','L3:4','R5:6','L5:6','R7:8','L7:8' /" )
        #         tn1 = "TSYS " + ', '.join(self.antcodetrue) + " INDEX='R1:2','L1:2','R3:4','L3:4','R5:6','L5:6','R7:8','L7:8' /" + '\n'
        #         self.logger.log("\nUsing the following INDEX line:\n" + ''.join(tn2))
        #         tn2.append(tn1); del tn1
        # self.index = tn2
        
        tn2 = []
        for i in range(len(self.antcodetrue)):
            antcode = self.antcodetrue[i]  # Get the antcode for the current iteration
            if self.interactive_mode:
                self.logger.log("\n--TSYS & INDEX-- Checking the order and number of receivers (also Bandwidth & Central Freqs.)")
                for i34 in self.tsys1[i][:18]:
                    self.logger.log(i34)
                self.logger.log("\nAntenna -INDEX- %.0f /" % int(i+1), len(self.antcodetrue), " -->", antcode)
                self.logger.log("*********************************")
                tn1 = input("<Q> In the output lines above, check the number of Tsys columns and the R/L order! \n (e.g., TSYS PT INDEX='R1','L1','R2','L2','R3','L3','R4','L4' /) \n (e.g., TSYS BR INDEX='R1:2','L1:2','R3:4','L3:4','R5:6','L5:6','R7:8','L7:8' /) \n **The frequencies should be increasing from top to bottom \n **End the line with '/' \n: ")
                tn1 = tn1 + '\n'
                tn2.append(tn1)
            else:
                tn1 = "TSYS " + antcode + " INDEX='R1:2','L1:2','R3:4','L3:4','R5:6','L5:6','R7:8','L7:8' /" + '\n'
                self.logger.log("\nThe index line for " + antcode + " is probably:\n" + tn1 )
                tn2.append(tn1)
        self.index = tn2

        for i in range(len(self.antcodetrue)):
            self.logger.log("\n<Q> How many tsys-columns in the data? (e.g., 16 or 8 or 2; for VLBA, it is normally 8) \n ")
            if self.interactive_mode:
                hd_legnum = input("Since it's " + self.antcodetrue[i] + " there are probably 8 tsys-columns, (Press ENTER to confirm or give the right number)\n: ")
                if hd_legnum.strip() == "":
                        hd_legnum = 8
                        hd_legnum = int(hd_legnum)
                        self.totcol = hd_legnum
                        self.logger.log("\nTotal", hd_legnum, "columns of tsys measurements.")
                else:
                        hd_legnum = int(hd_legnum)
                        self.totcol = hd_legnum
                        self.logger.log("\nTotal", hd_legnum, "columns of tsys measurements.")        
                        hd_legnum = int(hd_legnum)
            else:
                self.logger.log("Since it's " + self.antcodetrue[i] + " there are probably 8 tsys-columns\n  ")
                hd_legnum = 8
                hd_legnum = int(hd_legnum)
                self.totcol = hd_legnum
                self.logger.log("\nTotal", hd_legnum, "columns of tsys measurements.")

            lgst = []
            self.logger.log("\n<Q> For plot-legends later, see the following and select your data order by giving its number \n 1 -  R1,R2,R3,R4,L1,L2,L3,L4 \n 2 -  R1,R2,R3,R4,R5,R6,R7,R8,L1,L2,L3,L4,L5,L6,L7,L8  --> (usually EU ants.) \n 3 -  R1,L1,R2,L2,R3,L3,R4,L4                         ---> (VLBA as part of GMVA!!) \n 4 -  RCP,LCP                                      ------> (usually those Special ants.) \n (..not there? then Give it manually here: e.g., R1,L1,L2, ...)")
            if self.interactive_mode:
                while True:
                    provided_indleg = input("Since it's " + self.antcodetrue[i] + " the data order is probably \n3 -  R1,L1,R2,L2,R3,L3,R4,L4 , \n(Press ENTER to confirm or give the right number)\n:  ")
                    if provided_indleg.strip() == "":
                        legtem = 'R1(1-2),L1(1-2),R2(3-4),L2(3-4),R3(5-6),L3(5-6),R4(7-8),L4(7-8)'
                        lgst = legtem.split(',')
                        break
                    else:
                        legtem = input(": ")
                        if len(legtem) == 1:
                            legtem = int(legtem)
                            if legtem == 1:
                                legtem = 'R1,R2,R3,R4,L1,L2,L3,L4'
                            elif legtem == 2:
                                legtem = 'R1,R2,R3,R4,R5,R6,R7,R8,L1,L2,L3,L4,L5,L6,L7,L8'
                            elif legtem == 3:
                                if self.freq1 == '3mm':
                                    legtem = 'R1(1-2),L1(1-2),R2(3-4),L2(3-4),R3(5-6),L3(5-6),R4(7-8),L4(7-8)'
                                elif self.freq1 == '7mm':
                                    legtem = 'R1,L1,R2,L2,R3,L3,R4,L4'
                            elif legtem == 4:
                                legtem = 'RCP,LCP'
                        lgst = legtem.split(',')
            else:
                if self.freq1 == '3mm':
                    legtem = 'R1(1-2),L1(1-2),R2(3-4),L2(3-4),R3(5-6),L3(5-6),R4(7-8),L4(7-8)'
                elif self.freq1 == '7mm':
                    legtem = 'R1,L1,R2,L2,R3,L3,R4,L4'
                self.logger.log("Since it's " + self.antcodetrue[i] + " the data order is probably \n " + legtem)
                lgst = legtem.split(',')
            self.legend = lgst


    def tsysread(self):

        self.logger.log("\n@Task: tsysread")
        self.logger.log("===========================")
        self.arg = self.antabfile
        q1 = open(self.antabfile)
        q2 = q1.readlines()
        antcode = []
        tsysdat1 = []
        tsysdat2 = []
        for i,j in enumerate(q2, 0):
            if (j.startswith("T")) and ("TSYS" in j):
                antcode.append( j.split()[1] )
                for k in q2[i:]:
                    tsysdat1.append( k )
                    if (k.startswith("/")) and (len(k.split())==1):
                        # print(antcode)
                        tsysdat2.append(tsysdat1)
                        del tsysdat1
                        tsysdat1 = []
                        break
                    else:
                        continue

        self.antcode = antcode
        self.dall = copy.deepcopy(tsysdat2)
        self.logger.log("The following antennas are in this session", "\n :", self.antcode)
        fqlit = []
        for i in self.dall[0]:
            if (i.startswith('!')) and ('RCP' in i) and ('MHz' in i) and (i.split()[2] not in fqlit):
                fqlit.append(i.split()[2])
        self.freq2 = fqlit
        self.logger.log("\nObserving wavelengths in the data \n :", self.freq2)


        if self.interactive_mode:
            self.freq = input("<Q> Observing wavelength? (e.g., 3mm) \n:")
            self.session = input("\n<Q> Give the session code of the current input VLBA data? (e.g., 'a' if it is 'c222acal.vlba') \n: ")
            self.logger.log(" ")
        else:
            self.logger.log("\nAn observing wavelength of "+ self.wavelength +" is assumed\n")
            self.freq1 = self.wavelength
            pattern = r'^\w{4}(\w)'
            match = re.search(pattern, self.antabfile)
            if match:
                self.session = match.group(1)
            else:
                return None
            
            self.logger.log("\n<Q> Give the session code of the current input VLBA data (e.g., 'a' if it is 'c222acal.vlba') \n since it's " + self.antabfile + " the session code is probably:\n" + self.session)


        fq99 = [s for s in self.freq2 if self.freq1 != s]
        fitsys1 = []
        fitsys2 = []
        nofreql = []
        for i in range(len(self.antcode)):
            self.logger.log("..working on.. ", self.antcode[i])#; time.sleep(0.6)
            abt = 0
            for o,p in enumerate(tsysdat2[i], 0):
                if self.freq1 in p:
                    for q in range(int(1e6)):
                        #% if the ongoing cycle is the end of the data
                        if '/\n' == tsysdat2[i][o+4]:    #% ending part of each antenna
                            for r in tsysdat2[i][o:o+4]:
                                fitsys1.append(r)
                            fitsys2.append(fitsys1)
                            del fitsys1; fitsys1 = []; abt = 1
                            break
                        fitsys1.append(tsysdat2[i][o])
                        del tsysdat2[i][o]
                        #% still in the middle of the data
                        if ([s for s in fq99 if s in tsysdat2[i][o+1]]) or ([s for s in fq99 if s in tsysdat2[i][o+2]]):
                            break
                        else:
                            continue
                elif (p.startswith("/")) and (len(p.split())==1) and (p == tsysdat2[i][-1]):  # at Last
                    if fitsys1 != []:           # ends with other freq.
                        fitsys2.append(fitsys1)
                        del fitsys1; fitsys1 = []
                        break
                    elif (fitsys1 == []) and (abt == 0):     # no data at this antenna
                        self.logger.log(" !!No", self.freq1, "data found at --> ", self.antcode[i])#; time.sleep(1.0)
                        nofreql.append(self.antcode[i])
                        break
                else:
                    continue
        self.antcodetrue = copy.deepcopy(self.antcode)
        for nnn in nofreql:
            for ooo in self.antcodetrue:
                if nnn == ooo:
                    self.antcodetrue.remove(nnn)
                    self.logger.log("\nAntenna without", self.freq1, "data have been removed from the antcode list..")#; time.sleep(0.6)
            self.tsys0 = fitsys2
        if self.freq1 == '7mm':
            self.logger.log('\n     ***For GMVA Obs., 7mm (43 GHz) band usually has multiple bandwidth: Cal. & Sci. scans***')
            input("\npress 'Enter'")
            wchM = []
            for i in range(len(self.antcode)):
                for k,j in enumerate(self.dall[i], 0):
                    if ('RCP' in j) and ('MHz' in j) and ("M " in j) and (self.freq1 in j):
                        if j.split()[8] in wchM:
                            continue
                        else:
                            wchM.append(j.split()[8])
                            self.logger.log("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
                            self.logger.log("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
                            self.logger.log("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  *Check below*")
                            if k <= len(self.dall[i]):
                                for l in self.dall[i][k-1:k+20]:
                                    self.logger.log(l)
                                self.logger.log("\n*********************************")
                                self.logger.log("Check the scans above and find which one is the Scientific band.")
                                tn85 = input("<Q> Move on? (press 'Enter' to continue)")
            self.bwidth0 = wchM
            self.logger.log("\n<Q> Which one is the Sci. band? (normally the longest one; e.g., '128M' or '64M')")
            self.bwidth1 = input(":")
            self.logger.log("<Q> How many Tsys data columns? (e.g., 8 for 128MHz)")
            thecoln = input(":")
            thecoln = int(thecoln)
    #     data fintering 2
        fint0 = []
        if (self.freq1 == '7mm') and (thecoln == 8):
            fint1 = []
            for i in range(len(self.antcodetrue)):
                for o,p in enumerate(fitsys2[i], 0):
                    if not p.startswith('!') and (len(p.split()) >= 10):
                        if '!' in p:
                            uselespart = p.index('!')
                            p99 = p[:uselespart] + '\n'
                            fint0.append(p99)
                        else:
                            fint0.append(p)
                fint1.append(fint0)
                del fint0; fint0 = []
            self.tsys1 = fint1
        elif (self.freq1 == '7mm') and (thecoln == 4):
            fint1 = []
            for i in range(len(self.antcodetrue)):
                for o,p in enumerate(fitsys2[i], 0):
                    if not p.startswith('!') and (len(p.split()) <= 8):
                        if '!' in p:
                            uselespart = p.index('!')
                            p99 = p[:uselespart] + '\n'
                            fint0.append(p99)
                        else:
                            fint0.append(p)
                fint1.append(fint0)
                del fint0; fint0 = []
            self.tsys1 = fint1
        elif self.freq1 == '3mm':
            fint1 = []
            for i in range(len(self.antcodetrue)):
                for o,p in enumerate(fitsys2[i], 0):
                    if not p.startswith('!'):
                        if '!' in p:
                            uselespart = p.index('!')
                            p99 = p[:uselespart] + '\n'
                            fint0.append(p99)
                        else:
                            fint0.append(p)
                fint1.append(fint0)
                del fint0; fint0 = []
            self.tsys1 = fint1
    #     Give the other information
        z1 = open(self.gainfile)
        print("data:", self.gainfile)
        self.z2 = z1.readlines()
        self.gn1 = []
        for i in self.z2:
            if ('BAND' in i) and ('TIMERANG' in i) and (i.split()[2] not in self.gn1):
                self.gn1.append(i.split()[2])

        self.logger.log("\nSetting some data format!")
        for i,j in enumerate(self.gn1, 0):
            print(" ",i,"-",j)

        if self.interactive_mode:
            self.gn2 = input("\n<Q> Select the number (0~14) of your observing band from above \n:")
            self.gn2 = int(self.gn2)
        else:
            if self.wavelength == '7mm':
                self.logger.log("\nAn observing wavelength of "+ self.wavelength +" is assumed\n")
                self.gn2= 13
            elif self.wavelength == '3mm':
                self.logger.log("\nAn observing wavelength of "+ self.wavelength +" is assumed\n")
                self.gn2 = 14
        #start parser
        self.parserVLBA()

    #     DONE.
        self.logger.log("\nOUTPUT: self.arg, self.antcode, self.antcodetrue, self.dall, self.freq2, self.freq1, self.tsys0, self.tsys1, self.gain, self.index, self.totcol, self.legend, self.bwidth0, self.bwidth1")
        self.logger.log(" **just-in-case**..if Something went wrong (incorrect inputs among the above Pars.)")
        self.logger.log("   ,one can give it manually (e.g., self.xx = !!) or repeat this function again.")






    def dbcheck(self):
        """
        To check if the data contains unreliable values.
        run it like e.g., self.dbcheck()
        output: self.tarr
        """
        self.logger.log("\n@Task: dbcheck")
        self.logger.log("===========================")
        chers = []       # results of the checking
        time99 = []
        realdata = []
        if self.interactive_mode:
            sep1 = input("<Q> Give the order of RCP/LCP in the Tsys columns from left to right, either 1 or 2? \n  1 (EU ant. form: = R1 R2 R3 R4 ... L1 L2 L3 L4 ...) \n  2 (VLBA ant. form: = R1 L1 R2 L2 R3 L3 R4 L4 ...) \n: ")
            sep1 = int(sep1)
            if sep1 == 1:
                self.columnorder = "EU normally --> R1 R2 R3 R4 ... L1 L2 L3 L4 ..."
            elif sep1 == 2:
                self.columnorder = "VLBA normally --> R1 L1 R2 L2 R3 L3 R4 L4 ..."
        else:
            self.logger.log("\nSince it's VLBA the order of RCP/LCP in the Tsys columns from left to right are probably\nR1 L1 R2 L2 R3 L3 R4 L4")
            sep1 = 2
            self.columnorder = "VLBA normally --> R1 L1 R2 L2 R3 L3 R4 L4 ..."
        for i77 in range( len(self.antcodetrue) ):
            self.logger.log("\nChecking Reliability  %.0f /" % int(i77+1), len(self.antcodetrue), " -->", self.antcodetrue[i77])
            self.logger.log("*************************************")
            new_t = []
            simpletest = self.tsys1[i77][0].split()[1]
            if simpletest.count(":") == 2:
                for l in self.tsys1[i77]:
                    tem_d = l.split()[0]        #% day array
                    tem_hms = l.split()[1]       #% hh:mm:ss array
                    #% converting all of them to a day with decimal numbers
                    div_h = tem_hms.split(':')[0]
                    div_m = tem_hms.split(':')[1]
                    div_s = tem_hms.split(':')[2]
                    dayform = (int(div_h)/24.) + (float(div_m)/60. /24.) + (float(div_s)/60. /60. /24.) + int(tem_d)
                    new_t.append(dayform)

            elif simpletest.count(":") == 1:
                for l in self.tsys1[i77]:
                    tem_d = l.split()[0]        #% day array
                    tem_hm = l.split()[1]       #% hh:mm array
                    div_h = tem_hm.split(':')[0]
                    div_m = tem_hm.split(':')[1]
                    dayform = (int(div_h)/24.) + (float(div_m)/60. /24.) + int(tem_d)
                    new_t.append(dayform)
                    #% --> actual x-array for interpolation

            new_t = np.asarray(new_t)
            time99.append(new_t)
            tsyslist = []
            for kk1 in range(len(self.legend)):
                eachcol = np.genfromtxt(self.tsys1[i77], usecols=kk1+2, unpack=True, invalid_raise=False, missing_values='', usemask=False, dtype=float)
                tsyslist.append(eachcol)
                del eachcol

            realdata.append(tsyslist)
            self.logger.log("\n**Now, checking any values of [999, >9999, 0, negative, btw 0 and 10]")
            aa = 0
            if sep1 == 1:       #% first half -> rcp, and then lcp (e.g., EU antennas)
                for i in range(len(tsyslist)):
                    if i < int( len(tsyslist)/2 ):
                        s2 = tsyslist[i]
                        self.logger.log("..First, RCP !")
                        for k1,k2 in enumerate(s2,0):
                            if (k2 == 999) or (k2 == 999.0):
                                self.logger.log("-- 999 found!!")
                                aa = 1
                                continue
                            elif k2 >= 9999:
                                self.logger.log("-- higher than 9999 found!!")
                                aa = 1
                                continue
                            elif k2 == 0:
                                self.logger.log("-- 0 found!!")
                                aa = 1
                                continue
                            elif k2 < 0:
                                self.logger.log("-- negative value found!!")
                                aa = 1
                                continue
                            elif 0 < k2 < 10:
                                self.logger.log("-- 0K < tsys < 10K, unreasonably low value found!!")
                                aa = 1
                                continue
                    else:
                        s2 = tsyslist[i]
                        self.logger.log("..Second, LCP !")
                        for k1,k2 in enumerate(s2,0):
                            if (k2 == 999) or (k2 == 999.0):
                                self.logger.log("-- 999 found!!")
                                aa = 1
                                continue
                            elif k2 >= 9999:
                                self.logger.log("-- higher than 9999 found!!")
                                aa = 1
                                continue
                            elif k2 == 0:
                                self.logger.log("-- 0 found!!")
                                aa = 1
                                continue
                            elif k2 < 0:
                                self.logger.log("-- negative value found!!")
                                aa = 1
                                continue
                            elif 0 < k2 < 10:
                                self.logger.log("-- 0K < tsys < 10K, unreasonably low value found!!")
                                aa = 1
                                continue
            elif sep1 == 2:       #% rcp -> even, and lcp -> odd (e.g., VLBA)
                for i in range(len(tsyslist)):
                    if (i % 2) == 0:     #% even, thus RCP
                        s2 = tsyslist[i]
                        self.logger.log("..First, RCP !")
                        for k1,k2 in enumerate(s2,0):
                            if (k2 == 999) or (k2 == 999.0):
                                self.logger.log("-- 999 found!!")
                                aa = 1
                                continue
                            elif k2 >= 9999:
                                self.logger.log("-- higher than 9999 found!!")
                                aa = 1
                                continue
                            elif k2 == 0:
                                self.logger.log("-- 0 found!!")
                                aa = 1
                                continue
                            elif k2 < 0:
                                self.logger.log("-- negative value found!!")
                                aa = 1
                                continue
                            elif 0 < k2 < 10:
                                self.logger.log("-- 0K < tsys < 10K, unreasonably low value found!!")
                                aa = 1
                                continue
                    else:     #% odd, thus LCP
                        s2 = tsyslist[i]
                        self.logger.log("..Second, LCP !")
                        for k1,k2 in enumerate(s2,0):
                            if (k2 == 999) or (k2 == 999.0):
                                self.logger.log("-- 999 found!!")
                                aa = 1
                                continue
                            elif k2 >= 9999:
                                self.logger.log("-- higher than 9999 found!!")
                                aa = 1
                                continue
                            elif k2 == 0:
                                self.logger.log("-- 0 found!!")
                                aa = 1
                                continue
                            elif k2 < 0:
                                self.logger.log("-- negative value found!!")
                                aa = 1
                                continue
                            elif 0 < k2 < 10:
                                self.logger.log("-- 0K < tsys < 10K, unreasonably low value found!!")
                                aa = 1
                                continue




            if aa == 0:
                self.logger.log("\nResult of the antenna ", self.antcodetrue[i77], "(%s/%s)" % (int(i77+1), len(self.antcodetrue)))
                self.logger.log("--------------------------------")
                self.logger.log("***All measurements seem FINE! Ready to be used in the calibration!***")
                self.logger.log(" -Find an output file: xx_Tsys0.dat, in your working directory.")
                self.logger.log(" -But, STRONGLY RECOMMEND to plot for double-checking the data with showsys().")
                self.tsysout = copy.deepcopy(self.tsys1[i77])
                self.tsysout.append("/\n")
                self.tsysout.insert(0, self.index[i77])
                self.tsysout.insert(0, self.gain[i77])
                outn = self.logger.pref+"_" + self.antcodetrue[i77] + "_" + self.session + "_" + self.freq1 +"_Tsys0.dat"
                with open(outn, 'w') as ft26:
                    ft26.writelines(self.tsysout)
                self.logger.log(" -If the plot looks ok, then, use the output file: xx_Tsys0.dat, here.")
                chers.append("FINE")
            elif aa == 1:
                self.logger.log("\nResult of the antenna ", self.antcodetrue[i77], "(%s/%s)" % (int(i77+1), len(self.antcodetrue)))
                self.logger.log("--------------------------------")
                self.logger.log("***Unreliable value found!***")
                self.logger.log(" -Find an output file: xx_Tsys0.dat, in your working directory.")
                self.logger.log(" -Fix it further with showsys() and intpsys()")
                self.logger.log(" -Then, use self.intp with writetxt()")
                chers.append("WRONG")
                self.tsysout = copy.deepcopy(self.tsys1[i77])
                self.tsysout.append("/\n")
                self.tsysout.insert(0, self.index[i77])
                self.tsysout.insert(0, self.gain[i77])
                outn = self.logger.pref+"_" + self.antcodetrue[i77] + "_" + self.session + "_" + self.freq1 +"_Tsys0.dat"
                with open(outn, 'w') as ft26:
                    ft26.writelines(self.tsysout)
            # ntant = input("\nMove on the next Ant.?  (**just press 'Enter key'**)")
        # for i77 in range(len(self.antcodetrue)):
        #     # outn = self.logger.pref + self.antcodetrue[i77] + "_" + self.session + "_" + self.freq1 + "Tsys0.dat"
        #     destination = os.path.join(self.folder_name, os.path.basename(outn))
        #     os.rename(outn, destination)
        #     self.logger.log("\nMoved file", outn, "to", destination)



        self.exmdat = chers
        self.tarr = time99
        self.rlcp = realdata
        self.logger.log("\nOUTPUT: self.tarr, self.columnorder, self.tsysout, self.rlcp, and self.exmdat")
        self.showsys(antnam=self.antcodetrue[i77], auto=True, Ymax=None, Ymin=None) # ADDED HERE 

        if aa == 1:
            self.logger.log("\nSince unreliable values found, starting intpsys")
            # time.sleep(4)
            self.logger.log("\nStarting intpsys with: " + self.antcodetrue[i77])
            #auto = True in order to plot all antennas 
            self.intpsys(antnam=self.antcodetrue[i77], auto=True, hnot=9999, lnot=10, sepif=False)




 

    def showsys(self, auto=False, antnam=None, Ymax=None, Ymin=None):
        """
        Plot the initial tsys measurements in antab file
        For VLBA, a number of antennas there, set 'auto=True' to plot all Ants. automatically
        e.g.,
        all ant:          showsys(self, auto=True)
        individual ant:   showsys(self, antnam='PT')
        
        auto: a closer look at individual antenna -> True
        antnam: if auto is True, then give here antenna code name (e.g., 'EF')
        Ymax: Y-axis range, Maximum value.
        Ymin: Y-axis range, Minimum value.
        
        """
        self.logger.log("\n@Task: showsys")
        self.logger.log("===========================")
        if auto:
            self.logger.log("*Find a png file in your working directory")
            matplotlib.use('Agg')       # turn the GUI-mode off
            for ao in self.antcodetrue:
                ta11 = self.antcodetrue.index(ao)
                new_t = self.tarr[ta11]
                thedata = self.rlcp[ta11]
                legen8 = self.legend
                colord = self.columnorder
                antcode = ao
                NOFIG = False
                self.plotn1(legen8, thedata, colord, new_t, antcode, self.session, self.freq1, Ymax, Ymin, NOFIG)
                self.logger.log(" ", ao, "is DONE!")
            plt.close()
            # matplotlib.use('Qt5Agg')     # turn the GUI-mode on
        else:
            self.logger.log("*Find a png file in your working directory")
            ta11 = self.antcodetrue.index(antnam)
            new_t = self.tarr[ta11]
            thedata = self.rlcp[ta11]
            legen8 = self.legend
            colord = self.columnorder
            antcode = antnam
            if Ymax:
                NOFIG = True
                self.plotn1(legen8, thedata, colord, new_t, antcode, self.session, self.freq1, Ymax, Ymin, NOFIG)
            else:
                NOFIG = False
                self.plotn1(legen8, thedata, colord, new_t, antcode, self.session, self.freq1, Ymax, Ymin, NOFIG)

        self.logger.log("\n*If the plot looks OK and dbcheck did not report suspicious tsys values,")
        self.logger.log(" --> just use the output file of dbcheck. That's it!")
        self.logger.log("*But if not, then go 'intpsys' and use self.intp instead of self.tsysout")
        self.logger.log("\nOUTPUT: A figure file (i.e., xxxx_Tsys0.png)")




    def intpsys(self, auto=False, antnam=None, hnot=9999, lnot=10, sepif=False):
        """
        Linear interpolation to replace those unreliable tsys measurements.
        For VLBA, a number of antennas there, set 'auto=True' to plot all Ants. automatically
        e.g.,
        all ant:          showsys(self, auto=True)
        individual ant:   showsys(self, antnam='PT')
        """
        self.logger.log("\n@task: intpsys")
        self.logger.log("===========================")
        recd1 = []
        recd2 = []
        if auto:
            self.logger.log("*Find a png file in your working directory")
            matplotlib.use('Agg')       # turn the GUI-mode off
            for ao in self.antcodetrue:
                ta11 = self.antcodetrue.index(ao)
                new_t = self.tarr[ta11]
                thedata = self.rlcp[ta11]
                legen8 = self.legend
                colord = self.columnorder
                antcode = ao
                time0 = self.tsys1[ta11]
                outsys = self.plotn2(legen8, thedata, colord, new_t, antcode, time0, hnot, lnot, sepif, self.session, self.freq1)
                self.logger.log(" ", ao, "is DONE!")
                #! save it
                recd1.append(outsys)
                outsys.append("/\n")
                outsys.insert(0, self.index[ta11])
                outsys.insert(0, self.gain[ta11])
                recd2.append(outsys)
            plt.close()
            # matplotlib.use('Qt5Agg')     # turn the GUI-mode on
        else:
            self.logger.log("*Find a png file in your working directory")
            print(antnam) # BR 
            ta11 = self.antcodetrue.index(antnam) # 0
            print(ta11) # 0
            new_t = self.tarr[ta11]
            thedata = self.rlcp[ta11]
            legen8 = self.legend
            colord = self.columnorder
            antcode = antnam #self.antcode = antnam
            time0 = self.tsys1[ta11]
            outsys = self.plotn2(legen8, thedata, colord, new_t, antcode, time0, hnot, lnot, sepif, self.session, self.freq1)
            #        ! save it
            recd1.append(outsys)
            outsys.append("/\n")
            outsys.insert(0, self.index[ta11])
            outsys.insert(0, self.gain[ta11])
            recd2.append(outsys)
            # print(outsys)
            # should print outsys for every antenna 
        self.intp = copy.deepcopy(recd1)       #% interpolated tsys
        if ((antnam) and (hnot != 9999)) or (sepif == True):
            outn = self.logger.pref+"_" + antnam + "_" + self.session + "_" + self.freq1 +"_Tsys3.dat"
            with open(outn, 'w') as ft7:
                ft7.writelines(self.intp[0])
        else:
            #  add option for new folder here: move files only for aa==1

            self.intpout = copy.deepcopy(recd2)    #% record the whole Tsys columns
            for ao in self.antcodetrue:
                outn = self.logger.pref+"_" + ao + "_" + self.session + "_" + self.freq1 +"_Tsys1.dat"
                ta11 = self.antcodetrue.index(ao)
                with open(outn, 'w') as ft8:
                    ft8.writelines(self.intpout[ta11])
        # self.showsys()
        self.logger.log("\n*Check if the results look OK with those replaced Tsys values,")
        self.logger.log("*--> The results are saved in your working directory automatically!")
        self.logger.log("\nOUTPUT: Figure/Data files (i.e., xxxx_Tsys1.png/dat), self.intp, and self.intpout")



# for show sys 
    def plotn1(self, legen8, thedata, colord, new_t, antcode, sen99, fq99, Ymax, Ymin, NOFIG):
        """
        A function for the class type_A
        Plotting tsys
        """
        fig_pol=plt.figure(figsize=(12,9))
        plt.rcParams['legend.frameon'] = 'False'
        plt.rcParams['xtick.labelsize'] = 11
        plt.rcParams['ytick.labelsize'] = 11
        for i in range(len(legen8)):
            k1 = legen8[i]
            k2 = thedata[i]
            if colord.startswith('EU'):
                if i < int( len(legen8)/2 ):
                    ccc = '#1f77b4'
                else:
                    ccc = '#ff7f0e'
            elif colord.startswith('VLBA'):
                if (i % 2) == 0:
                    ccc = '#1f77b4'
                else:
                    ccc = '#ff7f0e'
            if len(legen8) == 16:
                plt.subplot(4, 4, i+1)
            elif len(legen8) == 8:
                plt.subplot(4, 2, i+1)
            plt.plot(new_t, k2, markersize=1.5, marker='o', linewidth=0.8, label=k1, color=ccc)
            fts2 = 11
            plt.legend(loc=2, fontsize=fts2, columnspacing=0.5, markerscale=1.2, labelspacing=0.3, handletextpad=0.3, borderaxespad=0.2, handlelength=1.5, handleheight=0.5)
            vax = np.round(max(k2), 1)
            vin = np.round(min(k2), 1)
            vav = np.round(np.median(k2), 1)
            vst = np.round(np.std(k2), 1)
            if Ymax:
                if Ymin:
                    plt.ylim([Ymin, Ymax])
                else:
                    plt.ylim([0, Ymax])
            else:
                opz1 = np.linspace(plt.ylim()[0], plt.ylim()[1], 7)
                thegap = np.abs(np.abs(opz1[-1])-np.abs(opz1[-2]))
                plt.ylim(opz1[0]-thegap*0.0, opz1[-1]+thegap*1.3)
            plt.title("*"+antcode+"* : Max={}    Min={}    Med={}    Std={}".format(vax,vin,vav,vst), color='r', fontsize=12)
            plt.xlabel("Decimal day", fontsize=13)
            plt.ylabel("Tsys [K]", fontsize=13)
            plt.minorticks_on()
            plt.tick_params(axis='both', which='major', length=6.5, direction='in', pad=2, color='k')
            plt.tick_params(axis='both', which='minor', length=3, direction='in', pad=2, color='k')
        plt.show()
        plt.tight_layout()
        #  plt.subplots_adjust(left=0.045, bottom=0.050, right=0.975, top=0.965, wspace=0.190, hspace=0.345)
        if NOFIG:
            self.logger.log("\nNo output .png file.")
        else:
            plt.savefig(self.logger.pref+"_" + antcode + "_" + sen99 + "_" + fq99 +"_Tsys0.png", format='png', dpi=150, transparent=False)





# for intpsys 
    def plotn2(self, legen8, thedata, colord, new_t, antcode, time0, hnot, lnot, sepif, sen99, fq99):
        """
        A function for the class type_A
        Interpolating tsys
        """
        fig_pol=plt.figure(figsize=(12,9))
        plt.rcParams['legend.frameon'] = 'False'
        plt.rcParams['xtick.labelsize'] = 11
        plt.rcParams['ytick.labelsize'] = 11
        dali = []
        for i in range(len(legen8)):
            k1 = legen8[i]
            s2 = thedata[i]
            if colord.startswith('EU'):
                if i < int( len(legen8)/2 ):
                    ccc = '#1f77b4'
                else:
                    ccc = '#ff7f0e'
            elif colord.startswith('VLBA'):
                if (i % 2) == 0:
                    ccc = '#1f77b4'
                else:
                    ccc = '#ff7f0e'
            if len(legen8) == 16:
                plt.subplot(4, 4, i+1)
            elif len(legen8) == 8:
                plt.subplot(4, 2, i+1)
            #% linear interpolation
            temp1 = copy.deepcopy(new_t)
            temp2 = copy.deepcopy(s2)
            wrong_index = []
            wrong_times = []   
            if sepif:
                self.logger.log("\n**Set high/low Tsys cutoff for each RCP/LCP IF**")
                self.logger.log("\n-->", legen8[i])
                hgo = input("<Q> Tsys shouldn't EXCEED (press 'Enter' to skip): ")
                if len(hgo) == 0:
                    hgo = 9999
                else:
                    hgo = float(hgo)
                lgo = input("<Q> Tsys shouldn't be BELOW (press 'Enter' to skip): ")
                if len(lgo) == 0:
                    lgo = 10
                else:
                    lgo = float(lgo)
                self.logger.log("\n!! any tsys values of 999, >=YOUR CHOICE, <0, 0, and 0<tsys<YOUR CHOICE will be removed !!")
                for aa,bb in enumerate(temp2, 0):
                    if (bb == 999) or (bb == 999.0):
                        wrong_index.append(aa)
                        wrong_times.append(temp1[aa])
                        continue
                    elif bb >= hgo:
                        wrong_index.append(aa)
                        wrong_times.append(temp1[aa])
                        continue
                    elif bb < 0:
                        wrong_index.append(aa)
                        wrong_times.append(temp1[aa])
                        continue
                    elif bb == 0:
                        wrong_index.append(aa)
                        wrong_times.append(temp1[aa])
                        continue
                    elif 0 < bb < lgo:
                        wrong_index.append(aa)
                        wrong_times.append(temp1[aa])
                        continue
            else:
                self.logger.log("\n!! any tsys values of 999, >=9999, <0, 0, and 0<tsys<10 will be removed !!")
                for aa,bb in enumerate(temp2, 0):
                    if (bb == 999) or (bb == 999.0):
                        wrong_index.append(aa)
                        wrong_times.append(temp1[aa])
                        continue
                    elif bb >= hnot:
                        wrong_index.append(aa)
                        wrong_times.append(temp1[aa])
                        continue
                    elif bb < 0:
                        wrong_index.append(aa)
                        wrong_times.append(temp1[aa])
                        continue
                    elif bb == 0:
                        wrong_index.append(aa)
                        wrong_times.append(temp1[aa])
                        continue
                    elif 0 < bb < lnot:
                        wrong_index.append(aa)
                        wrong_times.append(temp1[aa])
                        continue

            #% x & y array to be used for interpolation
            temp3 = np.asarray([tg for ii, tg in enumerate(temp1) if ii not in wrong_index])
            temp4 = np.asarray([tg for ii, tg in enumerate(temp2) if ii not in wrong_index])
            #=================================
            #% in case the first and/or last measurements are not reasonable!!
            t1 = list(copy.deepcopy(new_t)); t2 = list(copy.deepcopy(temp3))
            z1 = list(copy.deepcopy(s2)); z2 = list(copy.deepcopy(temp4))
            if temp3[0] != new_t[0]:
                self.logger.log("tsys starts with wrong values.. to fill in: using Mean of 15% of the data starting from left to right!")
                fil1 = t1.index(t2[0])
                avg_front = np.mean(temp4[:round(len(t1)*0.15+0.5)])
                for j in range(fil1):
                    z2.insert(0, avg_front)
                    t2.insert(j, t1[j])

            else:
                self.logger.log("Front: good to go (the first data point looks ok)")

            if temp3[-1] != new_t[-1]:
                self.logger.log("tsys ends with wrong values.. to fill in: using Mean of 15% of the data starting from right to left!")
                for k in np.arange(-1, np.negative(len(t1)), -1):
                    if t2[-1] == t1[k]:
                        fil2 = k
                        break

                avg_end = np.mean(temp4[-round(len(t1)*0.15+0.5):])
                for m in np.arange(fil2+1, 0, +1):
                    z2.append(avg_end)
                    t2.append(t1[m])

            else:
                self.logger.log("End: good to go (the last data point looks ok)")
            #=================================
            #% Linear interpolation
            gtg_time = np.asarray(t2)
            gtg_tsys = np.asarray(z2)
            linv = interp1d(gtg_time, gtg_tsys, kind='linear')

            plt.plot(new_t, linv(new_t), markersize=1.5, marker='o', linewidth=0.8, alpha=1.0, zorder=0, color=ccc, label=k1)
            if wrong_index:
                plt.plot(wrong_times, linv(wrong_times), linewidth=0, marker='*', markersize=9, color='#2ca02c', label='Lin. Interp.', alpha=0.8, mew=0, mec='w')
                self.logger.log("\nInterpolation has been performed!!")
            else:
                self.logger.log("\nAll measurements are reasonable!!")

            self.logger.log(" ",k1)
            self.logger.log("  --> Number of unreasonable tsys values:", len(wrong_index))
            self.logger.log("---------------------------------------------")
            resu = np.asarray([round(g, 2) for g in linv(new_t)])   #% for recording later..
            dali.append(resu); del(resu)
            #% basic statistical properties
            vax = np.round(max(linv(new_t)), 1)
            vin = np.round(min(linv(new_t)), 1)
            vav = np.round(np.median(linv(new_t)), 1)
            vst = np.round(np.std(linv(new_t)), 1)
            #% plot layout
            fts3 = 11
            plt.legend(loc=2, fontsize=fts3, ncol=2, columnspacing=0.5, markerscale=1.2, labelspacing=0.3, handletextpad=0.3, borderaxespad=0.2, handlelength=1.5, handleheight=0.5)
            opz1 = np.linspace(plt.ylim()[0], plt.ylim()[1], 7)
            thegap = np.abs(np.abs(opz1[-1])-np.abs(opz1[-2]))
            plt.ylim(opz1[0]-thegap*0.0, opz1[-1]+thegap*1.3)
            plt.title("*"+antcode+"* : Max={}    Min={}    Med={}    Std={}".format(vax,vin,vav,vst), color='r', fontsize=12)
            plt.xlabel("Decimal day", fontsize=13)
            plt.ylabel("Tsys [K]", fontsize=13)
            plt.minorticks_on()
            plt.tick_params(axis='both', which='major', length=6.5, direction='in', pad=2, color='k')
            plt.tick_params(axis='both', which='minor', length=3, direction='in', pad=2, color='k')
        # plt.show()
        plt.tight_layout()
        #  plt.subplots_adjust(left=0.045, bottom=0.050, right=0.975, top=0.965, wspace=0.190, hspace=0.345)
        for blk in range(2):
            self.logger.log("     ")

        self.logger.log(">>>> NEXT ANT!")
        if (hnot != 9999) or (sepif == True):
            plt.savefig(self.logger.pref+"_" + antcode + "_" + sen99 + "_" + fq99 +"_Tsys3.png", format='png', dpi=150, transparent=False)
        else:
            plt.savefig(self.logger.pref+"_" + antcode + "_" + sen99 + "_" + fq99 +"_Tsys1.png", format='png', dpi=150, transparent=False)
        #% make a new list with the final results
        td,tmh = np.genfromtxt(time0, usecols=(0,1), unpack=True, invalid_raise=False, missing_values='', usemask=False, dtype=str)
        finl = []
        for n1 in range(len(td)):
            fin99 = str(td[n1]) + " " + str(tmh[n1])
            for n2 in dali:
                fin99 = fin99 + " " + str(n2[n1])
                if len(fin99.split()) == len(legen8) + 2:
                    fin99 = fin99 + " \n"
                else:
                    continue
            finl.append(fin99)
            del fin99
        return finl





    def somegains(self):
        """
        To check a recent Gain Info. of some special antennas (i.e., NN, GL, PV, and GB; GBT).
        run it like e.g., self.somegains()
        output: print the Info. in terminal
        """
        self.logger.log("\n@Task: somegains")
        self.logger.log("===========================")
        self.logger.log("The following Gain Info. were checked in 2022.")
        self.logger.log("(..easy to find the Info. for other antennas)")
        self.logger.log("\nGAIN NN ALTAZ DPFU=0.414,0.415 POLY=1.0 /")
        self.logger.log("\nGAIN PV ELEV DPFU=0.141,0.142 POLY=0.91437,0.005815,-0.00014101,1.5319e-6,-7.5242e-9 /")
        self.logger.log("\nGAIN GL ALTAZ DPFU=0.032,0.032 POLY=1.0 /")
        self.logger.log("\nGAIN GB ELEV DPFU = 0.85,0.85 POLY=3.20720839e-01,1.74591876e-03,-1.84532794e-05 /")
        self.logger.log("\nOUTPUT: n/a")



    def smth(self, antenna=None, whatdat=None, siglev=2):
        """
        To smooth out the curve by removing OUTLIERS!
        <Param>
            antenna: antenna name code
            whatdat: 1 - raw data
                    2 - the basic interpolation
                    3 - interpolation with a specific cutoff or used 'sepif'
            siglev: confidence level of the sampling intervals; e.g., 2 --> 2*std(the intervals). Default=2.
        output: for all VLBA antennas, ~Tsys2.dat, ~Tsys2.png.
        """
        self.logger.log("\n@task: smth")
        self.logger.log("===========================")
        self.logger.log("*Find a png file of the resultant plot in the current path.")
        fig_pol=plt.figure(figsize=(12,9))
        plt.rcParams['legend.frameon'] = 'False'
        plt.rcParams['xtick.labelsize'] = 11
        plt.rcParams['ytick.labelsize'] = 11
        antindex = self.antcodetrue.index(antenna)
        timearr = self.tarr[antindex]
        dali = []
        for i in range(len(self.legend)):
            k1 = self.legend[i]
            if whatdat == 1:
                s2 = self.rlcp[antindex][i]
            elif whatdat == 2:
                s2 = np.genfromtxt(self.intpout[antindex][2:-1], usecols=i+2)
            elif whatdat == 3:
                s2 = np.genfromtxt(self.intp[0][2:-1], usecols=i+2)
            if self.columnorder.startswith('EU'):
                if i < int( len(self.legend)/2 ):
                    ccc = '#1f77b4'
                else:
                    ccc = '#ff7f0e'
            elif self.columnorder.startswith('VLBA'):
                if (i % 2) == 0:
                    ccc = '#1f77b4'
                else:
                    ccc = '#ff7f0e'
            if len(self.legend) == 16:
                plt.subplot(4, 4, i+1)
            elif len(self.legend) == 8:
                plt.subplot(4, 2, i+1)
            #% linear interpolation
            temp1 = copy.deepcopy(timearr)
            temp2 = copy.deepcopy(s2)
            cutC = siglev * np.std(self.interv(temp2)[1])
            wrong_index = []
            wrong_times = []   
            self.logger.log("\n!! Gonna take the Ourliers out from the data !!")
            for aa,bb in enumerate(temp2, 0):
                if aa == len(temp2)-1:
                    break
                if np.abs(temp2[aa+1]-bb) > cutC:
                    if aa < round(len(temp2)/2):
                        half1 = np.mean(temp2[aa:aa+6])  # get a local mean value and judge the data.
                        temarr9 = np.asarray([bb, temp2[aa+1]])
                        inddex = np.abs(temarr9 - half1).argmax()
                        if temarr9[inddex] == bb:
                            if aa in wrong_index:
                                continue
                            else:
                                wrong_index.append(aa)
                                wrong_times.append(temp1[aa])
                                continue
                        else:
                            if aa+1 in wrong_index:
                                continue
                            else:
                                wrong_index.append(aa+1)
                                wrong_times.append(temp1[aa+1])
                                continue
                    else:
                        half2 = np.mean(temp2[aa-5:aa+1])
                        temarr9 = np.asarray([bb, temp2[aa+1]])
                        inddex = np.abs(temarr9 - half2).argmax()
                        if temarr9[inddex] == bb:
                            if aa in wrong_index:
                                continue
                            else:
                                wrong_index.append(aa)
                                wrong_times.append(temp1[aa])
                                continue
                        else:
                            if aa+1 in wrong_index:
                                continue
                            else:
                                wrong_index.append(aa+1)
                                wrong_times.append(temp1[aa+1])
                                continue
            #% x & y array to be used for interpolation
            temp3 = np.asarray([tg for ii, tg in enumerate(temp1) if ii not in wrong_index])
            temp4 = np.asarray([tg for ii, tg in enumerate(temp2) if ii not in wrong_index])
            #=================================
            #% in case the first and/or last measurements are not reasonable!!
            t1 = list(copy.deepcopy(timearr)); t2 = list(copy.deepcopy(temp3))
            z1 = list(copy.deepcopy(s2)); z2 = list(copy.deepcopy(temp4))
            if temp3[0] != timearr[0]:
                self.logger.log("tsys starts with wrong values.. to fill in: using Mean of 15% of the data starting from left to right!")
                fil1 = t1.index(t2[0])
                avg_front = np.mean(temp4[:round(len(t1)*0.15+0.5)])
                for j in range(fil1):
                    z2.insert(0, avg_front)
                    t2.insert(j, t1[j])

            else:
                self.logger.log("Front: good to go (the first data point looks ok)")

            if temp3[-1] != timearr[-1]:
                self.logger.log("tsys ends with wrong values.. to fill in: using Mean of 15% of the data starting from right to left!")
                for k in np.arange(-1, np.negative(len(t1)), -1):
                    if t2[-1] == t1[k]:
                        fil2 = k
                        break

                avg_end = np.mean(temp4[-round(len(t1)*0.15+0.5):])
                for m in np.arange(fil2+1, 0, +1):
                    z2.append(avg_end)
                    t2.append(t1[m])

            else:
                self.logger.log("End: good to go (the last data point looks ok)")
            #=================================
            #% Linear interpolation
            gtg_time = np.asarray(t2)
            gtg_tsys = np.asarray(z2)
            linv = interp1d(gtg_time, gtg_tsys, kind='linear')

            plt.plot(timearr, linv(timearr), markersize=1.5, marker='o', linewidth=0.8, alpha=1.0, zorder=0, color=ccc, label=k1)
            if wrong_index:
                plt.plot(wrong_times, linv(wrong_times), linewidth=0, marker='*', markersize=9, color='#2ca02c', label='Lin. Interp.', alpha=0.8, mew=0, mec='w')
                self.logger.log("\nInterpolation has been performed!!")
            else:
                self.logger.log("\nAll measurements are reasonable!!")

            self.logger.log(" ",k1)
            self.logger.log("  --> Number of OUTLYING tsys values:", len(wrong_index))
            self.logger.log("---------------------------------------------")
            resu = np.asarray([round(g, 2) for g in linv(timearr)])   #% for recording later..
            dali.append(resu); del(resu)
            #% basic statistical properties
            vax = np.round(max(linv(timearr)), 1)
            vin = np.round(min(linv(timearr)), 1)
            vav = np.round(np.median(linv(timearr)), 1)
            vst = np.round(np.std(linv(timearr)), 1)
            #% plot layout
            fts3 = 11
            plt.legend(loc=2, fontsize=fts3, ncol=2, columnspacing=0.5, markerscale=1.2, labelspacing=0.3, handletextpad=0.3, borderaxespad=0.2, handlelength=1.5, handleheight=0.5)
            opz1 = np.linspace(plt.ylim()[0], plt.ylim()[1], 7)
            thegap = np.abs(np.abs(opz1[-1])-np.abs(opz1[-2]))
            plt.ylim(opz1[0]-thegap*0.0, opz1[-1]+thegap*1.3)
            plt.title("*"+antenna+"* : Max={}    Min={}    Med={}    Std={}".format(vax,vin,vav,vst), color='r', fontsize=12)
            plt.xlabel("Decimal day", fontsize=13)
            plt.ylabel("Tsys [K]", fontsize=13)
            plt.minorticks_on()
            plt.tick_params(axis='both', which='major', length=6.5, direction='in', pad=2, color='k')
            plt.tick_params(axis='both', which='minor', length=3, direction='in', pad=2, color='k')
        plt.show()
        plt.tight_layout()
        plt.savefig(self.logger.pref + antenna + "_" + self.session + "_" + self.freq1 +"Tsys2"+"_"+str(siglev)+"-sig"+".png", format='png', dpi=150, transparent=False)
        #% make a new list with the final results
        td,tmh = np.genfromtxt(self.tsys1[antindex], usecols=(0,1), unpack=True, invalid_raise=False, missing_values='', usemask=False, dtype=str)
        finl = []
        for n1 in range(len(td)):
            fin99 = str(td[n1]) + " " + str(tmh[n1])
            for n2 in dali:
                fin99 = fin99 + " " + str(n2[n1])
                if len(fin99.split()) == len(self.legend) + 2:
                    fin99 = fin99 + " \n"
                else:
                    continue
            finl.append(fin99)
            del fin99
        finl.append("/\n")
        finl.insert(0, self.index[antindex])
        finl.insert(0, self.gain[antindex])
        self.smip = copy.deepcopy(finl)
        outn = self.logger.pref + antenna + "_" + self.session + "_" + self.freq1 +"Tsys2"+"_"+str(siglev)+"-sig"+".dat"
        with open(outn, 'w') as ft57:
            ft57.writelines(finl)
        self.logger.log("\n*If the result looks OK with the smoothed out data,")
        self.logger.log("*--> Use the output 'xxx_Tsys2.dat'. That's it!")
        self.logger.log("\nOUTPUT: Tsys2.dat file, self.smip, and a png file (i.e., xxxx_Tsys2.png)")






    def cpif(self, antenna=None, whatdat=None):
        """
        Copy and paste from a single IF with good Tsys to the IFs with bad Tsys.
        <Param>
            antenna: antenna name code
            whatdat: 1 - raw data
                    2 - the basic interpolation
                    3 - interpolation with a specific cutoff or used 'sepif'
                    4 - smoothed data
        output: ~Tsys4.dat & ~Tsys4.png
        """
        self.logger.log("\n@task: cpif")
        self.logger.log("===========================")
        self.logger.log("*Find a png file of the resultant plot in the current path.")
        plt.rcParams['legend.frameon'] = 'False'
        plt.rcParams['xtick.labelsize'] = 11
        plt.rcParams['ytick.labelsize'] = 11
        for i,j in enumerate(self.legend, 1):
            self.logger.log(i," --> " ,j)
        colfrom = input("\n<Q> Copy-and-Paste from (see above and select a number): ")
        colfrom = int(colfrom) -1
        colto1 = input("\n<Q> Copy-and-Paste to (e.g., single IF: 3 / multiple IFs: 1,2,3): ")
        colto = []
        if len(colto1) == 1:
            colto.append( int(colto1) -1 )
        else:
            for i in colto1.split(','):
                colto.append( int(i) -1 )
        fig_pol=plt.figure(figsize=(12,9))
        antindex = self.antcodetrue.index(antenna)
        timearr = self.tarr[antindex]
        dali = []
        for i in range(len(self.legend)):
            k1 = self.legend[i]
            if i in colto:
                refdat = colfrom
                if whatdat == 1:
                    s2 = self.rlcp[antindex][refdat]
                    self.logger.log("\n *Raw data*")
                elif whatdat == 2:
                    s2 = np.genfromtxt(self.intpout[antindex][2:-1], usecols=refdat+2)
                    self.logger.log("\n *Data from the Basic interpolation*")
                elif whatdat == 3:
                    s2 = np.genfromtxt(self.intp[0][2:-1], usecols=refdat+2)
                    self.logger.log("\n *Data from another interpolation*")
                elif whatdat == 4:
                    s2 = np.genfromtxt(self.smip[2:-1], usecols=refdat+2)
                    self.logger.log("\n *Data from the Smoothing approach*")
            else:
                if whatdat == 1:
                    s2 = self.rlcp[antindex][i]
                    self.logger.log("\n *Raw data*")
                elif whatdat == 2:
                    s2 = np.genfromtxt(self.intpout[antindex][2:-1], usecols=i+2)
                    self.logger.log("\n *Data from the Basic interpolation*")
                elif whatdat == 3:
                    s2 = np.genfromtxt(self.intp[0][2:-1], usecols=i+2)
                    self.logger.log("\n *Data from another interpolation*")
                elif whatdat == 4:
                    s2 = np.genfromtxt(self.smip[2:-1], usecols=i+2)
                    self.logger.log("\n *Data from the Smoothing approach*")
            if self.columnorder.startswith('EU'):
                if i < int( len(self.legend)/2 ):
                    ccc = '#1f77b4'
                else:
                    ccc = '#ff7f0e'
            elif self.columnorder.startswith('VLBA'):
                if (i % 2) == 0:
                    ccc = '#1f77b4'
                else:
                    ccc = '#ff7f0e'
            if len(self.legend) == 16:
                plt.subplot(4, 4, i+1)
            elif len(self.legend) == 8:
                plt.subplot(4, 2, i+1)

            if i in colto:
                plt.plot(timearr, s2, markersize=1.5, marker='o', linewidth=0.8, alpha=1.0, zorder=0, color='r', label=k1+str("   *Copied*"))
            elif i == colfrom:
                plt.plot(timearr, s2, markersize=1.5, marker='o', linewidth=0.8, alpha=1.0, zorder=0, color='m', label=k1+str("   *FROM*"))
            else:
                plt.plot(timearr, s2, markersize=1.5, marker='o', linewidth=0.8, alpha=1.0, zorder=0, color=ccc, label=k1)
            self.logger.log(" ",k1)
            self.logger.log("---------------------------------------------")
            dali.append(s2)
            #% basic statistical properties
            vax = np.round(max(s2), 1)
            vin = np.round(min(s2), 1)
            vav = np.round(np.median(s2), 1)
            vst = np.round(np.std(s2), 1)
            #% plot layout
            fts3 = 11
            plt.legend(loc=2, fontsize=fts3, ncol=2, columnspacing=0.5, markerscale=1.2, labelspacing=0.3, handletextpad=0.3, borderaxespad=0.2, handlelength=1.5, handleheight=0.5)
            opz1 = np.linspace(plt.ylim()[0], plt.ylim()[1], 7)
            thegap = np.abs(np.abs(opz1[-1])-np.abs(opz1[-2]))
            plt.ylim(opz1[0]-thegap*0.0, opz1[-1]+thegap*1.3)
            plt.title("*"+antenna+"* : Max={}    Min={}    Med={}    Std={}".format(vax,vin,vav,vst), color='r', fontsize=12)
            plt.xlabel("Decimal day", fontsize=13)
            plt.ylabel("Tsys [K]", fontsize=13)
            plt.minorticks_on()
            plt.tick_params(axis='both', which='major', length=6.5, direction='in', pad=2, color='k')
            plt.tick_params(axis='both', which='minor', length=3, direction='in', pad=2, color='k')
        plt.show()
        plt.tight_layout()
        plt.savefig(self.logger.pref + antenna + "_" + self.session + "_" + self.freq1 +"Tsys4.png", format='png', dpi=150, transparent=False)
        #% make a new list with the final results
        td,tmh = np.genfromtxt(self.tsys1[antindex], usecols=(0,1), unpack=True, invalid_raise=False, missing_values='', usemask=False, dtype=str)
        finl = []
        for n1 in range(len(td)):
            fin99 = str(td[n1]) + " " + str(tmh[n1])
            for n2 in dali:
                fin99 = fin99 + " " + str(n2[n1])
                if len(fin99.split()) == len(self.legend) + 2:
                    fin99 = fin99 + " \n"
                else:
                    continue
            finl.append(fin99)
            del fin99
        finl.append("/\n")
        finl.insert(0, self.index[antindex])
        finl.insert(0, self.gain[antindex])
        self.cpdat = copy.deepcopy(finl)
        outn = self.logger.pref + antenna + "_" + self.session + "_" + self.freq1 +"Tsys4.dat"
        with open(outn, 'w') as ft57:
            ft57.writelines(finl)
        self.logger.log("*--> Use the output 'xxx_Tsys4.dat'. That's it!")
        self.logger.log("\nOUTPUT: Tsys4.dat file, self.cpdat, and a png file (i.e., xxxx_Tsys4.png)")





    def vlbacon(self):
        """
        An extra task only for VLBA.
            --A single VLBA antab file contains the data of multiple antennas.
        After the data examination, do this to connect all sessions (e.g., c211a, c211b, ..c, ..d)

        prefx = prefix of the output file name, use the session code (e.g., 2021A --> 'c211')
        frefx = observing frequency of the data (e.g., '3mm')
        sesfx = id of observing sessions (e.g., 'abc' for three sessions with ~a,~b,~c files
        OUTPUT: An antab file in a new form - Gain first, and then Tsys measurements.
        """
        prefx = input("<Q> Give a prefix of the output filename (e.g., 'c222') \n: ")
        frefx = input("<Q> Observing wavelength? (e.g., '3mm') \n: ")
        sesfx = input("<Q> Observing sessions? (e.g., 'abc' for three sessions; ~a,~b,~c.antabs) \n: ")
        con1 = ['BR','FD','HN','KP','LA','MK','NL','OV','PT','SC', 'GB']
        con2 = ['01','02','03','04','05','06','07','08','09','10', '11']
        for i in con1:
            con3 = 0
            con4 = []
            self.logger.log("\nSearching for..  **%s**" % i)
            for j in listdir():
                if (i in j) and ('png' not in j) and ('dat' in j):
                    self.logger.log(" ", j)
                    con3 = 1
                else:
                    continue
            if con3 == 0:
                self.logger.log("------------------------")
                self.logger.log("No data was found for", i)
                con4 = input("\n<Q> Press 'Enter' key if this is correct, OR say no \n:")
                if con4 == '':
                    continue
                else:
                    self.logger.log("\nSomething went wrong!")
                    return
            else:
                self.logger.log("----------------------")
                self.logger.log("Check the above output")
                con4 = input("\n<Q> List the file names in alphabetical order \n (e.g., C211a_BR_Tsys1.dat C211b_BR_Tsys1.dat C211c_BR_Tsys1.dat ..)   **Space btw. files \n Use '~Tsys1.dat' when both '~Tsys1.dat' and '~Tsys0.dat' exist for a telescope! \n: ")
            con5 = con4.split()
            if len(con5) == 1:
                self.logger.log("\n-----!!!Only one session, just use it!!!-----"); time.sleep(2)
                continue
            else:
                COMB = []
                for k in con5:
                    if k == con5[0]:
                        GG1 = open(k); GG2 = GG1.readlines()
                        COMB = COMB + GG2[:-1]
                        del GG1; del GG2
                        continue
                    elif (k != con5[0]) and (k != con5[-1]):
                        GG1 = open(k); GG2 = GG1.readlines()
                        COMB = COMB + GG2[2:-1]
                        del GG1; del GG2
                        continue
                    elif k == con5[-1]:
                        GG1 = open(k); GG2 = GG1.readlines()
                        COMB = COMB + GG2[2:]
                        del GG1; del GG2
        #      filnam = 'FIN_' + prefx +"_"+ i +"_"+ frefx +"Tsys.dat"
            filnam = "MERG_" + prefx + 'VLBA' + "-" + i + "_" + sesfx + "_" + frefx + "Tsys.dat"
            with open(filnam, 'w') as ft44:
                ft44.writelines(COMB)
        self.logger.log("\n--FINISHED!--")


class TypeB:
    """
    MULTIPLE TSYS COLUMNS: suitable for European antennas e.g., EF, ON, MH, ...
    output: self.dall, self.tsys, self.arg, self.antcode, self.gain, self.index, self.totcol, and self.legend
    """
    def __init__(self, antcode, antfile, session_codes, interactive_mode, wavelength):
        self.wavelength= wavelength
        self.interactive_mode = interactive_mode
        self.antcode = antcode
        self.antfile = antfile
        self.session = session_codes
        self.logger = Logger()
        self.logger.genlog()
        self.tsysread()
        self.aa = None    # Initialize aa as None
        self.tarr= None
        self.columnorder =None 



    def parser(self):
        if self.interactive_mode:
            arg = []
            reciv = input("\n<Q> Give antab file(s) of this station/session \n (e.g., 'c211aef.antabfs c211bef.antabfs c211cef.antabfs' --> space between files) \n: ")
            arg = reciv.split()
            self.session = input("\n<Q> Give all the session codes of the input data? (e.g., 'abc' for three a,b,c antab files) \n: ")
            self.logger.log("\n Total %.0f sessions \n" % len(arg))
        else:
            arg = self.antfile
            # reciv = self.antfile
            # arg = reciv.split()
            # self.session =[]
            # for filename in reciv.split():
            #     pattern = r'^\w{4}(\w)'                
            #     match = re.search(pattern, filename)
            #     if match:
            #         self.session.append(match.group(1))
            # self.session = ','.join(self.session)
            # self.logger.log("\n since you selected " + reciv + " , your session code(s) is (are) probably: \n" + self.session)

        self.arg = arg
        N = len(arg)
        dset = []
        for i in range(N):
            ts1 = open(arg[i])
            ts2 = ts1.readlines()
            dset.append(ts2)

        dat = dset
        if len(dat) == 1:
            datfull = dat[0]
        else:
            datfull = dat[0]
            for i in range(len(dat)):
                if i == 0:
                    continue
                else:
                    datfull = datfull + dat[i]

        for i, j in enumerate(datfull[:30]):
            print(str(i) + ")", j)
        ################################################################        
        gain_line = ''
        poly_line = ''
        index_line = ''
        combined_data = '\n'.join(datfull)

        gain_match = re.search(r'GAIN[^/]+POLY.*?/', combined_data, re.DOTALL)
        if gain_match:
            gain_line = gain_match.group(0).split('POLY')[0].strip()

        poly_match = re.search(r'POLY=.*?/', combined_data, re.DOTALL)
        if poly_match:
            poly_line = poly_match.group(0).strip().replace('\n', '')

        index_match = re.search(r'TSYS\s+[^/]+INDEX\s*=\s*[^/]+/', combined_data, re.IGNORECASE)
        if index_match:
            index_line = index_match.group(0).strip().replace('\n', '')
        else:
            print("INDEX line not available, please check this information")

        ################################################################

        self.logger.log("\n<Q> Find the *GAIN LINE* \n(e.g., GAIN AA ELEV DPFU=1.0) \n (e.g., GAIN EF ELEV DPFU=0.14,0.14 FREQ=84000,95500) \n **FREQ is optional and you can skip it \n**Do NOT end this line with '/'\n")
        if self.interactive_mode:
            while True:
                provided_gain = input("\nIs this the right GAIN line? Do not forget the station code!\n(Press ENTER to confirm or input a different line): \n" + gain_line + "\n: ")
                if provided_gain.strip() == "":
                    gaininfo = gain_line
                    break
                else:
                    gaininfo = provided_gain
                    if not gaininfo.startswith('GAIN'):
                        print("Please provide a line starting with 'GAIN'.")
                    else:
                        break
        else:
            if self.antcode in ['KU', 'KT', 'KY']:
                if self.antcode == 'KU' and self.wavelength =='3mm':
                    gaininfo = "GAIN KU ELEV FREQ=84000,117000 DPFU=0.04842397"
                if self.antcode == 'KT' and self.wavelength =='3mm':
                    gaininfo = "GAIN KT ELEV FREQ=84000,117000 DPFU=0.06290550"
                if self.antcode == 'KY' and self.wavelength =='3mm':
                    gaininfo = "GAIN KY ELEV FREQ=84000,117000 DPFU=0.05205973"
            else:
                self.logger.log("\nGAIN line): \n" + gain_line + "\n: ")
                gaininfo = gain_line
        self.logger.log("\n<Q> Find the *POLY LINE* \n(e.g., POLY=1.0 /) \n(e.g., POLY=0.727089,0.947364E−02,-0.822152E−04 /) \n**End this line with '/' \n ")
        if self.interactive_mode:
            while True:
                provided_poly = input("\nIs this the right POLY line? \n(Press ENTER to confirm or input a different line): \n" + poly_line + "\n: ")
                if provided_poly.strip() == "":
                    gainpoly = poly_line
                    break
                else:
                    gainpoly = provided_poly
                    if not gainpoly.startswith('POLY'):
                        print("Please provide a line starting with 'POLY'.")
                    else:
                        break
        else:
            if self.antcode in ['KU', 'KT', 'KY']:
                if self.antcode == 'KU' and self.wavelength =='3mm':
                    gainpoly = "POLY=0.58343484,0.02216809,-0.00019488 /"
                if self.antcode == 'KT' and self.wavelength =='3mm':
                    gainpoly = "POLY=1.04279778,0.00545253,-0.00006587 /"
                if self.antcode == 'KY' and self.wavelength =='3mm':
                    gainpoly = "POLY=1.04279778,0.00545253,-0.00006587 /"
            self.logger.log("\nPOLY line: \n" + poly_line + "\n: ")
            gainpoly = poly_line
        self.logger.log("\n<Q> Check the number of Tsys columns and the R/L order! \n (e.g., TSYS EF INDEX='R1','R2','R3','R4','R5','R6','R7','R8','L1','L2','L3','L4','L5','L6','L7','L8' /) \n **End the line with '/'\n ")
        if self.interactive_mode:
            while True:
                provided_index = input("\nIs this the right TSYS INDEX line? \n(Press ENTER to confirm or input a different line): \n" + index_line + "\n: ")
                if provided_index.strip() == "":
                    indexinfo = index_line
                    break
                else:
                    indexinfo = provided_index
                    if not indexinfo.startswith('TSYS'):
                        print("Please provide a line starting with 'TSYS'.")
                    else:
                        break
        else:
            self.logger.log("\nTSYS INDEX line: \n" + index_line + "\n: ")
            indexinfo = index_line
        gaininfo = gaininfo + " " + gainpoly
        self.gain = gaininfo + '\n'
        self.index = indexinfo + '\n'
        self.datfull = datfull
        ################################################################
        # end parser for GAIN, POLY, INDEX
        ################################################################

        #now: parser for tsys columns
        if self.antcode in ['EF', 'ON', 'MH', 'YS', 'KT', 'KU', 'KY']:
            self.logger.log("\n<Q> How many tsys-columns in the data? (give a number e.g., 16 or 8) \n ")
            if self.interactive_mode:
                hd_legnum = input("Since it's " + self.antcode + " there are probably 16 tsys-columns, (Press ENTER to confirm or give the right number)\n: ")
                if hd_legnum.strip() == "":
                    hd_legnum = 16
                    hd_legnum = int(hd_legnum)
                    self.totcol = hd_legnum
                    self.logger.log("\nTotal", hd_legnum, "columns of tsys measurements.")
                else:
                    hd_legnum = int(hd_legnum)
                    self.totcol = hd_legnum
                    self.logger.log("\nTotal", hd_legnum, "columns of tsys measurements.")
            else:
                self.logger.log("Since it's " + self.antcode + " there are probably 16 tsys-columns\n: ")
                hd_legnum = 16
                hd_legnum = int(hd_legnum)
                self.totcol = hd_legnum
                self.logger.log("\nTotal", hd_legnum, "columns of tsys measurements.")
        else:
            hd_legnum = input("\n<Q> How many tsys-columns in the data? (give a number e.g., 16 or 8) \n: ")
            hd_legnum = int(hd_legnum)
            self.totcol = hd_legnum
            self.logger.log("\nTotal", hd_legnum, "columns of tsys measurements.")

        #now: parser for legend
        lgst = []
        self.logger.log("\n<Q> For plot-legends later, see the following and select your data order by giving its number \n 1 -  R1,R2,R3,R4,L1,L2,L3,L4 \n 2 -  R1,R2,R3,R4,R5,R6,R7,R8,L1,L2,L3,L4,L5,L6,L7,L8  --> (usually EU ants.) \n 3 -  R1,R,L,R,L,R,L,R,L                         ---> (usually VLBA ants.) \n 4 -  R1,L1,R2,L2,R3,L3,R4,L4,R5,L5,R6,L6,R7,L7,R8,L8 \n 5 -  RCP,LCP                                      ------> (usually those Special ants.) \n (..not there? then Give it manually here: e.g., R1,L1,L2, ...)\n")
        if self.interactive_mode:
            while True:
                if self.antcode in ['EF', 'ON', 'MH', 'YS']:
                    provided_indleg = input("Since it's " + self.antcode + " the data order is probably \n2 -  R1,R2,R3,R4,R5,R6,R7,R8,L1,L2,L3,L4,L5,L6,L7,L8, \n(Press ENTER to confirm or give the right number)\n:  ")
                    if provided_indleg.strip() == "":
                        legtem = 'R1,R2,R3,R4,R5,R6,R7,R8,L1,L2,L3,L4,L5,L6,L7,L8'
                        lgst = legtem.split(',')
                        break
                else:
                    legtem = input(": ")
                    if len(legtem) == 1:
                        legtem = int(legtem)
                        if legtem == 1:
                            legtem = 'R1,R2,R3,R4,L1,L2,L3,L4'
                        elif legtem == 2:
                            legtem = 'R1,R2,R3,R4,R5,R6,R7,R8,L1,L2,L3,L4,L5,L6,L7,L8'
                        elif legtem == 3:
                            legtem = 'R1,R,L,R,L,R,L,R,L'
                        elif legtem == 4:
                            legtem = 'R1,L1,R2,L2,R3,L3,R4,L4,R5,L5,R6,L6,R7,L7,R8,L8'
                        elif legtem == 5:
                            legtem = 'RCP,LCP'
                            lgst = legtem.split(',')
                            break
        else:
            self.logger.log("Since it's " + self.antcode + " the data order is \n2 -  R1,R2,R3,R4,R5,R6,R7,R8,L1,L2,L3,L4,L5,L6,L7,L8)\n  ")
            legtem = 'R1,R2,R3,R4,R5,R6,R7,R8,L1,L2,L3,L4,L5,L6,L7,L8'
            lgst = legtem.split(',')
        self.legend = lgst


    def tsysread(self):
        print("===========================")
        print('Starting tsysread for '+ self.antcode)
        print("===========================")
        self.logger.log("\n@Task: tsysread")

        print(" ")
        print("'ufo()' has some useful information!")
        # self.freq = input("<Q> Observing wavelength? (e.g., 3mm) \n:")
        # self.logger.log("\nAn observing wavelength of"+ self.wavelength +"is assumed\n")
        # self.freq = self.wavelength

        if self.interactive_mode:
            self.freq = input("<Q> Observing wavelength? (e.g., 3mm) \n:")
            
        else:
            self.logger.log("\nAn observing wavelength of "+ self.wavelength +" is assumed\n")
            self.freq= self.wavelength
        for i in os.listdir():
            if i.endswith('.antab') or i.endswith('.antabfs'):
                self.logger.log(i)

        #start parser
        self.parser()

        datfull = self.datfull 
        self.dall = datfull
        temdat = copy.deepcopy(datfull)
        for j in temdat[:]:
            if j[0].isdigit():
                continue
            else:
                temdat.remove(j)

        self.tsys = temdat      # only actual data for plotting
        self.logger.log("\n***How it looks?")
        for zz in range(15):
            self.logger.log(temdat[zz])

        self.logger.log("\nOUTPUT: self.dall, self.tsys, self.arg, self.antcode, self.freq, self.gain, self.index, self.totcol, self.session, and self.legend")
        self.logger.log(" **just-in-case**..if Something went wrong (incorrect inputs among the above Pars.)")
        self.logger.log("   ,one can give it manually (e.g., self.xx = !!) or repeat this function again.")


    def dbcheck(self):
        """
        To check if the data contains unreliable values.
        run it like e.g., self.dbcheck()
        output: self.tarr
        """
        self.logger.log("===========================")
        self.logger.log("\n@Task: dbcheck")
        self.logger.log("===========================")
        new_t = []
        simpletest = self.tsys[0].split()[1]
        if simpletest.count(":") == 2:
            for l in self.tsys:
                tem_d = l.split()[0]        #% day array
                tem_hms = l.split()[1]       #% hh:mm:ss array
                #% converting all of them to a day with decimal numbers
                div_h = tem_hms.split(':')[0]
                div_m = tem_hms.split(':')[1]
                div_s = tem_hms.split(':')[2]
                dayform = (int(div_h)/24.) + (float(div_m)/60. /24.) + (float(div_s)/60. /60. /24.) + int(tem_d)
                new_t.append(dayform)

        elif simpletest.count(":") == 1:
            for l in self.tsys:
                tem_d = l.split()[0]        #% day array
                tem_hm = l.split()[1]       #% hh:mm array
                div_h = tem_hm.split(':')[0]
                div_m = tem_hm.split(':')[1]
                dayform = (int(div_h)/24.) + (float(div_m)/60. /24.) + int(tem_d)
                new_t.append(dayform)
                #% --> actual x-array for interpolation

        new_t = np.asarray(new_t)
        self.tarr = new_t
    #      log("\nEuropean antennas usually have 8 rcp columns and 8 lcp columns, respectively. \n")

        if self.interactive_mode:
            sep1 = input("\n<Q> Give the order of RCP/LCP in the Tsys columns from left to right, either 1 or 2? \n  1 (EU ant. form: = R1 R2 R3 R4 ... L1 L2 L3 L4 ...) \n  2 (VLBA ant. form: = R1 L1 R2 L2 R3 L3 R4 L4 ...) \n: ")
            sep1 = int(sep1)
            if sep1 == 1:
                self.columnorder = "EU normally --> R1 R2 R3 R4 ... L1 L2 L3 L4 ..."
            elif sep1 == 2:
                self.columnorder = "VLBA normally --> R1 L1 R2 L2 R3 L3 R4 L4 ..."
        else:
            self.logger.log("\nOrder of RCP/LCP in the Tsys columns from left to right is (EU ant. form: = R1 R2 R3 R4 ... L1 L2 L3 L4 ...) ")
            sep1 = 1
            self.columnorder = "EU normally --> R1 R2 R3 R4 ... L1 L2 L3 L4 ..."
            
        tsyslist = []
        for kk1 in range(len(self.legend)):
            eachcol = np.genfromtxt(self.tsys, usecols=kk1+2, unpack=True, invalid_raise=False, missing_values='', usemask=False, dtype=float)
            tsyslist.append(eachcol)
            del eachcol

        self.rlcp = tsyslist
        self.logger.log("\n**Now, checking any values of [999, >9999, 0, negative, btw 0 and 10]")
        aa = 0
        if sep1 == 1:       #% first half -> rcp, and then lcp (e.g., EU antennas)
            for i in range(len(tsyslist)):
                if i < int( len(tsyslist)/2 ):
                    s2 = tsyslist[i]
                    self.logger.log("..First, RCP !")
                    for k1,k2 in enumerate(s2,0):
                        if (k2 == 999) or (k2 == 999.0):
                            self.logger.log("-- 999 found!!")
                            aa = 1
                            continue
                        elif k2 >= 9999:
                            self.logger.log("-- higher than 9999 found!!")
                            aa = 1
                            continue
                        elif k2 == 0:
                            self.logger.log("-- 0 found!!")
                            aa = 1
                            continue
                        elif k2 < 0:
                            self.logger.log("-- negative value found!!")
                            aa = 1
                            continue
                        elif 0 < k2 < 10:
                            self.logger.log("-- 0K < tsys < 10K, unreasonably low value found!!")
                            aa = 1
                            continue
                else:
                    s2 = tsyslist[i]
                    self.logger.log("..Second, LCP !")
                    for k1,k2 in enumerate(s2,0):
                        if (k2 == 999) or (k2 == 999.0):
                            self.logger.log("-- 999 found!!")
                            aa = 1
                            continue
                        elif k2 >= 9999:
                            self.logger.log("-- higher than 9999 found!!")
                            aa = 1
                            continue
                        elif k2 == 0:
                            self.logger.log("-- 0 found!!")
                            aa = 1
                            continue
                        elif k2 < 0:
                            self.logger.log("-- negative value found!!")
                            aa = 1
                            continue
                        elif 0 < k2 < 10:
                            self.logger.log("-- 0K < tsys < 10K, unreasonably low value found!!")
                            aa = 1
                            continue
        elif sep1 == 2:       #% rcp -> even, and lcp -> odd (e.g., VLBA)
            for i in range(len(tsyslist)):
                if (i % 2) == 0:     #% even, thus RCP
                    s2 = tsyslist[i]
                    self.logger.log("..First, RCP !")
                    for k1,k2 in enumerate(s2,0):
                        if (k2 == 999) or (k2 == 999.0):
                            self.logger.log("-- 999 found!!")
                            aa = 1
                            continue
                        elif k2 >= 9999:
                            self.logger.log("-- higher than 9999 found!!")
                            aa = 1
                            continue
                        elif k2 == 0:
                            self.logger.log("-- 0 found!!")
                            aa = 1
                            continue
                        elif k2 < 0:
                            self.logger.log("-- negative value found!!")
                            aa = 1
                            continue
                        elif 0 < k2 < 10:
                            self.logger.log("-- 0K < tsys < 10K, unreasonably low value found!!")
                            aa = 1
                            continue
                else:     #% odd, thus LCP
                    s2 = tsyslist[i]
                    self.logger.log("..Second, LCP !")
                    for k1,k2 in enumerate(s2,0):
                        if (k2 == 999) or (k2 == 999.0):
                            self.logger.log("-- 999 found!!")
                            aa = 1
                            continue
                        elif k2 >= 9999:
                            self.logger.log("-- higher than 9999 found!!")
                            aa = 1
                            continue
                        elif k2 == 0:
                            self.logger.log("-- 0 found!!")
                            aa = 1
                            continue
                        elif k2 < 0:
                            self.logger.log("-- negative value found!!")
                            aa = 1
                            continue
                        elif 0 < k2 < 10:
                            self.logger.log("-- 0K < tsys < 10K, unreasonably low value found!!")
                            aa = 1
                            continue

        self.showsys() # ADDED HERE     
        if aa == 0:
            self.logger.log("\n---------------------------------------------")
            self.logger.log("***All measurements seem FINE! Ready to be used in the calibration!***")
            self.logger.log(" -Find an output file: xx_Tsys0.dat, in your working directory.")
            self.logger.log(" -But, STRONGLY RECOMMEND to plot for double-checking the data with showsys().")
            self.tsysout = copy.deepcopy(self.tsys)
            self.tsysout.append("/\n")
            self.tsysout.insert(0, self.index)
            self.tsysout.insert(0, self.gain)
            outn = self.logger.pref +"_"+self.antcode + "_"+ "_".join(self.session)+ "_" + self.freq+ "_"+"Tsys0.dat"
            with open(outn, 'w') as ft26:
                ft26.writelines(self.tsysout)
            self.logger.log(" -If the plot looks ok, then, use the output file: xx_Tsys0.dat, here.")
        elif aa == 1:
            self.logger.log("\n---------------------------------------------")
            self.logger.log("***Unreliable value found!***")
            self.logger.log(" -Find an output file: xx_Tsys0.dat, in your working directory.")
            self.logger.log(" -Fix it further with showsys() and intpsys()")
            self.logger.log(" -Then, use self.intp with writetxt()")
            self.tsysout = copy.deepcopy(self.tsys)
            self.tsysout.append("/\n")
            self.tsysout.insert(0, self.index)
            self.tsysout.insert(0, self.gain)
            outn = self.logger.pref+"_" +self.antcode+"_"+ "_".join(self.session) + "_" + self.freq+ "_"+"Tsys0.dat"
            with open(outn, 'w') as ft26:
                ft26.writelines(self.tsysout)
 
            self.intpsys()
        self.logger.log("\nOUTPUT: self.tarr, self.columnorder, self.tsysout, and self.rlcp")

    def showsys(self, Ymax=None, Ymin=None):

        """
        Plot tsys measurements
        Ymax: Y-axis range, Maximum value.
        Ymin: Y-axis range, Minimum value.
        """
        self.logger.log("===========================")
        self.logger.log("\n@Task: showsys")
        self.logger.log("===========================")
        new_t = self.tarr
        self.logger.log("*Find a png file in your working directory")
        fig_pol=plt.figure(figsize=(17,10.5))
        plt.rcParams['legend.frameon'] = 'False'
        plt.rcParams['xtick.labelsize'] = 11
        plt.rcParams['ytick.labelsize'] = 11
        for i in range(len(self.legend)):
            k1 = self.legend[i]
            k2 = self.rlcp[i]
            if self.columnorder.startswith('EU'):
                if i < int( len(self.legend)/2 ):
                    ccc = '#1f77b4'
                else:
                    ccc = '#ff7f0e'
            elif self.columnorder.startswith('VLBA'):
                if (i % 2) == 0:
                    ccc = '#1f77b4'
                else:
                    ccc = '#ff7f0e'

            if len(self.legend) == 16:
                plt.subplot(4, 4, i+1)
            elif len(self.legend) == 8:
                plt.subplot(4, 2, i+1)
            plt.plot(new_t, k2, markersize=1.5, marker='o', linewidth=0.8, label=k1, color=ccc)
            fts2 = 9.5
            plt.legend(loc=2, fontsize=fts2)
            vax = int(np.round(max(k2)))
            vin = int(np.round(min(k2)))
            vav = int(np.round(np.median(k2)))
            vst = int(np.round(np.std(k2)))
            fts1 = 10
            if Ymax:
                if Ymin:
                    plt.ylim([Ymin, Ymax])
                else:
                    plt.ylim([0, Ymax])
            else:
                opz1 = np.linspace(plt.ylim()[0], plt.ylim()[1], 7)
                thegap = np.abs(np.abs(opz1[-1])-np.abs(opz1[-2]))
                plt.ylim(opz1[0]-thegap*0.0, opz1[-1]+thegap*1.3)
            plt.title("*"+self.antcode+"*: Max={}, Min={}, Med={}, Std={}".format(vax,vin,vav,vst), color='r', fontsize=fts1)
            plt.xlabel("Decimal day")
            plt.ylabel("Tsys [K]")
            plt.minorticks_on()
            plt.tick_params(axis='both', which='major', length=5, direction='in', pad=2, color='k')
            plt.tick_params(axis='both', which='minor', length=3, direction='in', pad=2, color='k')
        plt.tight_layout()
        plt.subplots_adjust(left=0.045, bottom=0.050, right=0.975, top=0.965, wspace=0.190, hspace=0.345)
        if Ymax:
            self.logger.log("\nNo output .png file.")
        else:
            plt.savefig(self.logger.pref+"_"+self.antcode+"_" + "_".join(self.session) + "_" + self.freq+ "_"+"Tsys0.png", format='png', dpi=150, transparent=False)
        # plt.show()
        self.logger.log("\n*If the plot looks OK and dbcheck did not report suspicious tsys values,")
        self.logger.log("*--> just use the output file of dbcheck. That's it!")
        self.logger.log("\nOUTPUT: A figure file (i.e., xxxx_Tsys0.png)")

    def intpsys(self, hnot=9999, lnot=10, sepif=False):
        """
        Linear interpolation to replace those unreliable tsys values.
        One can modify the cutoff of the high and low Tsys values.
        """
        print("===========================")
        print('Starting intpsys')
        print("===========================")
        self.logger.log("\n@task: intpsys")
        self.logger.log("===========================")
        timearr = self.tarr
        self.logger.log("*Find a png file of the resultant plot.")
        fig_pol=plt.figure(figsize=(17,10.5))
        plt.rcParams['legend.frameon'] = 'False'
        plt.rcParams['xtick.labelsize'] = 11
        plt.rcParams['ytick.labelsize'] = 11
        dali = []
    #      statrl = []
        for i in range(len(self.legend)):
            k1 = self.legend[i]
            s2 = self.rlcp[i]
            if self.columnorder.startswith('EU'):
                if i < int( len(self.legend)/2 ):
                    ccc = '#1f77b4'
                else:
                    ccc = '#ff7f0e'
            elif self.columnorder.startswith('VLBA'):
                if (i % 2) == 0:
                    ccc = '#1f77b4'
                else:
                    ccc = '#ff7f0e'

            if len(self.legend) == 16:
                plt.subplot(4, 4, i+1)
            elif len(self.legend) == 8:
                plt.subplot(4, 2, i+1)

            #% linear interpolation
            temp1 = copy.deepcopy(timearr)
            temp2 = copy.deepcopy(s2)
            wrong_index = []
            wrong_times = []
            if sepif:
                self.logger.log("\n**Set high/low Tsys cutoff for each RCP/LCP IF**")
                self.logger.log("\n-->", k1)
                hgo = input("<Q> Tsys shouldn't EXCEED (press 'Enter' to skip): ")
                if len(hgo) == 0:
                    hgo = hnot
                else:
                    hgo = float(hgo)
                lgo = input("<Q> Tsys shouldn't be BELOW (press 'Enter' to skip): ")
                if len(lgo) == 0:
                    lgo = lnot
                else:
                    lgo = float(lgo)
                self.logger.log("\n!! any tsys values of 999, >=YOUR CHOICE, <0, 0, and 0<tsys<YOUR CHOICE will be removed !!")
                for aa,bb in enumerate(temp2, 0):
                    if (bb == 999) or (bb == 999.0):
                        wrong_index.append(aa)
                        wrong_times.append(temp1[aa])
                        continue
                    elif bb >= hgo:
                        wrong_index.append(aa)
                        wrong_times.append(temp1[aa])
                        continue
                    elif bb < 0:
                        wrong_index.append(aa)
                        wrong_times.append(temp1[aa])
                        continue
                    elif bb == 0:
                        wrong_index.append(aa)
                        wrong_times.append(temp1[aa])
                        continue
                    elif 0 < bb < lgo:
                        wrong_index.append(aa)
                        wrong_times.append(temp1[aa])
                        continue
            else:
                self.logger.log("\n!! any tsys values of 999, >=9999, <0, 0, and 0<tsys<10 will be removed !!")
                for aa,bb in enumerate(temp2, 0):
                    if (bb == 999) or (bb == 999.0):
                        wrong_index.append(aa)
                        wrong_times.append(temp1[aa])
                        continue
                    elif bb >= hnot:
                        wrong_index.append(aa)
                        wrong_times.append(temp1[aa])
                        continue
                    elif bb < 0:
                        wrong_index.append(aa)
                        wrong_times.append(temp1[aa])
                        continue
                    elif bb == 0:
                        wrong_index.append(aa)
                        wrong_times.append(temp1[aa])
                        continue
                    elif 0 < bb < lnot:
                        wrong_index.append(aa)
                        wrong_times.append(temp1[aa])
                        continue

            #% x & y array to be used for interpolation
            temp3 = np.asarray([tg for ii, tg in enumerate(temp1) if ii not in wrong_index])
            temp4 = np.asarray([tg for ii, tg in enumerate(temp2) if ii not in wrong_index])
            #=================================
            #% in case the first and/or last measurements are not reasonable!!
            t1 = list(copy.deepcopy(timearr)); t2 = list(copy.deepcopy(temp3))
            z1 = list(copy.deepcopy(s2)); z2 = list(copy.deepcopy(temp4))
            if temp3[0] != timearr[0]:
                self.logger.log("tsys starts with wrong values.. to fill in: using Mean of 15% of the data starting from left to right!")
                fil1 = t1.index(t2[0])
                avg_front = np.mean(temp4[:round(len(t1)*0.15+0.5)])
                for j in range(fil1):
                    z2.insert(0, avg_front)
                    t2.insert(j, t1[j])

            else:
                self.logger.log("Front: good to go (the first data point looks ok)")

            if temp3[-1] != timearr[-1]:
                self.logger.log("tsys ends with wrong values.. to fill in: using Mean of 15% of the data starting from right to left!")
                for k in np.arange(-1, np.negative(len(t1)), -1):
                    if t2[-1] == t1[k]:
                        fil2 = k
                        break

                avg_end = np.mean(temp4[-round(len(t1)*0.15+0.5):])
                for m in np.arange(fil2+1, 0, +1):
                    z2.append(avg_end)
                    t2.append(t1[m])

            else:
                self.logger.log("End: good to go (the last data point looks ok)")
            #=================================
            #% Linear interpolation
            gtg_time = np.asarray(t2)
            gtg_tsys = np.asarray(z2)
            linv = interp1d(gtg_time, gtg_tsys, kind='linear')

            plt.plot(timearr, linv(timearr), markersize=1.5, marker='o', linewidth=0.8, alpha=1.0, zorder=0, color=ccc, label=k1)
            if wrong_index:
                plt.plot(wrong_times, linv(wrong_times), linewidth=0, marker='*', markersize=9, color='#2ca02c', label='Lin. Interp.')
                self.logger.log("\nInterpolation has been performed!!")
            else:
                self.logger.log("\nAll measurements are reasonable!!")

            self.logger.log(" ",k1)
            self.logger.log("  --> Number of unreasonable tsys values:", len(wrong_index))
            self.logger.log("---------------------------------------------")
            resu = np.asarray([round(g, 1) for g in linv(timearr)])   #% for recording later..
            dali.append(resu); del(resu)
            #% basic statistical properties
            vax = int(np.round(max(linv(timearr))))
            vin = int(np.round(min(linv(timearr))))
            vav = int(np.round(np.median(linv(timearr))))
            vst = int(np.round(np.std(linv(timearr))))
            #% plot layout
            fts3 = 9.5
            plt.legend(loc=2, fontsize=fts3, ncol=1, columnspacing=0.9, markerscale=1.2)
            fts1 = 10
            opz1 = np.linspace(plt.ylim()[0], plt.ylim()[1], 7)
            thegap = np.abs(np.abs(opz1[-1])-np.abs(opz1[-2]))
            plt.ylim(opz1[0]-thegap*0.0, opz1[-1]+thegap*1.3)
            plt.title("*"+self.antcode+"*: Max={}, Min={}, Med={}, Std={}".format(vax,vin,vav,vst), color='r', fontsize=fts1)
            plt.xlabel("Decimal day")
            plt.ylabel("Tsys [K]")
            plt.minorticks_on()
            plt.tick_params(axis='both', which='major', length=5, direction='in', pad=2, color='k')
            plt.tick_params(axis='both', which='minor', length=3, direction='in', pad=2, color='k')
        # plt.tight_layout()
        # plt.show()
        plt.subplots_adjust(left=0.045, bottom=0.050, right=0.975, top=0.965, wspace=0.190, hspace=0.345)
        for blk in range(4):
            self.logger.log("     ")

        self.logger.log(">>>> NEXT ANT!")
        self.logger.log(">>>> Resultant plot is shown (outliers are now replaced with interpolared values)")
        if (hnot != 9999) or (sepif == True):
            plt.savefig(self.logger.pref+"_"+self.antcode+"_"+self.session+"_" + self.freq+ "_"+"Tsys3.png", format='png', dpi=150, transparent=False)
        else:
            plt.savefig(self.logger.pref+"_"+self.antcode+"_"+"_".join(self.session)+"_" + self.freq+ "_"+"Tsys1.png", format='png', dpi=150, transparent=False)
        #% make a new list with the final results
        td,tmh = np.genfromtxt(self.tsys, usecols=(0,1), unpack=True, invalid_raise=False, missing_values='', usemask=False, dtype=str)
        finl = []
        for n1 in range(len(td)):
            fin99 = str(td[n1]) + " " + str(tmh[n1])
            for n2 in dali:
                fin99 = fin99 + " " + str(n2[n1])
                if len(fin99.split()) == self.totcol + 2:
                    fin99 = fin99 + "  \n"
                continue
            finl.append(fin99)
            del fin99
            if n1 == len(td)-1:
                finl.append("/\n")
            else:
                continue
        finl.insert(0, self.index)
        finl.insert(0, self.gain)
        self.intp = copy.deepcopy(finl)
        if (hnot != 9999) or (sepif == True):
            outn = self.logger.pref+"_" +self.antcode+"_"+ "_".join(self.session) + "_" + self.freq+ "_"+"Tsys3.dat"
        else:
            #try:
            # outn = self.pref.prefLog +"_"+ self.session + "_" + self.freq+"Tsys1.dat"
            outn = self.logger.pref+"_" +self.antcode+"_"+ "_".join(self.session) + "_" + self.freq+ "_"+"Tsys1.dat"
        with open(outn, 'w') as ft24:
            ft24.writelines(self.intp)
        # plt.show()
        self.logger.log("\n*If the result looks OK with those replaced Tsys values,")
        self.logger.log("*--> Use the output 'xxx_Tsys1.dat'. That's it!")
        self.logger.log("\nOUTPUT: self.intp and a figure file (i.e., xxxx_Tsys1.png)")

    def interv(a1):     # checking the sampling interval of data.
        a=[]
        for i in np.arange(0, len(a1)):
            b = a1[i+1] - a1[i]
            a = np.append(a, b)
            if i == len(a1)-2:
                print(a)
                break
        return np.mean(a), a

    def smth(self, whatdat=None, siglev=2):

        """
        To smooth out the curve by removing OUTLIERS!
        Para. whatdat: 1 if raw data, 2 if interpolated data.
        Para. siglev: confidence level of the sampling interval for the cutoff; e.g., 2 --> 2*std(the intervals). Default=2.
        output: ~Tsys2.dat, ~Tsys2.png
        """
        self.logger.log("\n@task: smth")
        self.logger.log("===========================")
        timearr = self.tarr
        self.logger.log("*Find a png file of the resultant plot in the current path.")
        fig_pol=plt.figure(figsize=(17,10.5))
        plt.rcParams['legend.frameon'] = 'False'
        plt.rcParams['xtick.labelsize'] = 11
        plt.rcParams['ytick.labelsize'] = 11
        dali = []
    #      statrl = []
        for i in range(len(self.legend)):
            k1 = self.legend[i]
            if whatdat == 1:
                s2 = self.rlcp[i]
            elif whatdat == 2:
                s2 = np.genfromtxt(self.intp[2:-1], usecols=i+2)
            if self.columnorder.startswith('EU'):
                if i < int( len(self.legend)/2 ):
                    ccc = '#1f77b4'
                else:
                    ccc = '#ff7f0e'
            elif self.columnorder.startswith('VLBA'):
                if (i % 2) == 0:
                    ccc = '#1f77b4'
                else:
                    ccc = '#ff7f0e'

            if len(self.legend) == 16:
                plt.subplot(4, 4, i+1)
            elif len(self.legend) == 8:
                plt.subplot(4, 2, i+1)

            #% linear interpolation
            temp1 = copy.deepcopy(timearr)
            temp2 = copy.deepcopy(s2)
            print(temp2)
            print(self.interv(temp2)[1])
            cutC = siglev * np.std(self.interv(temp2)[1])
            wrong_index = []
            wrong_times = []
            self.logger.log("\n!! Gonna take the Outliers out from the data !!")
            for aa,bb in enumerate(temp2, 0):
                if aa == len(temp2)-1:
                    break
                if np.abs(temp2[aa+1]-bb) > cutC:
                    if aa < round(len(temp2)/2):
                        half1 = np.mean(temp2[aa:aa+6])  # get a local mean value and judge the data.
                        temarr9 = np.asarray([bb, temp2[aa+1]])
                        inddex = np.abs(temarr9 - half1).argmax()
                        if temarr9[inddex] == bb:
                            if aa in wrong_index:
                                continue
                            else:
                                wrong_index.append(aa)
                                wrong_times.append(temp1[aa])
                                continue
                        else:
                            if aa+1 in wrong_index:
                                continue
                            else:
                                wrong_index.append(aa+1)
                                wrong_times.append(temp1[aa+1])
                                continue
                else:
                    half2 = np.mean(temp2[aa-5:aa+1])
                    temarr9 = np.asarray([bb, temp2[aa+1]])
                    inddex = np.abs(temarr9 - half2).argmax()
                    if temarr9[inddex] == bb:
                        if aa in wrong_index:
                            continue
                        else:
                            wrong_index.append(aa)
                            wrong_times.append(temp1[aa])
                            continue
                    else:
                        if aa+1 in wrong_index:
                            continue
                        else:
                            wrong_index.append(aa+1)
                            wrong_times.append(temp1[aa+1])
                            continue
            #% x & y array to be used for interpolation
            temp3 = np.asarray([tg for ii, tg in enumerate(temp1) if ii not in wrong_index])
            temp4 = np.asarray([tg for ii, tg in enumerate(temp2) if ii not in wrong_index])
            #=================================
            #% in case the first and/or last measurements are not reasonable!!
            t1 = list(copy.deepcopy(timearr)); t2 = list(copy.deepcopy(temp3))
            z1 = list(copy.deepcopy(s2)); z2 = list(copy.deepcopy(temp4))
            if temp3[0] != timearr[0]:
                self.logger.log("tsys starts with wrong values.. to fill in: using Mean of 15% of the data starting from left to right!")
                fil1 = t1.index(t2[0])
                avg_front = np.mean(temp4[:round(len(t1)*0.15+0.5)])
                for j in range(fil1):
                    z2.insert(0, avg_front)
                    t2.insert(j, t1[j])

            else:
                self.logger.log("Front: good to go (the first data point looks ok)")

            if temp3[-1] != timearr[-1]:
                self.logger.log("tsys ends with wrong values.. to fill in: using Mean of 15% of the data starting from right to left!")
                for k in np.arange(-1, np.negative(len(t1)), -1):
                    if t2[-1] == t1[k]:
                        fil2 = k
                        break

                avg_end = np.mean(temp4[-round(len(t1)*0.15+0.5):])
                for m in np.arange(fil2+1, 0, +1):
                    z2.append(avg_end)
                    t2.append(t1[m])

            else:
                self.logger.log("End: good to go (the last data point looks ok)")
            #=================================
            #% Linear interpolation
            gtg_time = np.asarray(t2)
            gtg_tsys = np.asarray(z2)
            linv = interp1d(gtg_time, gtg_tsys, kind='linear')

            plt.plot(timearr, linv(timearr), markersize=1.5, marker='o', linewidth=0.8, alpha=1.0, zorder=0, color=ccc, label=k1)
            if wrong_index:
                plt.plot(wrong_times, linv(wrong_times), linewidth=0, marker='*', markersize=9, color='#2ca02c', label='Lin. Interp.')
                self.logger.log("\nInterpolation has been performed!!")
            else:
                self.logger.og("\nAll measurements are reasonable!!")

            self.logger.log(" ",k1)
            self.logger.log("  --> Number of OUTLYING tsys values:", len(wrong_index))
            self.logger.log("---------------------------------------------")
            resu = np.asarray([round(g, 1) for g in linv(timearr)])   #% for recording later..
            dali.append(resu); del(resu)
            #% basic statistical properties
            vax = int(np.round(max(linv(timearr))))
            vin = int(np.round(min(linv(timearr))))
            vav = int(np.round(np.median(linv(timearr))))
            vst = int(np.round(np.std(linv(timearr))))
            #% plot layout
            fts3 = 9.5
            plt.legend(loc=2, fontsize=fts3, ncol=1, columnspacing=0.9, markerscale=1.2)
            fts1 = 10
            opz1 = np.linspace(plt.ylim()[0], plt.ylim()[1], 7)
            thegap = np.abs(np.abs(opz1[-1])-np.abs(opz1[-2]))
            plt.ylim(opz1[0]-thegap*0.0, opz1[-1]+thegap*1.3)
            plt.title("*"+self.antcode+"*: Max={}, Min={}, Med={}, Std={}".format(vax,vin,vav,vst), color='r', fontsize=fts1)
            plt.xlabel("Decimal day")
            plt.ylabel("Tsys [K]")
            plt.minorticks_on()
            plt.tick_params(axis='both', which='major', length=5, direction='in', pad=2, color='k')
            plt.tick_params(axis='both', which='minor', length=3, direction='in', pad=2, color='k')
        # plt.tight_layout()
        plt.show()
        plt.subplots_adjust(left=0.045, bottom=0.050, right=0.975, top=0.965, wspace=0.190, hspace=0.345)
        plt.savefig(self.logger.pref+"_"+self.session+"_" + self.freq+"Tsys2"+"_"+str(siglev)+"-sig"+".png", format='png', dpi=150, transparent=False)
        #% make a new list with the final results
        td,tmh = np.genfromtxt(self.tsys, usecols=(0,1), unpack=True, invalid_raise=False, missing_values='', usemask=False, dtype=str)
        finl = []
        for n1 in range(len(td)):
            fin99 = str(td[n1]) + " " + str(tmh[n1])
            for n2 in dali:
                fin99 = fin99 + " " + str(n2[n1])
                if len(fin99.split()) == self.totcol + 2:
                    fin99 = fin99 + "  \n"
                continue
            finl.append(fin99)
            del fin99
            if n1 == len(td)-1:
                finl.append("/\n")
            else:
                continue
        finl.insert(0, self.index)
        finl.insert(0, self.gain)
        self.smip = copy.deepcopy(finl)
        outn = self.logger.pref +"_"+ self.session + "_" + self.freq+"Tsys2"+"_"+str(siglev)+"-sig"+".dat"
        with open(outn, 'w') as ft24:
            ft24.writelines(self.smip)
        self.logger.log("\n*If the result looks OK with the smoothed out data,")
        self.logger.log("*--> Use the output 'xxx_Tsys2.dat'. That's it!")
        self.logger.log("\nOUTPUT: Tsys2.dat file, self.smip, and a png file (i.e., xxxx_Tsys2.png)")


class TypeC:
    """
    SINGLE TSYS COLUMNS: suitable for special antennas PV, NN, GL, ...
    output: self.dall, self.tsys, self.arg, self.antcode, self.gain, self.index, self.totcol, self.freq, and self.legend
    """
    def __init__(self, antcode, antfile, session_codes, interactive_mode, wavelength):
        self.wavelength= wavelength
        self.interactive_mode = interactive_mode
        self.antcode = antcode
        self.antfile = antfile
        self.session = session_codes
        self.logger = Logger()
        self.logger.genlog()
        # self.session = None
        self.dall = None
        self.tsysout = None
        self.tarr = None
        self.rcp = None
        self.lcp = None
        self.aa = None    
        self.intp = None
        self.gain = None
        self.arg = None
        self.index = None
        self.totcol = None
        self.columnorder = None
        self.freq = None
        self.legend = None
        self.smip = None
        self.tsysread()

    def parser(self):
        if self.interactive_mode:
            arg = []
            reciv = input("\n<Q> Give antab file(s) of this station/session \n (e.g., 'c211aef.antabfs c211bef.antabfs c211cef.antabfs' --> space between files) \n: ")
            arg = reciv.split()
            self.session = input("\n<Q> Give all the session codes of the input data? (e.g., 'abc' for three a,b,c antab files) \n: ")
            self.logger.log("\n Total %.0f sessions \n" % len(arg))
        else:
            arg = self.antfile
        self.arg = arg
        N = len(arg)
        dset = []
        for i in range(N):
            ts1 = open(arg[i])
            ts2 = ts1.readlines()
            dset.append(ts2)

        dat = dset
        if len(dat) == 1:
            datfull = dat[0]
        else:
            datfull = dat[0]
            for i in range(len(dat)):
                if i == 0:
                    continue
                else:
                    datfull = datfull + dat[i]

        for i, j in enumerate(datfull[:30]):
            print(str(i) + ")", j)
        ################################################################        
        gain_line = ''
        poly_line = ''
        index_line = ''
        combined_data = '\n'.join(datfull)

        gain_match = re.search(r'GAIN[^/]+POLY.*?/', combined_data, re.DOTALL)
        if gain_match:
            gain_line = gain_match.group(0).split('POLY')[0].strip()
            if self.antcode == 'PV':
                gain_line = gain_line.replace('GAIN', f'GAIN {self.antcode}')

        poly_match = re.search(r'POLY=.*?/', combined_data, re.DOTALL)
        if poly_match:
            poly_line = poly_match.group(0).strip().replace('\n', '')

        index_match = re.search(r'[Tt][Ss][Yy][Ss].*[Ii][Nn][Dd][Ee][X].*?/', combined_data, re.IGNORECASE)
        if index_match:
            index_line = index_match.group(0).strip().replace('\n', '')
        else:
            print("INDEX line not available, please check this information")

        ################################################################

        self.logger.log("\n<Q> Find the *GAIN LINE* \n(e.g., GAIN AA ELEV DPFU=1.0) \n (e.g., GAIN EF ELEV DPFU=0.14,0.14 FREQ=84000,95500) \n **FREQ is optional and you can skip it \n**Do NOT end this line with '/'\n")
        if self.interactive_mode:
            while True:
                if self.antcode in ['KT', 'KU', 'KY']:
                    print("..... oh... KVN...\nKVN data might contain more GAIN lines than needed. Please select GAIN line with the right frequency!\n")
                provided_gain = input("\nIs this the right GAIN line? Do not forget the station code!\n(Press ENTER to confirm or input a different line): \n" + gain_line + "\n: ")
                if provided_gain.strip() == "":
                    gaininfo = gain_line
                    break
                else:
                    gaininfo = provided_gain
                    if not gaininfo.startswith('GAIN'):
                        print("Please provide a line starting with 'GAIN'.")
                    else:
                        break
        else:
            self.logger.log("\nGAIN line): \n" + gain_line + "\n: ")
            gaininfo = gain_line
        self.logger.log("\n<Q> Find the *POLY LINE* \n(e.g., POLY=1.0 /) \n(e.g., POLY=0.727089,0.947364E−02,-0.822152E−04 /) \n**End this line with '/' \n ")
        if self.interactive_mode:
            while True:
                provided_poly = input("\nIs this the right POLY line? \n(Press ENTER to confirm or input a different line): \n" + poly_line + "\n: ")
                if provided_poly.strip() == "":
                    gainpoly = poly_line
                    break
                else:
                    gainpoly = provided_poly
                    if not gainpoly.startswith('POLY'):
                        print("Please provide a line starting with 'POLY'.")
                    else:
                        break
        else:
            self.logger.log("\nPOLY line: \n" + poly_line + "\n: ")
            gainpoly = poly_line
        self.logger.log("\n<Q> Check the number of Tsys columns and the R/L order! \n (e.g., TSYS EF INDEX='R1','R2','R3','R4','R5','R6','R7','R8','L1','L2','L3','L4','L5','L6','L7','L8' /) \n **End the line with '/'\n ")
        if self.interactive_mode:
            while True:
                provided_index = input("\nIs this the right TSYS INDEX line? \n(Press ENTER to confirm or input a different line): \n" + index_line + "\n: ")
                if provided_index.strip() == "":
                    indexinfo = index_line
                    break
                else:
                    indexinfo = provided_index
                    if not indexinfo.startswith('TSYS'):
                        print("Please provide a line starting with 'TSYS'.")
                    else:
                        break
        else:
            self.logger.log("\nTSYS INDEX line: \n" + index_line + "\n: ")
            indexinfo = index_line
        gaininfo = gaininfo + " " + gainpoly
        self.gain = gaininfo + '\n'
        self.index = indexinfo + '\n'
        self.datfull = datfull
        ################################################################
        # end parser for GAIN, POLY, INDEX
        ################################################################

        #now: parser for tsys columns
        if self.antcode == 'PV':
            self.logger.log("\n<Q> How many tsys-columns in the data? (give a number e.g., 16 or 8 or 2) \n ")
            # hd_legnum = input("Since it's " + self.antcode + " there are probably 2 tsys-columns, (Press ENTER to confirm or give the right number)\n: ")
            if self.interactive_mode:
                hd_legnum = input("Since it's " + self.antcode + " there are probably 2 tsys-columns, (Press ENTER to confirm or give the right number)\n: ")
                if hd_legnum.strip() == "":
                    hd_legnum = 2
                    hd_legnum = int(hd_legnum)
                    self.totcol = hd_legnum
                    self.logger.log("\nTotal", hd_legnum, "columns of tsys measurements.")
                else:
                    hd_legnum = int(hd_legnum)
                    self.totcol = hd_legnum
                    self.logger.log("\nTotal", hd_legnum, "columns of tsys measurements.")
            else:
                self.logger.log("Since it's " + self.antcode + " there are probably 2 tsys-columns\n: ")
                hd_legnum = 2
                hd_legnum = int(hd_legnum)
                self.totcol = hd_legnum
                self.logger.log("\nTotal", hd_legnum, "columns of tsys measurements.")
        elif self.antcode in ['KU', 'KY', 'KT']:
            self.logger.log("\n<Q> How many tsys-columns in the data? (give a number e.g., 16 or 8 or 2) \n ")
            if self.interactive_mode:
                hd_legnum = input("Since it's " + self.antcode + " there are probably 16 tsys-columns, (Press ENTER to confirm or give the right number)\n: ")
                if hd_legnum.strip() == "":
                    hd_legnum = 16
                    hd_legnum = int(hd_legnum)
                    self.totcol = hd_legnum
                    self.logger.log("\nTotal", hd_legnum, "columns of tsys measurements.")
                else:
                    hd_legnum = int(hd_legnum)
                    self.totcol = hd_legnum
                    self.logger.log("\nTotal", hd_legnum, "columns of tsys measurements.")
            else:
                self.logger.log("Since it's " + self.antcode + " there are probably 16 tsys-columns\n: ")
                hd_legnum = 16
                hd_legnum = int(hd_legnum)
                self.totcol = hd_legnum
                self.logger.log("\nTotal", hd_legnum, "columns of tsys measurements.")
        else:
            hd_legnum = input("\n<Q> How many tsys-columns in the data? (give a number e.g., 16 or 8 or 2) \n: ")
            hd_legnum = int(hd_legnum)
            self.totcol = hd_legnum
            self.logger.log("\nTotal", hd_legnum, "columns of tsys measurements.")

        #now: parser for legend
        lgst = []
        self.logger.log("\n<Q> For plot-legends later, see the following and select your data order by giving its number \n 1 -  R1,R2,R3,R4,L1,L2,L3,L4 \n 2 -  R1,R2,R3,R4,R5,R6,R7,R8,L1,L2,L3,L4,L5,L6,L7,L8  --> (usually EU ants.) \n 3 -  R1,R,L,R,L,R,L,R,L                         ---> (usually VLBA ants.) \n 4 -  R1,L1,R2,L2,R3,L3,R4,L4,R5,L5,R6,L6,R7,L7,R8,L8 \n 5 -  RCP,LCP                                      ------> (usually those Special ants.) \n (..not there? then Give it manually here: e.g., R1,L1,L2, ...)\n")
        if self.interactive_mode:
            while True:
                if self.antcode in ['KU', 'KY', 'KT']:
                    provided_indleg = input("Since it's " + self.antcode + " the data order is probably \n2 -  R1,R2,R3,R4,R5,R6,R7,R8,L1,L2,L3,L4,L5,L6,L7,L8, \n(Press ENTER to confirm or give the right number): \n ")
                    if provided_indleg.strip() == "":
                        legtem = 'R1,R2,R3,R4,R5,R6,R7,R8,L1,L2,L3,L4,L5,L6,L7,L8'
                        lgst = legtem.split(',')
                        break
                if self.antcode == 'PV':
                    provided_indleg = input("Since it's " + self.antcode + " the data order is probably \n5 -  'RCP,LCP' \n(Press ENTER to confirm or give the right number): \n ")
                    if provided_indleg.strip() == "":
                        legtem = 'RCP,LCP'
                        lgst = legtem.split(',')
                        break
                else:
                    legtem = input(": ")
                    if len(legtem) == 1:
                        legtem = int(legtem)
                        if legtem == 1:
                            legtem = 'R1,R2,R3,R4,L1,L2,L3,L4'
                        elif legtem == 2:
                            legtem = 'R1,R2,R3,R4,R5,R6,R7,R8,L1,L2,L3,L4,L5,L6,L7,L8'
                        elif legtem == 3:
                            legtem = 'R1,R,L,R,L,R,L,R,L'
                        elif legtem == 4:
                            legtem = 'R1,L1,R2,L2,R3,L3,R4,L4,R5,L5,R6,L6,R7,L7,R8,L8'
                        elif legtem == 5:
                            legtem = 'RCP,LCP'
                            lgst = legtem.split(',')
                            break
        else:
            if self.antcode in ['KU', 'KY', 'KT']:
                legtem = 'R1,R2,R3,R4,R5,R6,R7,R8,L1,L2,L3,L4,L5,L6,L7,L8'
                lgst = legtem.split(',')
            if self.antcode == 'PV':
                legtem = 'RCP,LCP'
                lgst = legtem.split(',')
        self.legend = lgst


    def tsysread(self):
        print("===========================")
        print('Starting tsysread for '+ self.antcode)
        print("===========================")
        self.logger.log("\n@Task: tsysread")

        print(" ")
        print("'ufo()' has some useful information!")
        if self.interactive_mode:
            self.freq = input("<Q> Observing wavelength? (e.g., 3mm) \n:")
            
        else:
            self.logger.log("\nAn observing wavelength of "+ self.wavelength +" is assumed\n")
            self.freq= self.wavelength
        if self.antcode in ['PV']:
            for i in os.listdir():
                if i.endswith('.antab') or i.endswith('.antabfs'):
                    self.logger.log(i)

            #start parser
            self.parser()

            datfull = self.datfull # ADDED THIS
            self.dall = datfull
            temdat = copy.deepcopy(datfull)
            for j in temdat[:]:
                if j[0].isdigit():
                    continue
                else:
                    temdat.remove(j)

            self.tsys = temdat      # only actual data for plotting
            self.logger.log("\n***How it looks?")
            for zz in range(15):
                self.logger.log(temdat[zz])

            #--------------------------------------- another antennas with a sad ANTAB data format
            #--------------------------------------- another antennas with a sad ANTAB data format
            if self.antcode == 'PV':
                self.logger.log("\noh, this is Pico Veleta..")
                pvday = np.genfromtxt(self.tsys, usecols=0, dtype='str')
                pvtim = np.genfromtxt(self.tsys, usecols=1, dtype='str')
                pvrcp = np.genfromtxt(self.tsys, usecols=2, dtype='str')
                pvlcp = np.genfromtxt(self.tsys, usecols=3, dtype='str')
                finl = []
                for y in range(len(pvday)):
                    tline = pvday[y] + " " + pvtim[y] + " " + pvrcp[y] + " " + pvlcp[y] + " \n"
                    finl.append(tline)
                    del tline
                    continue

                self.tsys = None
                self.tsys = finl

            self.logger.log("\nOUTPUT: self.dall, self.tsys, self.arg, self.antcode, self.gain, self.index, self.totcol, self.freq, self.session, and self.legend")
            self.logger.log(" **just-in-case**..if Something went wrong (incorrect inputs among the above Pars.)")
            self.logger.log("   ,one can give it manually (e.g., self.xx = !!) or repeat this function again.")

        elif self.antcode == 'NN':
            if self.interactive_mode:

                self.logger.log("\n........................Oh, this is Noema..")
                for i in listdir():
                    if i.endswith('.antab') or i.endswith('.antabfs') or i.endswith('.asc'):
                        self.logger.log(i)

                #since for NN there is no GAIN or POLY information given in the data, one needs to type it in manually
                arg = []
                reciv = input("\n<Q> Give antab file(s) in the order of a_RCP, a_LCP --> b_RCP, b_LCP \n (e.g., c211a-Nn-rcp.asc c211a-Nn-lcp.asc c211b-Nn-rcp.asc c211b-Nn-lcp.asc ... --> space between files) \n: ")
                arg = reciv.split()
                self.session = input("\n<Q> Give all the session codes of the input data? (e.g., 'abc' for three a,b,c antab files) \n: ")
                self.logger.log("\n Total %.0f sessions (for Noema, recorded RCP/LCP separately) \n" % int(len(arg)/2))
                self.arg = arg

                temn1 = open(arg[0])
                temn2 = temn1.readlines()
                for i,j in enumerate(temn2[:21]):
                    print(str(i)+")", j)
                gaininfo = input("\n<Q> Find the *GAIN LINE* from above output lines and re-write it down here \n!! A line for POLY is coming, so dont include it here (see example below) !! \n (e.g., GAIN AA ELEV DPFU=1.0) \n (e.g., GAIN EF ELEV DPFU=0.14,0.14 FREQ=84000,95500) \n\n**For some special antennas (i.e., PV,NN,GL), find their Gain Info. separately and give it here \n**FREQ is optional and you can skip it \n**Do NOT end this line with '/' \n: ")
                gainpoly = input("\n<Q> Find the *POLY LINE* from above output lines and re-write it down here \n(e.g., POLY=1.0 /) \n(e.g., POLY=0.727089,0.947364E−02,-0.822152E−04 /) \n\n**For some special antennas (i.e., PV,NN,GL), find their Gain Info. separately and give it here \n**End this line with '/' \n: ")
                gaininfo = gaininfo + " " + gainpoly
                indexinfo = input("\n<Q> In the output lines above, check the number of Tsys columns and the R/L order! \n (e.g., TSYS EF INDEX='R1','R2','R3','R4','R5','R6','R7','R8','L1','L2','L3','L4','L5','L6','L7','L8' /) \n (e.g., TSYS PV INDEX='R1:8','L1:8' /) \n (e.g., TSYS NN INDEX='R1:8','L1:8' /) \n **End the line with '/' \n: ")
                self.gain = gaininfo + '\n'
                self.index = indexinfo + '\n'
                hd_legnum = input("\n<Q> How many tsys-columns in the data? (give a number e.g., 16 or 8 or 2) \n: ")
                hd_legnum = int(hd_legnum)
                self.totcol = hd_legnum
                self.logger.log("\nTotal", hd_legnum, "columns of tsys measurements.")

                lgst = []
                self.logger.log("\n<Q> For plot-legends later, see the following and select your data order by giving its number \n 1 -  R1,R2,R3,R4,L1,L2,L3,L4 \n 2 -  R1,R2,R3,R4,R5,R6,R7,R8,L1,L2,L3,L4,L5,L6,L7,L8  --> (usually EU ants.) \n 3 -  R1,L1,R2,L2,R3,L3,R4,L4                         ---> (usually VLBA ants.) \n 4 -  RCP,LCP                                      ------> (usually those Special ants.) \n (..not there? then Give it manually here: e.g., R1,L1,L2, ...)")
                legtem = input(": ")
                if len(legtem) == 1:
                    legtem = int(legtem)
                    if legtem == 1:
                        legtem = 'R1,R2,R3,R4,L1,L2,L3,L4'
                    elif legtem == 2:
                        legtem = 'R1,R2,R3,R4,R5,R6,R7,R8,L1,L2,L3,L4,L5,L6,L7,L8'
                    elif legtem == 3:
                        legtem = 'R1,L1,R2,L2,R3,L3,R4,L4'
                    elif legtem == 4:
                        legtem = 'RCP,LCP'

                lgst = legtem.split(',')
                self.legend = lgst

                collecting = []
                for o in range(len(self.arg)):      # to select specific data column(s) from ANTAB
                    if ((o % 2) == 0) and ('rcp' in self.arg[o]):       #% if o is even (RCP)
                        rcv1 = open(self.arg[o]); rcv2 = rcv1.readlines()
                        lcv1 = open(self.arg[o+1]); lcv2 = lcv1.readlines()
                        rcv3 = []; lcv3 = []
                        for i in rcv2:
                            if (not i.startswith('!')) and (not i.startswith('/')):
                                rcv3.append(i)
                        for j in lcv2:
                            if (not j.startswith('!')) and (not j.startswith('/')):
                                lcv3.append(j)
                        for k in range(len(rcv3)):
                            itm1 = np.genfromtxt(rcv3, usecols=0, dtype='str')[k]
                            itm2 = np.genfromtxt(rcv3, usecols=1, dtype='str')[k]
                            itm3 = np.genfromtxt(rcv3, usecols=2, dtype='str')[k]
                            itm4 = np.genfromtxt(lcv3, usecols=2, dtype='str')[k]
                            itm99 = itm1 + " " + itm2 + " " + itm3 + " " + itm4 + " \n"
                            collecting.append(itm99)
                            del(itm99)
                    self.tsys = collecting

                N = len(self.arg)
                dset = []
                for i in range(N):
                    ts1 = open(self.arg[i])
                    ts2 = ts1.readlines()
                    dset.append(ts2)

                dat = dset
                if len(dat) == 1:
                    datfull = dat[0]
                else:
                    datfull = dat[0]
                    for i in range(len(dat)):
                        if i == 0:
                            continue
                        else:
                            datfull = datfull + dat[i]    # combine all sessions into one data group
                self.dall = datfull    # here "datfull" is just a collection of raw antab file(s).
                self.logger.log("\nOUTPUT: self.dall, self.tsys, self.arg, self.antcode, self.gain, self.index, self.totcol, self.freq, self.session, and self.legend")
                self.logger.log(" **just-in-case**..if Something went wrong (incorrect inputs among the above Pars.)")
                self.logger.log("   ,one can give it manually (e.g., self.xx = !!) or repeat this function again.")
            else:

                self.logger.log("\n........................Oh, this is Noema..")
                arg = []
                self.logger.log("\n<Q> Give antab file(s) in the order of a_RCP, a_LCP --> b_RCP, b_LCP \n (e.g., c211a-Nn-rcp.asc c211a-Nn-lcp.asc c211b-Nn-rcp.asc c211b-Nn-lcp.asc ... --> space between files) \n: ")
                reciv = ' '.join(self.antfile)
                print(reciv)
                arg = reciv.split()
                self.logger.log("\n<Q> Give all the session codes of the input data? (e.g., 'abc' for three a,b,c antab files) \n: ")
                
                self.session = self.session
                print(self.session)
                self.logger.log("\n Total %.0f sessions (for Noema, recorded RCP/LCP separately) \n" % int(len(arg)/2))
                self.arg = arg

                temn1 = open(arg[0])
                temn2 = temn1.readlines()
                for i,j in enumerate(temn2[:21]):
                    print(str(i)+")", j)
                self.logger.log('-------------------------------------------------------\n-----------------------------------------------------------------------\n-------------------------------------------------------------------')
                self.logger.log("The following GAIN, POLY and INDEX lines are assumed:\n")
                gaininfo ='GAIN NN ELEV DPFU = 0.149,0.151'
                gainpoly ='POLY= 0.91437, 0.005815, -0.00014101, +1.5319e-6, -7.5242e-9 /'
                self.logger.log(gaininfo, "\n", gainpoly)
                gaininfo = gaininfo + " " + gainpoly
                indexinfo = "Tsys NN   timeoff = 0.0  FT = 1.0 INDEX='R1:8','L1:8' /"
                self.logger.log(indexinfo)
                self.gain = gaininfo + '\n'
                self.index = indexinfo + '\n'

                hd_legnum = 2
                hd_legnum = int(hd_legnum)
                self.totcol = hd_legnum
                self.logger.log("\nTotal", hd_legnum, "columns of tsys measurements.")

                lgst = []
                legtem = 'RCP,LCP'
                lgst = legtem.split(',')
                lgst = legtem.split(',')
                self.legend = lgst

                collecting = []
                for o in range(len(self.arg)):      # to select specific data column(s) from ANTAB
                    if ((o % 2) == 0) and ('rcp' in self.arg[o]):       #% if o is even (RCP)
                        rcv1 = open(self.arg[o]); rcv2 = rcv1.readlines()
                        lcv1 = open(self.arg[o+1]); lcv2 = lcv1.readlines()
                        rcv3 = []; lcv3 = []
                        for i in rcv2:
                            if (not i.startswith('!')) and (not i.startswith('/')):
                                rcv3.append(i)
                        for j in lcv2:
                            if (not j.startswith('!')) and (not j.startswith('/')):
                                lcv3.append(j)
                        for k in range(len(rcv3)):
                            itm1 = np.genfromtxt(rcv3, usecols=0, dtype='str')[k]
                            itm2 = np.genfromtxt(rcv3, usecols=1, dtype='str')[k]
                            itm3 = np.genfromtxt(rcv3, usecols=2, dtype='str')[k]
                            itm4 = np.genfromtxt(lcv3, usecols=2, dtype='str')[k]
                            itm99 = itm1 + " " + itm2 + " " + itm3 + " " + itm4 + " \n"
                            collecting.append(itm99)
                            del(itm99)
                    self.tsys = collecting

                N = len(self.arg)
                dset = []
                for i in range(N):
                    ts1 = open(self.arg[i])
                    ts2 = ts1.readlines()
                    dset.append(ts2)

                dat = dset
                if len(dat) == 1:
                    datfull = dat[0]
                else:
                    datfull = dat[0]
                    for i in range(len(dat)):
                        if i == 0:
                            continue
                        else:
                            datfull = datfull + dat[i]    # combine all sessions into one data group
                self.dall = datfull    # here "datfull" is just a collection of raw antab file(s).
                self.logger.log("\nOUTPUT: self.dall, self.tsys, self.arg, self.antcode, self.gain, self.index, self.totcol, self.freq, self.session, and self.legend")
                self.logger.log(" **just-in-case**..if Something went wrong (incorrect inputs among the above Pars.)")
                self.logger.log("   ,one can give it manually (e.g., self.xx = !!) or repeat this function again.")
            #--------------------------------------- another antennas with a sad ANTAB data format
            #--------------------------------------- another antennas with a sad ANTAB data format
        elif self.antcode == 'GL':
            self.logger.log("\noh, this is Greenland (GLT)..")
            for i in listdir():
                if i.endswith('.antab') or i.endswith('.antabfs') or i.endswith('.txt'):
                    self.logger.log(i)
            #since for GLT there is no GAIN or POLY information given in the data, one needs to type it in manually

            arg = []
            reciv = input("\n<Q> Give antab file(s) in the order of alphabetical order (or RCP -> LCP) \n (e.g., gltsys_c211a.txt gltsys_c211b.txt gltsys_c211c.txt ... --> space between files) \n: ")
            arg = reciv.split()
            self.session = input("\n<Q> Give all the session codes of the input data? (e.g., 'abc' for three a,b,c antab files) \n: ")
            self.logger.log("\n Total %.0f sessions \n" % len(arg))
            self.arg = arg

            temn1 = open(arg[0])
            temn2 = temn1.readlines()
            for i,j in enumerate(temn2[:21]):
                print(str(i)+")", j)

            gaininfo = input("\n<Q> Find the *GAIN LINE* from above output lines and re-write it down here \n!! A line for POLY is coming, so dont include it here (see example below) !! \n (e.g., GAIN AA ELEV DPFU=1.0) \n (e.g., GAIN EF ELEV DPFU=0.14,0.14 FREQ=84000,95500) \n\n**For some special antennas (i.e., PV,NN,GL), find their Gain Info. separately and give it here \n**FREQ is optional and you can skip it \n**Do NOT end this line with '/' \n: ")
            gainpoly = input("\n<Q> Find the *POLY LINE* from above output lines and re-write it down here \n(e.g., POLY=1.0 /) \n(e.g., POLY=0.727089,0.947364E−02,-0.822152E−04 /) \n\n**For some special antennas (i.e., PV,NN,GL), find their Gain Info. separately and give it here \n**End this line with '/' \n: ")
            gaininfo = gaininfo + " " + gainpoly
            indexinfo = input("\n<Q> In the output lines above, check the number of Tsys columns and the R/L order! \n(e.g., TSYS EF INDEX='R1','R2','R3','R4','R5','R6','R7','R8','L1','L2','L3','L4','L5','L6','L7','L8' /) \n**End the line with '/' \n**A recent GL Index is, \n  --> TSYS GL INDEX='R1:8','L1:8' /       (but, double-check its antab file and data in it) \n: ")
            self.gain = gaininfo + '\n'
            self.index = indexinfo + '\n'
            hd_legnum = input("\n<Q> How many tsys-columns in the data? (give a number e.g., 16 or 8 or 2) \n: ")
            hd_legnum = int(hd_legnum)
            self.totcol = hd_legnum
            self.logger.log("\nTotal", hd_legnum, "columns of tsys measurements.")


            lgst = []
            self.logger.log("\n<Q> For plot-legends later, see the following and select your data order by giving its number \n 1 -  R1,R2,R3,R4,L1,L2,L3,L4 \n 2 -  R1,R2,R3,R4,R5,R6,R7,R8,L1,L2,L3,L4,L5,L6,L7,L8  --> (usually EU ants.) \n 3 -  R1,L1,R2,L2,R3,L3,R4,L4                         ---> (usually VLBA ants.) \n 4 -  R1,L1,R2,L2,R3,L3,R4,L4,R5,L5,R6,L6,R7,L7,R8,L8 \n 5 -  RCP,LCP                                      ------> (usually those Special ants.) \n (..not there? then Give it manually here: e.g., R1,L1,L2, ...)")
            legtem = input(": ")
            if len(legtem) == 1:
                legtem = int(legtem)
                if legtem == 1:
                    legtem = 'R1,R2,R3,R4,L1,L2,L3,L4'
                elif legtem == 2:
                    legtem = 'R1,R2,R3,R4,R5,R6,R7,R8,L1,L2,L3,L4,L5,L6,L7,L8'
                elif legtem == 3:
                    legtem = 'R1,L1,R2,L2,R3,L3,R4,L4'
                elif legtem == 4:
                    legtem = 'R1,L1,R2,L2,R3,L3,R4,L4,R5,L5,R6,L6,R7,L7,R8,L8'
                elif legtem == 5:
                    legtem = 'RCP,LCP'

            lgst = legtem.split(',')
            self.legend = lgst
            print(self.legend)
            N = len(self.arg)
            dset = []
            for i in range(N):
                ts1 = open(self.arg[i])
                ts2 = ts1.readlines()
                dset.append(ts2)

            dat = dset
            if len(dat) == 1:
                datfull = dat[0]
            else:
                datfull = dat[0]
                for i in range(len(dat)):
                    if i == 0:
                        continue
                    else:
                        datfull = datfull + dat[i]    # combine all sessions into one data group
            self.dall = datfull    # here "datfull" is just a collection of raw antab file(s).

            filtered = copy.deepcopy(datfull)
            for i in filtered[:]:
                if i[0].isdigit():
                    continue
                else:
                    filtered.remove(i)
            for hh,pp in enumerate(filtered,0):
                if len(pp.split()) < 3:
                    i99 = pp.split()[0] + " " + pp.split()[1] + " " + str(0) + " " + str(0) + " \n"
                    filtered[hh] = copy.deepcopy(i99)
                    del i99
                    continue
                elif len(pp.split()) == 3:
                    i99 = pp.split()[0] + " " + pp.split()[1] + " " + pp.split()[2].split("/")[0] + " " + pp.split()[2].split("/")[1] + " \n"
                    filtered[hh] = copy.deepcopy(i99)
                    del i99
                    continue
            self.tsys = filtered
            self.logger.log("\nOUTPUT: self.dall, self.tsys, self.arg, self.antcode, self.gain, self.index, self.totcol, self.freq, self.session, and self.legend")
            self.logger.log(" **just-in-case**..if Something went wrong (incorrect inputs among the above Pars.)")
            self.logger.log("   ,one can give it manually (e.g., self.xx = !!) or repeat this function again.")

    def dbcheck(self):
        """
        To check if the data contains unreliable values.
        run it like e.g., self.dbcheck()
        output: .dat file, self.tarr, self.rcp, self.lcp, self.tsysout
        """
        self.logger.log("===========================")
        self.logger.log("\n@Task: dbcheck")
        self.logger.log("===========================")
        new_t = []
        simpletest = self.tsys[0].split()[1]
        if simpletest.count(":") == 2:
            for l in self.tsys:
                tem_d = l.split()[0]        #% day array
                tem_hms = l.split()[1]       #% hh:mm:ss array
                #% converting all of them to a day with decimal numbers
                div_h = tem_hms.split(':')[0]
                div_m = tem_hms.split(':')[1]
                div_s = tem_hms.split(':')[2]
                dayform = (int(div_h)/24.) + (float(div_m)/60. /24.) + (float(div_s)/60. /60. /24.) + int(tem_d)
                new_t.append(dayform)

        elif simpletest.count(":") == 1:
            for l in self.tsys:
                tem_d = l.split()[0]        #% day array
                tem_hm = l.split()[1]       #% hh:mm array
                div_h = tem_hm.split(':')[0]
                div_m = tem_hm.split(':')[1]
                dayform = (int(div_h)/24.) + (float(div_m)/60. /24.) + int(tem_d)
                new_t.append(dayform)
                #% --> actual x-array for interpolation

        new_t = np.asarray(new_t)
        self.tarr = new_t

        if self.totcol == 2:      # --> spectral windows in a form of 'R1:8','L1:8'
            R1_8 = np.genfromtxt(self.tsys, usecols=2, unpack=True, invalid_raise=False, missing_values='', usemask=False, dtype=float)      # always RCP first!
            L1_8 = np.genfromtxt(self.tsys, usecols=3, unpack=True, invalid_raise=False, missing_values='', usemask=False, dtype=float)      # always LCP second!
            self.rcp = R1_8
            self.lcp = L1_8
            self.logger.log("\nFor NN/GL/KVN/PV, the INDEX of 'R1:8' & 'L1:8' is assumed.")
            self.logger.log("\n***Now, checking any values of [999, >9999, 0, negative, btw 0 and 10]***")
            aa = 0
            for i in range(2):
                if i == 0:
                    s2 = R1_8
                    self.logger.log("..First, RCP")
                elif i == 1:
                    s2 = L1_8
                    self.logger.log("..Second, LCP")

                for k1,k2 in enumerate(s2,0):
                    if (k2 == 999) or (k2 == 999.0):
                        self.logger.log("-- 999 found!!")
                        aa = 1
                        continue
                    elif k2 >= 9999:
                        self.logger.log("-- higher than 9999 found!!")
                        aa = 1
                        continue
                    elif k2 == 0:
                        self.logger.log("-- 0 found!!")
                        aa = 1
                        continue
                    elif k2 < 0:
                        self.logger.log("-- NEGATIVE value found!!")
                        aa = 1
                        continue
                    elif 0 < k2 < 10:
                        self.logger.log("-- 0K < tsys < 10K, unreasonably low value found!!")
                        aa = 1
                        continue
            #-------------------------------------------------------------------------------------
        else:          # --> spectral windows in a form of e.g., 'R1','R2'..'L1','L2'..
            if self.interactive_mode:
                sep1 = input("\n<Q> Give the order of RCP/LCP in the Tsys columns from left to right, either 1 or 2? \n  1 (EU ant. form: = R1 R2 R3 R4 ... L1 L2 L3 L4 ...) \n  2 (VLBA ant. form: = R1 L1 R2 L2 R3 L3 R4 L4 ...) \n: ")
                sep1 = int(sep1)
                if sep1 == 1:
                    self.columnorder = "EU normally --> R1 R2 R3 R4 ... L1 L2 L3 L4 ..."
                elif sep1 == 2:
                    self.columnorder = "VLBA normally --> R1 L1 R2 L2 R3 L3 R4 L4 ..."
            else:
                if self.antcode in ['KT', 'KU', 'KY']:
                    self.logger.log("\nOrder of RCP/LCP in the Tsys columns from left to right is (EU ant. form: = R1 R2 R3 R4 ... L1 L2 L3 L4 ...) ")
                    sep1 = 1
                    self.columnorder = "EU normally --> R1 R2 R3 R4 ... L1 L2 L3 L4 ..."
                else:
                    sep1 = input("\n<Q> Give the order of RCP/LCP in the Tsys columns from left to right, either 1 or 2? \n  1 (EU ant. form: = R1 R2 R3 R4 ... L1 L2 L3 L4 ...) \n  2 (VLBA ant. form: = R1 L1 R2 L2 R3 L3 R4 L4 ...) \n: ")
                    sep1 = int(sep1)
                    if sep1 == 1:
                        self.columnorder = "EU normally --> R1 R2 R3 R4 ... L1 L2 L3 L4 ..."
                    elif sep1 == 2:
                        self.columnorder = "VLBA normally --> R1 L1 R2 L2 R3 L3 R4 L4 ..."
            tsyslist = []
            for kk1 in range(len(self.legend)):
                eachcol = np.genfromtxt(self.tsys, usecols=kk1+2, unpack=True, invalid_raise=False, missing_values='', usemask=False, dtype=float)
                tsyslist.append(eachcol)
                del eachcol

            self.rlcp = tsyslist
            self.logger.log("\n**Now, checking any values of [999, >9999, 0, negative, btw 0 and 10]")
            aa = 0
            if sep1 == 1:       #% first half -> rcp, and then lcp (e.g., EU antennas)
                for i in range(len(tsyslist)):
                    if i < int( len(tsyslist)/2 ):
                        s2 = tsyslist[i]
                        self.logger.log("..First, RCP !")
                        for k1,k2 in enumerate(s2,0):
                            if (k2 == 999) or (k2 == 999.0):
                                self.logger.log("-- 999 found!!")
                                aa = 1
                                continue
                            elif k2 >= 9999:
                                self.logger.log("-- higher than 9999 found!!")
                                aa = 1
                                continue
                            elif k2 == 0:
                                self.logger.log("-- 0 found!!")
                                aa = 1
                                continue
                            elif k2 < 0:
                                self.logger.log("-- negative value found!!")
                                aa = 1
                                continue
                            elif 0 < k2 < 10:
                                self.logger.log("-- 0K < tsys < 10K, unreasonably low value found!!")
                                aa = 1
                                continue
                    else:
                        s2 = tsyslist[i]
                        self.logger.log("..Second, LCP !")
                        for k1,k2 in enumerate(s2,0):
                            if (k2 == 999) or (k2 == 999.0):
                                self.logger.log("-- 999 found!!")
                                aa = 1
                                continue
                            elif k2 >= 9999:
                                self.logger.log("-- higher than 9999 found!!")
                                aa = 1
                                continue
                            elif k2 == 0:
                                self.logger.log("-- 0 found!!")
                                aa = 1
                                continue
                            elif k2 < 0:
                                self.logger.log("-- negative value found!!")
                                aa = 1
                                continue
                            elif 0 < k2 < 10:
                                self.logger.log("-- 0K < tsys < 10K, unreasonably low value found!!")
                                aa = 1
                                continue
            elif sep1 == 2:       #% rcp -> even, and lcp -> odd (e.g., VLBA)
                for i in range(len(tsyslist)):
                    if (i % 2) == 0:     #% even, thus RCP
                        s2 = tsyslist[i]
                        self.logger.log("..First, RCP !")
                        for k1,k2 in enumerate(s2,0):
                            if (k2 == 999) or (k2 == 999.0):
                                self.logger.log("-- 999 found!!")
                                aa = 1
                                continue
                            elif k2 >= 9999:
                                self.logger.log("-- higher than 9999 found!!")
                                aa = 1
                                continue
                            elif k2 == 0:
                                self.logger.log("-- 0 found!!")
                                aa = 1
                                continue
                            elif k2 < 0:
                                self.logger.log("-- negative value found!!")
                                aa = 1
                                continue
                            elif 0 < k2 < 10:
                                self.logger.log("-- 0K < tsys < 10K, unreasonably low value found!!")
                                aa = 1
                                continue
                    else:     #% odd, thus LCP
                        s2 = tsyslist[i]
                        self.logger.log("..Second, LCP !")
                        for k1,k2 in enumerate(s2,0):
                            if (k2 == 999) or (k2 == 999.0):
                                self.logger.log("-- 999 found!!")
                                aa = 1
                                continue
                            elif k2 >= 9999:
                                self.logger.log("-- higher than 9999 found!!")
                                aa = 1
                                continue
                            elif k2 == 0:
                                self.logger.log("-- 0 found!!")
                                aa = 1
                                continue
                            elif k2 < 0:
                                self.logger.log("-- negative value found!!")
                                aa = 1
                                continue
                            elif 0 < k2 < 10:
                                self.logger.log("-- 0K < tsys < 10K, unreasonably low value found!!")
                                aa = 1
                                continue
        self.showsys() # ADDED HERE
        if aa == 0:
            self.logger.log("\n---------------------------------------------")
            self.logger.log("***All measurements seem FINE! Ready to be used in the calibration!***")
            self.logger.log(" -Find an output file: xx_Tsys0.dat, in your working directory.")
            self.logger.log(" -But, STRONGLY RECOMMEND to plot for double-checking the data with showsys().")
            self.tsysout = copy.deepcopy(self.tsys)
            self.tsysout.append("/\n")
            # print('self.index')
            # print(self.index)
            self.tsysout.insert(0, self.index)
            self.tsysout.insert(0, self.gain)
            outn = self.logger.pref+"_" + self.antcode + "_"+ "_".join(self.session) + "_" + self.freq+ "_"+"Tsys0.dat"
            with open(outn, 'w') as ft36:
                ft36.writelines(self.tsysout)
            self.logger.log(" -If the plot looks ok, then, use the output file: xx_Tsys0.dat, here.")
        elif aa == 1:
            self.logger.log("\n---------------------------------------------")
            self.logger.log("***Unreliable value found!***")
            self.logger.log(" -Find an output file: xx_Tsys0.dat, in your working directory.")
            self.logger.log(" -Fix it further with showsys() and intpsys()")
            self.logger.log(" -Then, use self.tsysout with writetxt()")
            self.tsysout = copy.deepcopy(self.tsys)
            self.tsysout.append("/\n")
            self.tsysout.insert(0, self.index)
            self.tsysout.insert(0, self.gain)
            outn = self.logger.pref+"_" +self.antcode+"_"+ "_".join(self.session) + "_" + self.freq+ "_"+"Tsys0.dat"
            with open(outn, 'w') as ft36:
                ft36.writelines(self.tsysout)
            self.intpsys()
        self.logger.log("\nOUTPUT: Tsys0.dat file, self.tarr, self.rcp, self.lcp, self.tsysout")

    def showsys(self, Ymax=None, Ymin=None):
        """
        Plot tsys measurements
        Ymax: Y-axis range, Maximum value.
        Ymin: Y-axis range, Minimum value.
        """
        self.logger.log("\n@Task: showsys")
        self.logger.log("===========================")
        antcode = self.antcode
        rlab = "RCP at " + self.freq
        llab = "LCP at " + self.freq
        timearr = self.tarr
        rdat = self.rcp
        ldat = self.lcp
        self.logger.log("*Find a png file of the resultant plot.")
        fig_pol=plt.figure(figsize=(7,5.5))
        plt.rcParams['legend.frameon'] = 'False'
        plt.rcParams['xtick.labelsize'] = 13
        plt.rcParams['ytick.labelsize'] = 13
        rav = np.round(np.median(rdat), 1)
        rst = np.round(np.std(rdat), 1)
        lav = np.round(np.median(ldat), 1)
        lst = np.round(np.std(ldat), 1)

        plt.subplot(2, 1, 1)
        plt.plot(timearr, rdat, markersize=10, marker='+', linewidth=1.5, label=rlab, color='#1f77b4')
        if Ymax:
            if Ymin:
                plt.ylim([Ymin, Ymax])
            else:
                plt.ylim([0, Ymax])
        else:
            opz1 = np.linspace(plt.ylim()[0], plt.ylim()[1], 7)
            thegap = np.abs(np.abs(opz1[-1])-np.abs(opz1[-2]))
            plt.ylim(opz1[0]-thegap*0.0, opz1[-1]+thegap*1.3)
        plt.legend(loc=2, fontsize=14, numpoints=1, markerscale=1.5, labelspacing=0.1, handletextpad=0.5, borderaxespad=0.2, handlelength=1.0, handleheight=0.5)
        plt.xticks([])
        plt.ylabel("Tsys [K]", fontsize=16)
        plt.tick_params(axis='both', which='major', length=5, direction='in', pad=2, color='k')
        plt.tick_params(axis='both', which='minor', length=3, direction='in', pad=2, color='k')
        plt.minorticks_on()
        #% median
        plt.text(plt.xlim()[0] + (plt.xlim()[1]-plt.xlim()[0])*0.01, plt.ylim()[-1]*1.01, r'R-med=%.1fK' % rav, color='r', fontsize=11)
        plt.text(plt.xlim()[0] + (plt.xlim()[1]-plt.xlim()[0])*0.47, plt.ylim()[-1]*1.01, r'L-med=%.1fK' % lav, color='m', fontsize=11)
        #% std
        plt.text(plt.xlim()[0] + (plt.xlim()[1]-plt.xlim()[0])*0.24, plt.ylim()[-1]*1.01, r'R-std=%.1fK' % rst, color='r', fontsize=11)
        plt.text(plt.xlim()[0] + (plt.xlim()[1]-plt.xlim()[0])*0.70, plt.ylim()[-1]*1.01, r'L-std=%.1fK' % lst, color='m', fontsize=11)
        #% antenna code
        for i in range(2):
            plt.text(plt.xlim()[0] + (plt.xlim()[1]-plt.xlim()[0])*0.93, plt.ylim()[-1]*1.01, antcode, color='g', fontsize=18)

        plt.subplot(2, 1, 2)
        plt.plot(timearr, ldat, markersize=10, marker='+', linewidth=1.5, label=llab, color='#ff7f0e')
        if Ymax:
            if Ymin:
                plt.ylim([Ymin, Ymax])
            else:
                plt.ylim([0, Ymax])
        else:
            opz1 = np.linspace(plt.ylim()[0], plt.ylim()[1], 7)
            thegap = np.abs(np.abs(opz1[-1])-np.abs(opz1[-2]))
            plt.ylim(opz1[0]-thegap*0.0, opz1[-1]+thegap*1.3)
        plt.legend(loc=2, fontsize=14, numpoints=1, markerscale=1.5, labelspacing=0.1, handletextpad=0.5, borderaxespad=0.2, handlelength=1.0, handleheight=0.5)
        plt.ylabel("Tsys [K]", fontsize=16)
        plt.xlabel("Decimal day", fontsize=16)
        plt.tick_params(axis='both', which='major', length=5, direction='in', pad=2, color='k')
        plt.tick_params(axis='both', which='minor', length=3, direction='in', pad=2, color='k')
        plt.minorticks_on()
        plt.subplots_adjust(left=0.120, bottom=0.107, right=0.940, top=0.930, wspace=0.200, hspace=0.0)
        if Ymax:
            self.logger.log("\nNo output .png file.")
        else:
            # plt.savefig(self.logger.pref+"_"+self.session+"_"+self.freq+"Tsys0.png", format='png', dpi=150, transparent=False)
            plt.savefig(self.logger.pref+"_"+self.antcode+"_" + "_".join(self.session) + "_" + self.freq+ "_"+"Tsys0.png", format='png', dpi=150, transparent=False)
        # plt.show()
        self.logger.log("\n*If the plot looks OK and dbcheck did not report suspicious tsys values,")
        self.logger.log("*--> just use the output file of dbcheck. That's it!")
        self.logger.log("\nOUTPUT: A figure file (i.e., xxxx_Tsys0.png)")

    def intpsys(self, hnot=9999, lnot=10, sepif=False):
        """
        Linear interpolation to replace those unreliable tsys values.
        One can modify the cutoff of the high and low Tsys values.
        output: self.intp, self.tsysout, Figure
        """
        self.logger.log("\n@task: intpsys")
        self.logger.log("===========================")
        antcode = self.antcode
        rlab = "RCP at " + self.freq
        llab = "LCP at " + self.freq
        timearr = self.tarr
        rdat = self.rcp
        ldat = self.lcp
        self.logger.log("Find a png file of the resultant plot.")
        fig_pol=plt.figure(figsize=(7,5.5))
        plt.rcParams['legend.frameon'] = 'False'
        plt.rcParams['xtick.labelsize'] = 13
        plt.rcParams['ytick.labelsize'] = 13
        dali = []
        statrl = []
        for i in range(2):
            if i == 0:
                s2 = rdat
                self.logger.log("..\nFirst, RCP !")
                nxt = 1
                ccc = '#1f77b4'
                k1 = rlab
                p1 = plt.subplot(2, 1, nxt)
            elif i == 1:
                s2 = ldat
                self.logger.log("..\nSecond, LCP !")
                nxt = 2
                ccc = '#ff7f0e'
                k1 = llab
                p2 = plt.subplot(2, 1, nxt)

            temp1 = copy.deepcopy(timearr)
            temp2 = copy.deepcopy(s2)
            wrong_index = []
            wrong_times = []
            if sepif:
                self.logger.log("\n**Set high/low Tsys cutoff for each RCP/LCP IF**")
                self.logger.log("\n-->", self.legend[i])
                hgo = input("<Q> Tsys shouldn't EXCEED (press 'Enter' to skip): ")
                if len(hgo) == 0:
                    hgo = hnot
                else:
                    hgo = float(hgo)
                lgo = input("<Q> Tsys shouldn't be BELOW (press 'Enter' to skip): ")
                if len(lgo) == 0:
                    lgo = lnot
                else:
                    lgo = float(lgo)
                self.logger.log("\n!! any tsys values of 999, >=YOUR CHOICE, <0, 0, and 0<tsys<YOUR CHOICE will be removed !!")
                for aa,bb in enumerate(temp2, 0):
                    if (bb == 999) or (bb == 999.0):
                        wrong_index.append(aa)
                        wrong_times.append(temp1[aa])
                        continue
                    elif bb >= hgo:
                        wrong_index.append(aa)
                        wrong_times.append(temp1[aa])
                        continue
                    elif bb < 0:
                        wrong_index.append(aa)
                        wrong_times.append(temp1[aa])
                        continue
                    elif bb == 0:
                        wrong_index.append(aa)
                        wrong_times.append(temp1[aa])
                        continue
                    elif 0 < bb < lgo:
                        wrong_index.append(aa)
                        wrong_times.append(temp1[aa])
                        continue
            else:
                self.logger.log("\n!! any tsys values of 999, >=9999, <0, 0, and 0<tsys<10 will be removed !!")
                for aa,bb in enumerate(temp2, 0):
                    if (bb == 999) or (bb == 999.0):
                        wrong_index.append(aa)
                        wrong_times.append(temp1[aa])
                        continue
                    elif bb >= hnot:
                        wrong_index.append(aa)
                        wrong_times.append(temp1[aa])
                        continue
                    elif bb < 0:
                        wrong_index.append(aa)
                        wrong_times.append(temp1[aa])
                        continue
                    elif bb == 0:
                        wrong_index.append(aa)
                        wrong_times.append(temp1[aa])
                        continue
                    elif 0 < bb < lnot:
                        wrong_index.append(aa)
                        wrong_times.append(temp1[aa])
                        continue

            #% x & y array to be used for interpolation
            temp3 = np.asarray([tg for ii, tg in enumerate(temp1) if ii not in wrong_index])
            temp4 = np.asarray([tg for ii, tg in enumerate(temp2) if ii not in wrong_index])
            #=================================
            #% in case the first and/or last measurements are not reasonable!!
            t1 = list(copy.deepcopy(timearr)); t2 = list(copy.deepcopy(temp3))
            z1 = list(copy.deepcopy(s2)); z2 = list(copy.deepcopy(temp4))
            if temp3[0] != timearr[0]:
                self.logger.log("tsys starts with wrong/no value.. to fill in: using the mean of the first 15% of the data; i.e., from left to right")
                fil1 = t1.index(t2[0])
                avg_front = np.mean(temp4[:round(len(t1)*0.15+0.5)])
                for j in range(fil1):
                    z2.insert(0, avg_front)
                    t2.insert(j, t1[j])

            else:
                self.logger.log("Front: good to go (the first data point looks ok)")

            if temp3[-1] != timearr[-1]:
                self.logger.log("tsys ends with wrong/no value.. to fill in: using the mean of the last 15% of the data; i.e., from right to left")
                for k in np.arange(-1, np.negative(len(t1)), -1):
                    if t2[-1] == t1[k]:
                        fil2 = k
                        break

                avg_end = np.mean(temp4[-round(len(t1)*0.15+0.5):])
                for m in np.arange(fil2+1, 0, +1):
                    z2.append(avg_end)
                    t2.append(t1[m])

            else:
                self.logger.log("End: good to go (the last data point looks ok)")
            #=================================
            #% Linear interpolation
            gtg_time = np.asarray(t2)
            gtg_tsys = np.asarray(z2)
            linv = interp1d(gtg_time, gtg_tsys, kind='linear')

            plt.plot(timearr, linv(timearr), markersize=10, marker='+', linewidth=1.5, alpha=1.0, zorder=0, color=ccc, label=k1)
            if wrong_index:
                plt.plot(wrong_times, linv(wrong_times), linewidth=0, marker='*', markersize=9, color='#2ca02c', label='Lin. Interp.')
                self.logger.log("\nInterpolation has been performed!!")
            else:
                self.logger.log("\nAll measurements are reasonable!!")

            self.logger.log(" ",k1)
            self.logger.log("  --> Number of unreasonable tsys values:", len(wrong_index))
            self.logger.log("---------------------------------------------")
            resu = np.asarray([round(g, 1) for g in linv(timearr)])     #% for recording later..
            dali.append(resu); del(resu)
            #% basic statistical properties
            vav = np.round(np.median(linv(timearr)), 1)
            statrl.append(vav)
            vst = np.round(np.std(linv(timearr)), 1)
            statrl.append(vst)
            #% plot layout
            opz1 = np.linspace(plt.ylim()[0], plt.ylim()[1], 7)
            thegap = np.abs(np.abs(opz1[-1])-np.abs(opz1[-2]))
            plt.ylim(opz1[0]-thegap*0.0, opz1[-1]+thegap*1.3)
            if i == 0:
                temp51 = copy.deepcopy(linv(timearr))
                plt.legend(loc=2, fontsize=14, numpoints=1, markerscale=1.5, labelspacing=0.1, handletextpad=0.5, borderaxespad=0.2, handlelength=1.0, handleheight=0.5, ncol=2, columnspacing=1.2)
                plt.xticks([])
                plt.ylabel("Tsys [K]", fontsize=16)
                plt.tick_params(axis='both', which='major', length=5, direction='in', pad=2, color='k')
                plt.tick_params(axis='both', which='minor', length=3, direction='in', pad=2, color='k')
                plt.minorticks_on()
            elif i == 1:
                plt.legend(loc=2, fontsize=14, numpoints=1, markerscale=1.5, labelspacing=0.1, handletextpad=0.5, borderaxespad=0.2, handlelength=1.0, handleheight=0.5, ncol=2, columnspacing=1.2)
                plt.ylabel("Tsys [K]", fontsize=16)
                plt.xlabel("Decimal day", fontsize=16)
                plt.tick_params(axis='both', which='major', length=5, direction='in', pad=2, color='k')
                plt.tick_params(axis='both', which='minor', length=3, direction='in', pad=2, color='k')
                plt.minorticks_on()
        plt.subplots_adjust(left=0.120, bottom=0.107, right=0.940, top=0.930, wspace=0.200, hspace=0.0)
        #% r-median
        p1.text(plt.xlim()[0] + (plt.xlim()[1]-plt.xlim()[0])*0.01, p1.set_ylim()[-1]*1.01, r'R-med=%.1fK' % statrl[0], color='r', fontsize=11)
        #% r-std
        p1.text(plt.xlim()[0] + (plt.xlim()[1]-plt.xlim()[0])*0.24, p1.set_ylim()[-1]*1.01, r'R-std=%.1fK' % statrl[1], color='r', fontsize=11)
        #% l-median
        p1.text(plt.xlim()[0] + (plt.xlim()[1]-plt.xlim()[0])*0.47, p1.set_ylim()[-1]*1.01, r'L-med=%.1fK' % statrl[2], color='m', fontsize=11)
        #% l-std
        p1.text(plt.xlim()[0] + (plt.xlim()[1]-plt.xlim()[0])*0.70, p1.set_ylim()[-1]*1.01, r'L-std=%.1fK' % statrl[3], color='m', fontsize=11)
        #% antcode
        for at in range(2):
            p1.text(plt.xlim()[0] + (plt.xlim()[1]-plt.xlim()[0])*0.93, p1.set_ylim()[-1]*1.01, antcode, color='g', fontsize=18)

        if (hnot != 9999) or (sepif == True):
            plt.savefig(self.logger.pref+"_"+self.antcode+"_"+self.session+"_" + self.freq+ "_"+"Tsys3.png", format='png', dpi=150, transparent=False)
        else:
            plt.savefig(self.logger.pref+"_"+self.antcode+"_"+"_".join(self.session)+"_" + self.freq+ "_"+"Tsys1.png", format='png', dpi=150, transparent=False)
        #% make a new list with the final results
        td,tmh = np.genfromtxt(self.tsys, usecols=(0,1), unpack=True, invalid_raise=False, missing_values='', usemask=False, dtype=str)
        finl = []
        for op in range(len(td)):
            fin99 = str(td[op]) + " " + str(tmh[op]) + " " + str(dali[0][op]) + " " + str(dali[1][op]) + " \n"
            finl.append(fin99)
            del fin99

        self.intp = copy.deepcopy(finl)
        self.intp.append("/\n")
        self.intp.insert(0, self.index)
        self.intp.insert(0, self.gain)
        if (hnot != 9999) or (sepif == True):
            outn = self.logger.pref+"_" +self.antcode+"_"+ "_".join(self.session) + "_" + self.freq+ "_"+"Tsys3.dat"
        else:
            outn = self.logger.pref+"_" +self.antcode+"_"+ "_".join(self.session) + "_" + self.freq+ "_"+"Tsys1.dat"
        with open(outn, 'w') as ft22:
            ft22.writelines(self.intp)
        # plt.show()
        self.logger.log("\n*If the result looks OK with those replaced Tsys values,")
        self.logger.log(" --> Use the output 'xxx_Tsys1.dat'. That's it!")
        self.logger.log("\nOUTPUT: Tsys1.dat file, self.intp, and a png file (i.e., xxxx_Tsys1.png)")

    def interv(a1):     # checking the sampling interval of data.
        a=[]
        for i in np.arange(0, len(a1)):
            b = a1[i+1] - a1[i]
            a = np.append(a, b)
            if i == len(a1)-2:
                print(a)
                break
        return np.mean(a), a

    def smth(self, whatdat=None, siglev=2):
        """
        To smooth out the curve by removing OUTLIERS!
        Para. whatdat: 1 if raw data, 2 if interpolated data.
        Para. siglev: confidence level of the sampling interval for the cutoff; e.g., 2 --> 2*std(the intervals). Default=2.
        output: ~Tsys2.dat, ~Tsys2.png
        """
        self.logger.log("\n@task: smth")
        self.logger.log("===========================")
        antcode = self.antcode
        rlab = "RCP at " + self.freq
        llab = "LCP at " + self.freq
        if whatdat == 1:
            rdat = self.rcp
            ldat = self.lcp
        elif whatdat == 2:
            rdat = np.genfromtxt(self.intp[2:-1], usecols=2)
            ldat = np.genfromtxt(self.intp[2:-1], usecols=3)
        timearr = self.tarr
        self.logger.log("Find a png file of the resultant plot in the current path.")
        fig_pol=plt.figure(figsize=(7,5.5))
        plt.rcParams['legend.frameon'] = 'False'
        plt.rcParams['xtick.labelsize'] = 13
        plt.rcParams['ytick.labelsize'] = 13
        dali = []
        statrl = []
        for i in range(2):
            if i == 0:
                s2 = rdat
                self.logger.log("..\nFirst, RCP !")
                nxt = 1
                ccc = '#1f77b4'
                k1 = rlab
                p1 = plt.subplot(2, 1, nxt)
            elif i == 1:
                s2 = ldat
                self.logger.log("..\nSecond, LCP !")
                nxt = 2
                ccc = '#ff7f0e'
                k1 = llab
                p2 = plt.subplot(2, 1, nxt)

            temp1 = copy.deepcopy(timearr)
            temp2 = copy.deepcopy(s2)
            cutC = siglev * np.std(self.interv(temp2)[1])
            wrong_index = []
            wrong_times = []
            self.logger.log("\n!! Gonna take the Ourliers out from the data !!")
            for aa,bb in enumerate(temp2, 0):
                if aa == len(temp2)-1:
                    break
                if np.abs(temp2[aa+1]-bb) > cutC:
                    if aa < round(len(temp2)/2):
                        half1 = np.mean(temp2[aa:aa+6])  # get a local mean value and judge the data.
                        temarr9 = np.asarray([bb, temp2[aa+1]])
                        inddex = np.abs(temarr9 - half1).argmax()
                        if temarr9[inddex] == bb:
                            if aa in wrong_index:
                                continue
                            else:
                                wrong_index.append(aa)
                                wrong_times.append(temp1[aa])
                                continue
                        else:
                            if aa+1 in wrong_index:
                                continue
                            else:
                                wrong_index.append(aa+1)
                                wrong_times.append(temp1[aa+1])
                                continue
                    else:
                        half2 = np.mean(temp2[aa-5:aa+1])
                        temarr9 = np.asarray([bb, temp2[aa+1]])
                        inddex = np.abs(temarr9 - half2).argmax()
                        if temarr9[inddex] == bb:
                            if aa in wrong_index:
                                continue
                            else:
                                wrong_index.append(aa)
                                wrong_times.append(temp1[aa])
                                continue
                        else:
                            if aa+1 in wrong_index:
                                continue
                            else:
                                wrong_index.append(aa+1)
                                wrong_times.append(temp1[aa+1])
                                continue
            #% x & y array to be used for interpolation
            temp3 = np.asarray([tg for ii, tg in enumerate(temp1) if ii not in wrong_index])
            temp4 = np.asarray([tg for ii, tg in enumerate(temp2) if ii not in wrong_index])
            #=================================
            #% in case the first and/or last measurements are not reasonable!!
            t1 = list(copy.deepcopy(timearr)); t2 = list(copy.deepcopy(temp3))
            z1 = list(copy.deepcopy(s2)); z2 = list(copy.deepcopy(temp4))
            if temp3[0] != timearr[0]:
                self.logger.log("tsys starts with wrong/no value.. to fill in: using the mean of the first 15% of the data; i.e., from left to right")
                fil1 = t1.index(t2[0])
                avg_front = np.mean(temp4[:round(len(t1)*0.15+0.5)])
                for j in range(fil1):
                    z2.insert(0, avg_front)
                    t2.insert(j, t1[j])

            else:
                self.logger.log("Front: good to go (the first data point looks ok)")

            if temp3[-1] != timearr[-1]:
                self.logger.log("tsys ends with wrong/no value.. to fill in: using the mean of the last 15% of the data; i.e., from right to left")
                for k in np.arange(-1, np.negative(len(t1)), -1):
                    if t2[-1] == t1[k]:
                        fil2 = k
                        break

                avg_end = np.mean(temp4[-round(len(t1)*0.15+0.5):])
                for m in np.arange(fil2+1, 0, +1):
                    z2.append(avg_end)
                    t2.append(t1[m])

            else:
                self.logger.log("End: good to go (the last data point looks ok)")
            #=================================
            #% Linear interpolation
            gtg_time = np.asarray(t2)
            gtg_tsys = np.asarray(z2)
            linv = interp1d(gtg_time, gtg_tsys, kind='linear')

            plt.plot(timearr, linv(timearr), markersize=10, marker='+', linewidth=1.5, alpha=1.0, zorder=0, color=ccc, label=k1)
            if wrong_index:
                plt.plot(wrong_times, linv(wrong_times), linewidth=0, marker='*', markersize=9, color='#2ca02c', label='Lin. Interp.')
                self.logger.log("\nInterpolation has been performed!!")
            else:
                self.logger.log("\nAll measurements are reasonable!!")

            self.logger.log(" ",k1)
            self.logger.log("  --> Number of OUTLYING tsys values:", len(wrong_index))
            self.logger.log("---------------------------------------------")
            resu = np.asarray([round(g, 1) for g in linv(timearr)])     #% for recording later..
            dali.append(resu); del(resu)

            #% basic statistical properties
            vav = np.round(np.median(linv(timearr)), 1)
            statrl.append(vav)
            vst = np.round(np.std(linv(timearr)), 1)
            statrl.append(vst)

            #% plot layout
            opz1 = np.linspace(plt.ylim()[0], plt.ylim()[1], 7)
            thegap = np.abs(np.abs(opz1[-1])-np.abs(opz1[-2]))
            plt.ylim(opz1[0]-thegap*0.0, opz1[-1]+thegap*1.3)
            if i == 0:
                temp51 = copy.deepcopy(linv(timearr))
                plt.legend(loc=2, fontsize=14, numpoints=1, markerscale=1.5, labelspacing=0.1, handletextpad=0.5, borderaxespad=0.2, handlelength=1.0, handleheight=0.5, ncol=2, columnspacing=1.2)
                plt.xticks([])
                plt.ylabel("Tsys [K]", fontsize=16)
                plt.tick_params(axis='both', which='major', length=5, direction='in', pad=2, color='k')
                plt.tick_params(axis='both', which='minor', length=3, direction='in', pad=2, color='k')
                plt.minorticks_on()

            elif i == 1:
                plt.legend(loc=2, fontsize=14, numpoints=1, markerscale=1.5, labelspacing=0.1, handletextpad=0.5, borderaxespad=0.2, handlelength=1.0, handleheight=0.5, ncol=2, columnspacing=1.2)
                plt.ylabel("Tsys [K]", fontsize=16)
                plt.xlabel("Decimal day", fontsize=16)
                plt.tick_params(axis='both', which='major', length=5, direction='in', pad=2, color='k')
                plt.tick_params(axis='both', which='minor', length=3, direction='in', pad=2, color='k')
                plt.minorticks_on()

        plt.subplots_adjust(left=0.120, bottom=0.107, right=0.940, top=0.930, wspace=0.200, hspace=0.0)
        #% r-median
        p1.text(plt.xlim()[0] + (plt.xlim()[1]-plt.xlim()[0])*0.01, p1.set_ylim()[-1]*1.01, r'R-med=%.1fK' % statrl[0], color='r', fontsize=11)
        #% r-std
        p1.text(plt.xlim()[0] + (plt.xlim()[1]-plt.xlim()[0])*0.24, p1.set_ylim()[-1]*1.01, r'R-std=%.1fK' % statrl[1], color='r', fontsize=11)
        #% l-median
        p1.text(plt.xlim()[0] + (plt.xlim()[1]-plt.xlim()[0])*0.47, p1.set_ylim()[-1]*1.01, r'L-med=%.1fK' % statrl[2], color='m', fontsize=11)
        #% l-std
        p1.text(plt.xlim()[0] + (plt.xlim()[1]-plt.xlim()[0])*0.70, p1.set_ylim()[-1]*1.01, r'L-std=%.1fK' % statrl[3], color='m', fontsize=11)
        #% antcode
        for at in range(2):
            p1.text(plt.xlim()[0] + (plt.xlim()[1]-plt.xlim()[0])*0.93, p1.set_ylim()[-1]*1.01, antcode, color='g', fontsize=18)

        plt.savefig(self.logger.pref+"_"+self.session+"_"+self.freq+"Tsys2"+"_"+str(siglev)+"-sig"+".png", format='png', dpi=150, transparent=False)
        #% make a new list with the final results
        td,tmh = np.genfromtxt(self.tsys, usecols=(0,1), unpack=True, invalid_raise=False, missing_values='', usemask=False, dtype=str)
        finl = []
        for op in range(len(td)):
            fin99 = str(td[op]) + " " + str(tmh[op]) + " " + str(dali[0][op]) + " " + str(dali[1][op]) + " \n"
            finl.append(fin99)
            del fin99

        self.smip = copy.deepcopy(finl)
        self.smip.append("/\n")
        self.smip.insert(0, self.index)
        self.smip.insert(0, self.gain)
        outn = self.logger.pref +"_"+ self.session + "_" +self.freq+"Tsys2"+"_"+str(siglev)+"-sig"+".dat"
        with open(outn, 'w') as ft22:
            ft22.writelines(self.smip)
        self.logger.log("\n*If the result looks OK with the smoothed out data,")
        self.logger.log(" --> Use the output 'xxx_Tsys2.dat'. That's it!")
        self.logger.log("\nOUTPUT: Tsys2.dat file, self.smip, and a png file (i.e., xxxx_Tsys2.png)")

    def somegains(self):
        """
        To check a recent Gain Info. of some special antennas (i.e., NN, GL, PV, and GB; GBT).
        run it like e.g., self.somegains()
        output: print the Info. in terminal
        """
        self.logger.log("\n@Task: somegains")
        self.logger.log("===========================")
        self.logger.log("The following Gain Info. were checked in 2022.")
        self.logger.log("(..easy to find the Info. for other antennas)")
        self.logger.log("\nGAIN NN ALTAZ DPFU=0.414,0.415 POLY=1.0 /")
        self.logger.log("\nGAIN PV ELEV DPFU=0.141,0.142 POLY=0.91437,0.005815,-0.00014101,1.5319e-6,-7.5242e-9 /")
        self.logger.log("\nGAIN GL ALTAZ DPFU=0.032,0.032 POLY=1.0 /")
        self.logger.log("\nGAIN GB ELEV DPFU = 0.85,0.85 POLY=3.20720839e-01,1.74591876e-03,-1.84532794e-05 /")
        self.logger.log("\nOUTPUT: n/a")






def allinone(logger):
    """
    Collect all to make a single antab file that contains everything.
    First, move all the examined output antab files in one empty directory and do this in there.

    prefx = prefix of the output file name, use the session code (e.g., 2021A --> 'c211')
    frefx = observing frequency of the data (e.g., '3mm')
    OUTPUT: An antab file in a new form - Gain first, and then Tsys measurements.
    """
    logger.log("There should only be the final output data of each antenna in YOUR CURRENT WORKING DIRECTORY!")
    prefx = input("\n<Q> Give a prefix of the output filename (e.g., 'c211ab', 'c222abcd', 'c221abc') \n: ")
    frefx = input("<Q> Array name + Observing wavelength (e.g., '3mmGMVA') \n: ")
    #   sesfx = input("Observing sessions? (e.g., 'abc' for three sessions; ~a,~b,~c.antabs) \n: ")
    # added here: check for .dat files in order to only merges these without logger files
    dat_files = [filename for filename in os.listdir() if filename.endswith('.dat')]

    big1 = []      # Gain
    big2 = []      # Tsys
    arg = []
    logger.log(" ")
    for i in dat_files:
        arg.append(i)
        logger.log(i)
    logger.log("--------------------------------------------\n-> Merge all the above data in one text file") #; time.sleep(0.5)
    # add here:
    # only use file with extension .dat 
    for j in arg:
        with open(j) as v1:
        # v1 = open(j)
            v2 = v1.readlines()
            if v2[-1] == '/':
                v2[-1] = '/\n'
            else:
                pass
            big1 = big1 + v2[:1]
            big2 = big2 + v2[1:]
    big3 = big1 + big2
    outn = "ALLINONE_" + prefx + "_" + frefx + ".antab"
    with open(outn, 'w') as ft54:
        ft54.writelines(big3)



if interactive_mode:
    #  ab hier semi-automatisch
    if sys.argv[1] in ['EF', 'ON', 'MH', 'YS', 'KT', 'KU', 'KY']:
        antcode = sys.argv[1]
        x = TypeB(antcode, interactive_mode, wavelength)
        x.dbcheck()                     #dbcheck includes inpsys if aa==1
        print('----------------------------\n')
        print('The antab generation is done. \nThere are more functions available for further analysis, such as: \n')
        print('> smth: To smooth out the curve by removing OUTLIERS ')
        while True:
            inputB = input("Do You want to continue? y/n: ")
            if inputB == 'y':
                inputB2 = input('Please type in the function shown before with which you want to continue: ')
                if inputB2 == 'smth':
                    x.smth()
            elif inputB == 'n':
                print('Done. Your antab file and corresponding plots for ' + str(antcode)+ ' have been generated. ') 
                break
            else:
                print("Invalid input. Please enter 'y' or 'n'.") 
    elif sys.argv[1] in ['PV', 'NN', 'GL']:
        antcode = sys.argv[1]
        z = TypeC(antcode, interactive_mode,wavelength)
        z.dbcheck()
        print('----------------------------\n')
        print('The antab generation is done. \nThere are more functions available for further analysis, such as: \n')
        print('> smth: To smooth out the curve by removing OUTLIERS \n> somegains: To check a recent Gain Info of some special antennas \n')
        while True:
            inputC = input("Do You want to continue? y/n: ")
            if inputC == 'y':
                inputC2 = input('Please type in the function shown before with which you want to continue: ')
                if inputC2 == 'smth':
                    z.smth()
                if inputC2 == 'somegains':
                    z.somegains()
            elif inputC == 'n':
                print('Done. Your antab file and corresponding plots for ' + str(antcode)+ ' have been generated. ') 
                break
            else:
                print("Invalid input. Please enter 'y' or 'n'.") 
    elif len(sys.argv) ==3 and sys.argv[2].endswith('.key') and sys.argv[1].endswith('.vlba'):
        antabfile = sys.argv[1]
        gainfile = sys.argv[2]
        i=0            
        folder_name = f"VLBA_{i}"
        while folder_name in os.listdir():
            print("The file name exists! Adding a higher number to the end")
            i += 1
            folder_name = f"VLBA_{i}"
        os.makedirs(folder_name, exist_ok=True)
        print("The directory " + folder_name + " has been generated. This folder will contain all .dat and .png files.")
        a = TypeA(antabfile, gainfile, interactive_mode, wavelength) #need antab file and gain file
        a.dbcheck() 
        print('----------------------------\n')
        print('The antab generation is done. \n')
        print("Now all generated files are moved to " + folder_name )
        import shutil
        vlba_stat= ['BR','FD','HN','KP','LA','MK','NL','OV','PT','SC', 'GB']
        #  add all files that contain the the names in vlba_stat in their file name
        search_vlbafile = re.compile('|'.join(vlba_stat))
        for file in os.listdir():
            if (file.endswith('.png') or file.endswith('.dat') or file.endswith('.txt')) and search_vlbafile.search(file):
                shutil.move(file, os.path.join(folder_name, file))
        folder_antab = os.path.join(folder_name, "VLBA_antab")
        os.makedirs(folder_antab, exist_ok=True)        
        print("For merging all .dat files "+ folder_antab + " has been generated")
        
        tsys_files = [f for f in os.listdir(folder_name) if re.search(r'Tsys1\.dat', f)]
        for tsys_file in tsys_files:
            shutil.copy(os.path.join(folder_name, tsys_file), os.path.join(folder_antab, tsys_file))
            
        print("All .dat files are copied to: " + folder_antab)



        # print('The antab generation is done. \nThere are more functions available for further analysis, such as: \n')
        # print('> smth: To smooth out the curve by removing OUTLIERS \n> somegains: To check a recent Gain Info of some special antennas \n> cpif: Copy and paste from a single IF with good Tsys to the IFs with bad Tsys')
        # while True:
        #     inputA = input("\nDo You want to continue? y/n: ")
        #     if inputA == 'y':
        #         inputA2 = input('Please type in the function shown before with which you want to continue: ')
        #         if inputA2 == 'smth':
        #             a.smth()
        #         elif inputA2 == 'somegains':
        #             a.somegains()
        #         elif inputA2 == 'cpif':
        #             a.cpif()
        #     elif inputA == 'n':
        #         print('Done. Your antab file and corresponding plots for VLBA antennas have been generated. ') 
        #         break
        #     else:
        #         print("Invalid input. Please enter 'y' or 'n'.")
        print('----------------------------\n')
        print('Next, you can compress all VLBA antab files to connect all sessions (e.g., c211a, c211b, ..c, ..d): \n')
        inputA = input("Do you want to continue? y/n: ")
        if inputA == 'y':
            a.vlbacon()
        else:
            print('Done. Your antab files and corresponding plots for VLBA antennas have been generated. ')
    elif len(sys.argv) != 2 and not sys.argv[1]:
        print("--------------------------------------------------------------------------")
        print("Invalid Syntax \nFor non-VLBA antennas use: ScriptTsysEdit.py Sationcode \nFor VLBA antennas use: ScriptTsysEdit.py xxx.vlba xxx.key")
        # print("Invalid Syntax for VLBA antennas\nUse: ScriptTsysEdit.py AntabFile GainFile")
        # print("Choose from the following files:")
        for i in os.listdir():
                if i.endswith('.key') or i.endswith('.vlba'):
                    print(i)
else:
    stat_pattern = r'([a-z]{2})\.antab(fs)?$'
    sess_pattern = r'^\w{4}(\w)' 
    station_data = {}
    antfile_dict = {}
    NN_rcp = None
    NN_lcp = None
    for filename in os.listdir():
        if filename.endswith('.antab') or filename.endswith('.antabfs'): 
            match_stat = re.search(stat_pattern, filename)
            statcode = match_stat.group(1).upper()
            session_match = re.search(sess_pattern, filename)
            session_code = session_match.group(1).lower()
            if statcode not in station_data:
                station_data[statcode] = []
            if session_code and session_code not in station_data[statcode]:
                station_data[statcode].append(session_code)
            if statcode not in antfile_dict:
                antfile_dict[statcode] = []
            antfile_dict[statcode].append(filename)   
        elif 'Nn' in filename and filename.endswith('.asc'):
            statcode = 'NN'
            session_code = filename[4].lower()  
            if statcode not in station_data:
                station_data[statcode] = []
            if session_code and session_code not in station_data[statcode]:
                station_data[statcode].append(session_code)
            if statcode not in antfile_dict:
                antfile_dict[statcode] = []
            antfile_dict[statcode].append(filename)
    for statcode, session_codes in station_data.items():
        print("Station code:", statcode)
        print("Session codes:", ", ".join(session_codes), "\n")
        antcode = statcode  
        if antcode in ['EF', 'ON', 'MH', 'YS', 'KT', 'KU', 'KY']:
            antfile = antfile_dict[antcode]
            x = TypeB(antcode, antfile,session_codes, interactive_mode, wavelength)
            x.dbcheck()                     #dbcheck includes inpsys if aa==1
            print("The antab generation for " + statcode + " is done. \n")
        if antcode in ['PV', 'NN','GL']:
            antfile = antfile_dict[antcode]
            print(antfile)
            if antcode == 'NN' and NN_lcp and NN_rcp:
                antfile.extend([NN_rcp, NN_lcp]) 
                x = TypeC(antcode, antfile, session_codes, interactive_mode, wavelength)
                x.dbcheck()  
            x = TypeC(antcode, antfile,session_codes, interactive_mode, wavelength)
            x.dbcheck()                    
            if interactive_mode:
                print('----------------------------\n')
                print('The antab generation is done. \nThere are more functions available for further analysis, such as: \n')
                print('> smth: To smooth out the curve by removing OUTLIERS ')
                while True:
                    inputB = input("Do You want to continue? y/n: ")
                    if inputB == 'y':
                        inputB2 = input('Please type in the function shown before with which you want to continue: ')
                        if inputB2 == 'smth':
                            x.smth()
                    elif inputB == 'n':
                        print('Done. Your antab file and corresponding plots for ' + str(antcode)+ ' have been generated. ') 
                        break
                    else:
                        print("Invalid input. Please enter 'y' or 'n'.") 
            else:
                print("The antab generation for " + statcode + " is done. \n")
    for filename in os.listdir():
        # gain: _gains.key    antabfile: cal.vlba
        if filename.endswith('.key') or filename.endswith('cal.vlba'):
            folder_name = f"VLBA"
            os.makedirs(folder_name, exist_ok=True)
            print("The directory " + folder_name + " has been generated. This folder will contain all .dat and .png files.")
            gainfile = next((f for f in os.listdir() if f.endswith('.key')), None)
            antabfile = next((f for f in os.listdir() if f.endswith('.vlba')), None)

            a = TypeA(antabfile, gainfile, interactive_mode, wavelength) 
            a.dbcheck() 
            print('----------------------------\n')
            print('The antab generation is done. \n')
            print("Now all generated files are moved to " + folder_name )
            import shutil
            vlba_stat= ['BR','FD','HN','KP','LA','MK','NL','OV','PT','SC', 'GB']
            search_vlbafile = re.compile('|'.join(vlba_stat))
            for file in os.listdir():
                if (file.endswith('.png') or file.endswith('.dat') or file.endswith('.txt')) and search_vlbafile.search(file):
                    shutil.move(file, os.path.join(folder_name, file))
            folder_antab = os.path.join(folder_name, "VLBA_antab")
            os.makedirs(folder_antab, exist_ok=True)        
            print("For merging all .dat files "+ folder_antab + " has been generated")
            tsys_files = [f for f in os.listdir(folder_name) if re.search(r'Tsys1\.dat', f)]
            for tsys_file in tsys_files:
                shutil.copy(os.path.join(folder_name, tsys_file), os.path.join(folder_antab, tsys_file))
                
            print("All .dat files are copied to: " + folder_antab)
    print('--------------------------------------------------------')
    print('----------ALL ANTAB FILES HAVE BEEN GENERATED ----------')
    print('--------------------------------------------------------')
    print('--------COPY ALL .DAT FILES INTO SEPARATE FOLDER--------')
    print('--------------------------------------------------------')

    os.makedirs("ALL_ANTAB", exist_ok=True) 
    station_codes = ['EF', 'ON', 'MH', 'YS', 'KT', 'KU', 'KY', 'PV', 'NN', 'GL']
    os.makedirs("ALL_ANTAB", exist_ok=True)
    station_tsys_files = {}
    for filename in os.listdir():
        if filename.endswith('.dat'):
            station_code = filename.split('_')[1]
            if station_code in station_codes:
                if station_code not in station_tsys_files:
                    station_tsys_files[station_code] = []
                station_tsys_files[station_code].append(filename)
    for station_code, tsys_files in station_tsys_files.items():
        tsys0_files = [f for f in tsys_files if '_Tsys0.dat' in f]
        tsys1_files = [f for f in tsys_files if '_Tsys1.dat' in f]
        if tsys1_files:
            for tsys1_file in tsys1_files:
                shutil.copy(tsys1_file, os.path.join("ALL_ANTAB", tsys1_file))
                print(f"Added {tsys1_file} to ALL_ANTAB folder for {station_code}.")
        elif tsys0_files:
            for tsys0_file in tsys0_files:
                shutil.copy(tsys0_file, os.path.join("ALL_ANTAB", tsys0_file))
                print(f"Added {tsys0_file} to ALL_ANTAB folder for {station_code}.")
    print('-------------------COPY VLBA FILES----------------------')
    vlba_antab_folder = "VLBA/VLBA_antab"
    for file in os.listdir(vlba_antab_folder):
        shutil.copy(os.path.join(vlba_antab_folder, file), os.path.join("ALL_ANTAB", file))
        print(f"Addded {file} to ALL_ANTAB.")
    print('--------------------------------------------------------')
    print('--------------MERGE ALL FILES INTO ONE FILE-------------')
    print('--------------------------------------------------------')

    os.chdir("ALL_ANTAB")

    dat_files = [filename for filename in os.listdir() if filename.endswith('.dat')]

    big1 = []      # Gain
    big2 = []      # Tsys
    arg = []
    # logger.log(" ")
    for i in dat_files:
        arg.append(i)
    for j in arg:
        with open(j) as v1:
        # v1 = open(j)
            v2 = v1.readlines()
            if v2[-1] == '/':
                v2[-1] = '/\n'
            else:
                pass
            big1 = big1 + v2[:1]
            big2 = big2 + v2[1:]
    big3 = big1 + big2
    outn = "ALLINONE.antab"
    # outn = "ALLINONE_" + prefx + "_" + frefx + ".antab"
    with open(outn, 'w') as ft54:
        ft54.writelines(big3)
    print(' ')
    print('--------------------------------------------------------')
    print('--------------------------------------------------------')
    print('-------------------------DONE---------------------------')
    print('--------------------------------------------------------')
    print('--------------------------------------------------------')





























