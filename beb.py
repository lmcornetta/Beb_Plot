"""
BEB cross section graphic generator from .txt file

Structure of the .txt file:
+-----------------------------------------------------
MO Binding energy B (eV)  MO Kinetic energy U (a.u.)
          .                         .
          .                         .
          .                         .

+-----------------------------------------------------

+-----------------------------------------------------
Example: 

> cat example.txt
-16.301000  3.340296
-16.286000  3.329908
-15.832000  3.445586
-14.456000  2.043151
-13.880000  2.542708
+-----------------------------------------------------

+-----------------------------------------------------
Run:

> python3 beb.py file.txt
+-----------------------------------------------------

 - Lucas Cornetta JUN/2021

"""

import sys
import numpy as np
import matplotlib.pyplot as plt

def sigma_beb(
    T,          # incident electron energy (float: positive)
    B,          # orbital binding energy (float: negative)
    U,          # orbital kinetic energy (float: positive)
    nocc=2      # occupation of the orbital (int: default = 2)
    ):

    """
    - BEB cross section for a particular MO
    in units of Bohr's radius squared

    """
    Ry = 13.6   # Rydberg constant
    t, u = -T/B, -U/B
    S = 4*np.pi*nocc*(Ry/B)**2

    if t > 1:
        return (S/(t+u+1))*((np.log(t)/2)*(1-1/(t**2))+1-1/t-(np.log(t)/(t+1)))
    
    return 0

class MolData:

    def __init__(
        self,
        filename,
        Tfinal=1e4,
        Npts=5e4
        ):

        """
        - Reading data from file
    
        """
        f = open(filename,"r")
        BU = [[float(line.split(  )[0]), 27.2103*float(line.split(    )[1])] for line in f]
        BUlist = list(zip(*BU))
        f.close()    

        """
        - Setting object attributes (verbatim constructor)

        """
        self.B, self.U = BUlist[0], BUlist[1]   
        self.n = len(self.B)        
        self.T = np.linspace(1.0,Tfinal,int(Npts))
        self.sigma_total = np.zeros(int(Npts))    

    def plot_orb_cs(
        self,
        orb
        ):

        """
        - Plotting BEB CS of the MO 'orb'

        """
        orb_cs = [sigma_beb(T,self.B[orb],self.U[orb]) for T in self.T]
        self.sigma_total += np.array(orb_cs)
        plt.plot(self.T,orb_cs,color='gray',linewidth=1.0,alpha=0.7)

    def plot_cs(
        self
        ):

        """
        - Calling plot_orb_cs(orb) method for all MOs
        - Plotting total BEB CS (sum over all MOs)

        """
        for orb in range(self.n):
            self.plot_orb_cs(orb)
        
        plt.plot(self.T,self.sigma_total,color='black',linewidth=2.5,zorder=self.n+1)
        plt.show()
        
if __name__=="__main__":
    
    plt.rcParams['legend.fontsize'] = 14
    plt.rcParams['xtick.labelsize'] = 13
    plt.rcParams['ytick.labelsize'] = 13
    plt.rcParams["legend.frameon"] = False
    plt.xscale("log")
    plt.xlabel("Incident energy (eV)",fontsize=13)
    plt.ylabel("Ionization cross section $(a_0^2)$",fontsize=13)    

    mol = MolData(sys.argv[1])
    mol.plot_cs()
