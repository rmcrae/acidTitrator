#intention: general n-protic acid titrator program:
#update Fall 2014

from pylab import *

class gen_acid:
    def __init__(self,K=[1.8e-5],ca0=0.100,V0=50.0):
        self.Kw = 1.0e-14
        self.K = K
        self.ca0 = ca0
        self.V0 = V0
        self.na0 = ca0 * V0
        self.nH = len(K)

    def titrate(self,cOH=0.10,npoints=200):
        self.cOH = cOH
        self.Veq = self.na0 / self.cOH
        delv = self.Veq / float(npoints)
        self.vlist = arange(0,(self.nH+0.5)*self.Veq,delv)
        self.Hlist = []
        for v in self.vlist:
            nOH = self.cOH * v
            V = v + self.V0
            cNa = nOH/V
            T = self.na0 / V
            coeff = [1,self.K[0] + cNa,self.K[0]*(cNa - self.Kw/self.K[0] - T)]
            iprod = self.K[0]
            for i in range(self.nH - 1):
                iprod = iprod * self.K[i+1]
                coeff[2+i] = coeff[2+i] + iprod
                coeff.append(iprod * (cNa - self.Kw/self.K[i+1] - (i+2)*T))

            coeff.append(-iprod*self.Kw)
            p = poly1d(coeff)
            H = max(p.r.real)
            self.Hlist.append(H)

        self.pH = -log10(array(self.Hlist))

if __name__ == "__main__":
    acid = gen_acid([6.31e-3,5.62e-5,2.14e-10])
    acid.titrate()
    plot(acid.vlist, acid.pH)
    title("Phosphoric Acid")
    xlabel("vol OH- added")
    ylabel("pH")
    show()
                     
            
            
